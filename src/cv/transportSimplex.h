/* RCSid $Id: transportSimplex.h,v 2.2 2015/08/01 23:27:04 greg Exp $ */
/*
transportSimplex.h

A C++ implementation of the transportation simplex algorithm. 

Last edit Sept 3 2006

Copyright (C) 2006 Darren MacDonald
dmacdonald@discover.uottawa.ca
www.site.uottawa.ca/~dmacd070/emd

except some data structures and interface adapted from 
http://ai.stanford.edu/~rubner/emd/emd.c, 
Copyright 1998 Yossi Rubner.

This algorithm solves the problem of finding the least-cost way to transport
goods from a set of source (supply) nodes to a set of sink (demand) nodes.

To use the code, simply include this file in the user code and add 'using namespace t_simplex;'. 
Organize the nodes into signatures using the TsSignature<TF> class, which contains a feature 
array (of type TF) and an array of the respective weights of the features.

Solve the system by calling
transportSimplex(&Sig1, &Sig2, grndDist);
where Sig1 is the source signature, containing an array of source features and an array of their respective
supplies, Sig2 is the sink signature, containing the sink features and and their respective demands, and 
grndDist is a pointer to a function which accepts two features pointers as arguments and returns the unit
cost of transporting goods between the two.

transportSimplex either returns the optimal transportation cost or throws a TsError.

For more information see the documentation at 
www.site.uottawa.ca/~dmacd070/emd/index.html 
and the example implemenation at
www.site.uottawa.ca/~dmacd070/emd/main.cpp.
*/

#ifndef _T_SIMPLEX_H
#define _T_SIMPLEX_H
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <new>

#define TSINFINITY       1e20
#define TSEPSILON        1e-6
#define TSPIVOTLIMIT     0.00

namespace t_simplex {

/* DECLARATION OF DATA TYPES */
enum TsError { TsErrBadInput };

//TsSignature is used for inputting the source and sink signatures
template <class TF>
class TsSignature {
public:
	int n;									// Number of features in the signature 
	TF *features;							// Pointer to the features vector 
	double *weights;						// Pointer to the weights of the features 
	TsSignature(int nin, TF *fin, double * win):n(nin), features(fin), weights(win){};
};

//TsFlow is used for outputting the final flow table
typedef struct TsFlow {
	int from;								// Feature number in signature 1 
	int to;									// Feature number in signature 2 
	double amount;							// Amount of flow from signature1.features[from] to signature2.features[to]
} TsFlow;

// TsBasic is used for 2D lists, allowing for easy navigation of the basic variables 
typedef struct TsBasic{
	int i, j;
	double val;
	TsBasic *nextSnk, *prevSnk;				// next/previous node in the column
	TsBasic *nextSrc, *prevSrc;				// next/previous node in the row
} TsBasic;

// TsStone is used for _BFS
typedef struct TsStone {
	struct TsStone *prev;
	struct TsBasic *node;
} TsStone;

// TsRusPen is used for 1D lists in _initRussel
typedef struct TsRusPen {
  int i;
  double val;
  struct TsRusPen *next, *prev;
} TsRusPen;

// TsVogPen is used for 1D lists in _initVogel
typedef struct TsVogPen {
  int i;
  struct TsVogPen *next, *prev;
  int one, two;
  double oneCost, twoCost;
} TsVogPen;

/* DECLARATION OF GLOBALS */
double ** _tsC = NULL;				// Cost matrix
double _tsMaxC;						// Maximum of all costs
double _tsMaxW;						// Maximum of all weights


/* INTERNAL FUNCTIONS */
double _pivot(TsBasic * basics, TsBasic ** srcBasics, TsBasic ** snkBasics, bool ** isBasic, int n1, int n2);
TsStone * _BFS(TsStone *stoneTree, TsBasic ** srcBasics, TsBasic ** snkBasics, bool complete = false);
void _initVogel(double *S, double *D, TsBasic * basicsEnd, TsBasic ** srcBasics, TsBasic ** snkBasics, bool ** isBasic, int n1, int n2);
void _initRussel(double *S, double *D, TsBasic * basicsEnd, TsBasic ** srcBasics, TsBasic ** snkBasics, bool ** isBasic, int n1, int n2);

/*
transportSimplex() - Program entry point. 

signature1 and signature2 define the source and sink sets. Unit costs between two features are computed from 
grndDist. Signature weights must be positive and all costs must be positive.
TsFlow is an output parameter which can be set to an array that will be filled with the final flow amounts. 
The array must be of size signature1->n + signature2->n - 1 . flowSize is a pointer to an integer which 
indicates the number of functional entries in Flow, because all spaces are not necessarily used. Flow and 
flowSize can be set to NULL or omitted from the argument list if this information is not important. 

The return value is the transportation cost. 

*/
template <class TF>
double transportSimplex(TsSignature<TF> *signature1, TsSignature<TF> *signature2, double (*grndDist)(TF *, TF *),
						TsFlow *flowTable = NULL, int *flowSize = NULL) {

	int n1, n2;								// number of features in signature1 and Signature 2
	int i, j;
	
	double totalCost;
	double w;
  
	double srcSum, snkSum, diff;  
		
	TF *P1, *P2;
	
	n1 = signature1->n;
	n2 = signature2->n;
	
	TsBasic * basics = NULL;					///Array of basic variables. 
	bool ** isBasic = NULL;						//Flag matrix. isBasic[i][j] is true there is flow between source i and sink j
	TsBasic **srcBasics = NULL;					//Array of pointers to the first basic variable in each row
	TsBasic **snkBasics = NULL;					//Array of pointers to the first basic variable in each column
	double * src = NULL;						//Array of source supplies
	double * snk =NULL;							//Array of sink demands
	
	
	// Equalize source and sink weights. A dummy source or sink may be added to equalize the total sink
	// and source weights. n1 = signature1->n + 1 if there is a dummy source, and n2 = signature2->n + 1
	// if there is a dummy sink.
	srcSum = 0.0;
	for(i=0; i < n1; i++)
		srcSum += signature1->weights[i];

	snkSum = 0.0;
	for(j=0; j < n2; j++)
		snkSum += signature2->weights[j];
	
	diff = srcSum - snkSum;
	if (fabs(diff) > TSEPSILON * srcSum) {
		if (diff < 0.0)
			n1++;
		else 
			n2++;
	}
	_tsMaxW = srcSum > snkSum ? srcSum : snkSum;
	w = srcSum < snkSum ? srcSum : snkSum;
    
	//Allocate memory for arrays
	try {
		
		basics = new TsBasic[n1 + n2];
	
		isBasic = new bool*[n1];
		for(i = 0; i < n1; i++) 
			isBasic[i] = NULL;

		for(i = 0; i < n1; i++) {
			isBasic[i] = new bool[n2];
			for (j=0; j < n2; j++)
				isBasic[i][j] = 0;
		}

		srcBasics = new TsBasic*[n1];
		for (i = 0; i < n1; i++)
			srcBasics[i] = NULL;
	
		snkBasics = new TsBasic*[n2];
		for (i = 0; i < n2; i++)
			snkBasics[i] = NULL;
		
		// Compute the cost matrix
		_tsMaxC = 0;
		_tsC = new double*[n1];
		for(i=0; i < n1; i++) 
			_tsC[i] = NULL;

		for(i=0, P1=signature1->features; i < n1; i++, P1++) {
			_tsC[i] = new double[n2];					//What happens if bad  here?
			for(j=0, P2=signature2->features; j < n2; j++, P2++) {
				if (i == signature1->n || j == signature2->n) {
					_tsC[i][j] = 0;			// cost is zero for flow to or from a dummy
				} else {
	         		_tsC[i][j] = grndDist(P1, P2);
					if (_tsC[i][j] < 0) throw TsErrBadInput;
				}

				
				if (_tsC[i][j] > _tsMaxC)
					_tsMaxC = _tsC[i][j];
			}
		}
		
		src = new double[n1];					//init the source array
		for (i = 0; i < signature1->n; i++) {
			src[i] = signature1->weights[i];
			if (src[i] < 0) throw TsErrBadInput;
		}
		if (n1 != signature1->n)
			src[signature1->n] = -diff;
	
		snk = new double[n2];					//init the sink array
		for (i = 0; i < signature2->n; i++) {
			snk[i] = signature2->weights[i];
			if (snk[i] < 0) throw TsErrBadInput;
		}
		if (n2 != signature2->n)
			snk[signature2->n] = diff;
			
		// Find the initail basic feasible solution. Use either _initRussel or _initVogel
		
		
		_initRussel(src, snk, basics, srcBasics, snkBasics, isBasic, n1, n2);
		//_initVogel(src, snk, basics, srcBasics, snkBasics, isBasic, n1, n2);
		
			
		// Enter the main pivot loop
		totalCost = _pivot(basics, srcBasics, snkBasics, isBasic, n1, n2);
		
	} catch (...) {
		for(i = 0; i < n1; i++)
			delete[] isBasic[i];
		delete[] isBasic;

		for (i = 0; i < n1; i++)
			delete[] _tsC[i];
		delete[] _tsC;
	
		delete[] src;
		delete[] snk;
		delete[] srcBasics;
		delete[] snkBasics;
		delete[] basics;
 
		throw;
	}

	// Fill the Flow data structure
	TsBasic * basicPtr;
	TsFlow * flowPtr = flowTable;
	if (flowTable != NULL) {
		for (i = 0; i < n1+n2 -1; i++) {
			basicPtr = basics + i;
			if (isBasic[basicPtr->i][basicPtr->j] && basicPtr->i != signature1->n && basicPtr->j != signature2->n && basicPtr->val != 0.0) {
				flowPtr->to = basicPtr->j;
				flowPtr->from = basicPtr->i;
				flowPtr->amount = basicPtr->val;
				flowPtr++;
			}
		}
	}
	if (flowSize != NULL) {
		*flowSize = (int)(flowPtr - flowTable);
	}
	
	
	for(i = 0; i < n1; i++)
		delete[] isBasic[i];
	delete[] isBasic;

	for (i = 0; i < n1; i++)
		delete[] _tsC[i];
	delete[] _tsC;
	
	delete[] src;
	delete[] snk;
	delete[] srcBasics;
	delete[] snkBasics;
	delete[] basics;
 
	return totalCost;
}

/*
Main pivot loop. 
Pivots until the system is optimal and return the optimal transportation cost. 
*/
double _pivot(TsBasic * basics, TsBasic ** srcBasics, TsBasic ** snkBasics, bool ** isBasic, int n1, int n2)  {
	
	double * srcDuals = NULL;
	double * snkDuals = NULL;
	TsStone * stonePath = NULL;
	
	TsStone * spitra, * spitrb, * leaving;
	TsBasic * XP;
	TsBasic * basicsEnd = basics + n1 + n2;
	TsBasic * entering = basicsEnd - 1 ;
	TsBasic dummyBasic;
		dummyBasic.i = -1;
		dummyBasic.j = 0;

	int i,j, lowI, lowJ;
	double objectiveValue = TSINFINITY, oldObjectiveValue = 0;
	double lowVal;
	int numPivots = 0;

	try {
		srcDuals = new double[n1];
		snkDuals = new double[n2];
		stonePath = new TsStone[n1 + n2];
	} catch (std::bad_alloc) {
		delete[] srcDuals;
		delete[] snkDuals;
		delete[] stonePath;
		throw;
	}

	while (1) {
		
		oldObjectiveValue = objectiveValue;
		objectiveValue = 0;
		
		for (XP = basics; XP != basicsEnd; XP++){
			if (XP != entering) 
				objectiveValue += _tsC[XP->i][XP->j] * XP->val;
		}

        // Compute the dual variables for each row and column. Begin by finding a spanning tree (stonepath)
		// of the basic variables using the breadth-first search routine seeded at an imaginary basic 
		// variable in the first column. The dual variables can then be computed incrementally by traversing 
		// the tree.
		stonePath[0].node = &dummyBasic;
		stonePath[0].prev = NULL;
		spitrb = _BFS(stonePath, srcBasics, snkBasics, true);

		spitra = stonePath;
		snkDuals[spitra->node->j] = 0;
		for (spitra++; spitra != spitrb; spitra++) {
			if (spitra->node->i == spitra->prev->node->i) {
				//node is in same row as parent
				snkDuals[spitra->node->j] = _tsC[spitra->node->i][spitra->node->j] - srcDuals[spitra->node->i];
			}  else if (spitra->node->j == spitra->prev->node->j) {
				srcDuals[spitra->node->i] = _tsC[spitra->node->i][spitra->node->j] - snkDuals[spitra->node->j];
			}
		}
		
		// After computing the duals, find the non-basic variable that has the greatest negative value of
		// delta = _tsC[i][j] - srcDuals[i] - snkDuals[j]. This is the entering variable
		lowVal = 0.0;
		for (i = 0; i < n1; i++) 
			for (j = 0; j < n2; j++) 
				if (!isBasic[i][j] &&  _tsC[i][j] - srcDuals[i] - snkDuals[j] < lowVal) {
						lowVal = _tsC[i][j] - srcDuals[i] - snkDuals[j];
						lowI = i;
						lowJ = j;
					}
		
		// If all delta values are non-negative, the table is optimal
		if (lowVal >= -TSEPSILON * _tsMaxC || (oldObjectiveValue - objectiveValue) < TSPIVOTLIMIT) {
			delete[] srcDuals;
			delete[] snkDuals;
			delete[] stonePath;
			//std::cout << numPivots << "\t";
			return objectiveValue;
		}	
		
		// Add the entering variable
		entering->i = lowI;
		entering->j = lowJ;
		isBasic[lowI][lowJ] = 1;
		entering->val = 0;

		entering->nextSrc = srcBasics[lowI];
		if (srcBasics[lowI] != NULL) srcBasics[lowI]->prevSrc = entering;
		entering->nextSnk = snkBasics[lowJ];
		if (snkBasics[lowJ] != NULL) snkBasics[lowJ]->prevSnk = entering;
		
		srcBasics[lowI] = entering;
		entering->prevSrc = srcBasics[lowI];
		snkBasics[lowJ] = entering;
		entering->prevSnk = snkBasics[lowJ];
		
		stonePath[0].node = entering;
		stonePath[0].prev = NULL;

		// Use breadth-first search to find a loop of basics.
		spitra = spitrb = _BFS(stonePath, srcBasics, snkBasics);
		lowVal = TSINFINITY;
		bool add = false;

		// Find the lowest flow along the loop (leaving variable)
		do {
			if (!add && spitrb->node->val < lowVal) {
				leaving = spitrb;
				lowVal = spitrb->node->val;
			}
			add = !add;
		} while ((spitrb = spitrb->prev));
		
		add = false;
		spitrb = spitra;

		// Alternately increase and decrease flow along the loop
		do {
			if (add) spitrb->node->val += lowVal;
			else spitrb->node->val -= lowVal;
			add = !add;
		} while ((spitrb = spitrb->prev));

		i = leaving->node->i;
		j = leaving->node->j;
		isBasic[i][j] = 0;
		

		if (srcBasics[i] == leaving->node) {
			srcBasics[i] = leaving->node->nextSrc;
			srcBasics[i]->prevSrc = NULL;
		} else {
			leaving->node->prevSrc->nextSrc = leaving->node->nextSrc;
			if (leaving->node->nextSrc != NULL)
				leaving->node->nextSrc->prevSrc = leaving->node->prevSrc;
		}
		
		if (snkBasics[j] == leaving->node) {
			snkBasics[j] = leaving->node->nextSnk;
			snkBasics[j]->prevSnk = NULL;
		} else {
			leaving->node->prevSnk->nextSnk = leaving->node->nextSnk;
			if (leaving->node->nextSnk != NULL)
				leaving->node->nextSnk->prevSnk = leaving->node->prevSnk;
		}
		entering = leaving->node;
		numPivots++;
	}
};


/*******************
		_BFS
Perform a breadth-first search of the basic varaibles. The search tree is stored in stoneTree, where each
'stone' contains a pointer to a basic variable, and a pointer to that stone's parent. srcBasics and snkBasics
are arrays of linked lists which allow for easy identification of a node's neighbours in the flow network.
If complete == true, find all basics in the table. Otherwise, terminate when a basic variable which completes
a loop has been found and return a pointer to the final stone in the loop.
*********************/
TsStone * _BFS(TsStone * stoneTree, TsBasic ** srcBasics, TsBasic ** snkBasics, bool complete) {
	bool column = true;
	int jumpoffset = 0;
	TsBasic * bitr;
	TsStone * sitra = &stoneTree[0], * sitrb = &stoneTree[1];
	do {
		if (column) {
			for (bitr = snkBasics[sitra->node->j]; bitr != NULL; bitr = bitr->nextSnk) {
				if (bitr != sitra->node){
					sitrb->node = bitr;
					sitrb->prev = sitra;
					sitrb++;				
				}
			}
		} else {
			for (bitr = srcBasics[sitra->node->i]; bitr != NULL; bitr = bitr->nextSrc) {
				if (bitr != sitra->node){
					sitrb->node = bitr;
					sitrb->prev = sitra;
					sitrb++;				
				}
			}
		}
	
		sitra++;	
		if (sitra == sitrb) //no cycle found and no cycles in tree
			return sitra;

		if (sitra->node->i == sitra->prev->node->i) 
			column = true;
		else 
			column = false;

			// cycle found
		if (!complete && sitra->node->i == stoneTree[0].node->i && sitra->node->j != stoneTree[0].node->j && column == false) 
			return sitra;
	} while(1);	
}


// Helper function for _initVogel
inline void addPenalty(TsVogPen * pitr, double cost, int i) {
	if (cost < pitr->oneCost) {
		pitr->twoCost = pitr->oneCost;
		pitr->two = pitr->one;
		pitr->oneCost = cost;
		pitr->one = i;
	} else if (cost < pitr->twoCost) {
		pitr->twoCost = cost;
		pitr->two = i;
	}
}


/**********************
    Vogel's initialization method
**********************/
void _initVogel(double *S, double *D, TsBasic * basicsEnd, TsBasic ** srcBasics, TsBasic ** snkBasics, bool ** isBasic, int n1, int n2) {
	int i, j;
	TsVogPen *srcPens = NULL;
	TsVogPen *snkPens = NULL;
	TsVogPen *pitra, *pitrb;  //iterators
	TsVogPen *maxPen;
	TsVogPen srcPenHead, snkPenHead;
	bool maxIsSrc;
	double lowVal;

	try {
		srcPens = new TsVogPen[n1];
		snkPens = new TsVogPen[n2];
	} catch (std::bad_alloc) {
		delete[] srcPens;
		delete[] snkPens;
		throw;	
	}

	srcPenHead.next = pitra = srcPens;
	for (i=0; i < n1; i++) {
		pitra->i = i;
		pitra->next = pitra+1;
		pitra->prev = pitra-1;
		pitra->one = pitra->two = 0;
		pitra->oneCost = pitra->twoCost = TSINFINITY;
		pitra++;	  
	}
	(--pitra)->next = NULL;
	srcPens[0].prev = &srcPenHead;

	snkPenHead.next = pitra = snkPens;
	for (i=0; i < n2; i++) {
		pitra->i = i;
		pitra->next = pitra+1;
		pitra->prev = pitra-1;
		pitra->one = pitra->two = 0;
		pitra->oneCost = pitra->twoCost = TSINFINITY;
		pitra++;	  
	}
	(--pitra)->next = NULL;
	snkPens[0].prev = &snkPenHead;


	for (pitra = srcPenHead.next, i=0; pitra != NULL; pitra = pitra->next, i++)
		for (pitrb = snkPenHead.next, j=0; pitrb != NULL; pitrb = pitrb->next, j++) {
		//initialize Source Penalties;
			addPenalty(pitra, _tsC[i][j], j);
			addPenalty(pitrb, _tsC[i][j], i);
	  }
  

	while (srcPenHead.next != NULL && snkPenHead.next != NULL) {
		maxIsSrc = true;
		for (maxPen = pitra = srcPenHead.next; pitra != NULL; pitra = pitra->next) 
			if ((pitra->twoCost - pitra->oneCost) > (maxPen->twoCost - maxPen->oneCost))
                maxPen = pitra;		
		
		for (pitra = snkPenHead.next; pitra != NULL; pitra = pitra->next) 
			if ((pitra->twoCost - pitra->oneCost) > (maxPen->twoCost - maxPen->oneCost)) {
				maxPen = pitra;
				maxIsSrc = false;
			}
		
		if (maxIsSrc) {
			i = maxPen->i;
			j = maxPen->one;
		} else {
			j = maxPen->i;
			i = maxPen->one;		
		}

		if (D[j] - S[i] > _tsMaxW * TSEPSILON || (srcPenHead.next->next != NULL && fabs(S[i] - D[j]) < _tsMaxW * TSEPSILON)) {
			//delete source
			lowVal = S[i];
			maxPen = srcPens + i;
			maxPen->prev->next = maxPen->next;
			if (maxPen->next != NULL)
				maxPen->next->prev = maxPen->prev;

			for (pitra = snkPenHead.next; pitra != NULL; pitra = pitra->next) {
				if (pitra->one == i || pitra->two == i){
					pitra->oneCost = TSINFINITY;
					pitra->twoCost = TSINFINITY;
					for (pitrb = srcPenHead.next; pitrb != NULL; pitrb = pitrb->next)
						addPenalty(pitra, _tsC[pitrb->i][pitra->i], pitrb->i);
				}
			}
		} else {
			//delete sink
			lowVal = D[j];
			maxPen = snkPens + j;
			maxPen->prev->next = maxPen->next;
			if (maxPen->next != NULL)
				maxPen->next->prev = maxPen->prev;

			for (pitra = srcPenHead.next; pitra != NULL; pitra = pitra->next) {
				if (pitra->one == j || pitra->two == j){
					pitra->oneCost = TSINFINITY;
					pitra->twoCost = TSINFINITY;
					for (pitrb = snkPenHead.next; pitrb != NULL; pitrb = pitrb->next)
						addPenalty(pitra, _tsC[pitra->i][pitrb->i], pitrb->i);
				}
			}
		}

		S[i] -= lowVal;
		D[j] -= lowVal;

		isBasic[i][j] = 1; 
		basicsEnd->val = lowVal;
		basicsEnd->i = i;
		basicsEnd->j = j;
	
		basicsEnd->nextSnk = snkBasics[j];
		if (snkBasics[j] != NULL) snkBasics[j]->prevSnk = basicsEnd;
		basicsEnd->nextSrc = srcBasics[i];
		if (srcBasics[i] != NULL) srcBasics[i]->prevSrc = basicsEnd;

		srcBasics[i] = basicsEnd;
		basicsEnd->prevSnk = NULL;
		snkBasics[j] = basicsEnd;
		basicsEnd->prevSrc = NULL;

		basicsEnd++;
		
	}
	delete[] srcPens;
	delete[] snkPens;
}


/**********************
    Russel's initialization method
**********************/
void _initRussel(double *S, double *D, TsBasic * basicsEnd, TsBasic ** srcBasics, TsBasic ** snkBasics, bool ** isBasic, int n1, int n2) {
	double ** Delta = NULL;
	int i, j, lowI, lowJ;
	TsRusPen *U  = NULL;
	TsRusPen *V = NULL;
	TsRusPen *Uhead, *Vhead;	
	TsRusPen *Uitr, *Vitr;
	double cost, lowVal;
	
	try {
		Delta = new double*[n1];
		for (i = 0; i < n1; i++)
			Delta[i] = new double[n2];

		U = new TsRusPen[n1];
		V = new TsRusPen[n2];
	} catch (std::bad_alloc) {
		for (i = 0; i < n1; i++)
			delete[] Delta[i];
		delete[] Delta;

		delete[] U;
		delete[] V;
		throw;
	}
	
	for (i = 0; i < n1; i++) {
		U[i].i = i;
		U[i].val = 0;
		U[i].next = &U[i+1];
		U[i].prev = &U[i-1];
	}
	U[n1-1].next = NULL;
	U[0].prev = NULL;
	
	for (i = 0; i < n2; i++) {
		V[i].i = i;
		V[i].val = 0;
		V[i].next = &V[i+1];
		V[i].prev = &V[i-1];
	}
	V[n2-1].next = NULL;
	V[0].prev = NULL;

	for (i = 0; i < n1; i++)
		for (j = 0; j < n2; j++) {
			cost = _tsC[i][j];
			if (cost > U[i].val)
				U[i].val = cost;
			if (cost > V[j].val)
				V[j].val = cost;
		}

	for (i = 0; i < n1; i++) 
		for (j = 0; j < n2; j++) 
			Delta[i][j] = _tsC[i][j] - U[i].val - V[j].val;
	
	Uhead = U;
	Vhead = V;
	while (Uhead != NULL && Vhead != NULL) {
		
		//Find lowest Delta
		lowVal = TSINFINITY;
		for (Uitr = Uhead; Uitr != NULL; Uitr = Uitr->next) 
			for (Vitr = Vhead; Vitr != NULL; Vitr = Vitr->next) 
				if (Delta[Uitr->i][Vitr->i] < lowVal) {
					lowI = Uitr->i;
					lowJ = Vitr->i;
					lowVal = Delta[Uitr->i][Vitr->i];
				}

		
		if (D[lowJ] - S[lowI] > _tsMaxW * TSEPSILON || (fabs(S[lowI] - D[lowJ]) < _tsMaxW * TSEPSILON && Uhead->next != NULL)) {
			//Delete Source
			if (&U[lowI] == Uhead) {	//Entering variable is first in list
				Uhead = Uhead->next;
				if (Uhead != NULL) Uhead->prev = NULL;
			} else {
				U[lowI].prev->next = U[lowI].next;	//Entering variable is in middle of  list;
				if (U[lowI].next != NULL)			//Entering variable is at the end of the list;
					U[lowI].next->prev = U[lowI].prev;
			}
			//See if this source was the maximum cost for any dest
			for (Vitr = Vhead; Vitr != NULL; Vitr = Vitr->next)  {
				if (Vitr->val == _tsC[lowI][Vitr->i]) {
					//it is; update the dest
					//find maximum cost in the dest
					Vitr->val = 0;
					for (Uitr = Uhead; Uitr != NULL; Uitr = Uitr->next)
						if (_tsC[Uitr->i][Vitr->i] > Vitr->val)
							Vitr->val = _tsC[Uitr->i][Vitr->i];
					//update Delta
					for (Uitr = Uhead; Uitr != NULL; Uitr = Uitr->next)
						Delta[Uitr->i][Vitr->i] = _tsC[Uitr->i][Vitr->i] - Uitr->val - Vitr->val;

				}
			}
			lowVal = S[lowI];
			
		} else {
			//Delete Dest
			if (&V[lowJ] == Vhead) {	//Entering variable is first in list
				Vhead = Vhead->next;
				if (Vhead != NULL) Vhead->prev = NULL;
			} else {
				V[lowJ].prev->next = V[lowJ].next;	//Entering variable is in middle of  list;
				if (V[lowJ].next != NULL)			//Entering variable is at the end of the list;
					V[lowJ].next->prev = V[lowJ].prev;
			}
			//See if this source was the maximum cost for any dest
			for (Uitr = Uhead; Uitr != NULL; Uitr = Uitr->next)  {
				if (Uitr->val == _tsC[Uitr->i][lowJ]) {
					//it is; update the dest
					//find maximum cost in the dest
					Uitr->val = 0;
					for (Vitr = Vhead; Vitr != NULL; Vitr = Vitr->next)
						if (_tsC[Uitr->i][Vitr->i] > Uitr->val)
							Uitr->val = _tsC[Uitr->i][Vitr->i];
					//update Delta
					for (Vitr = Vhead; Vitr != NULL; Vitr = Vitr->next)
						Delta[Uitr->i][Vitr->i] = _tsC[Uitr->i][Vitr->i] - Uitr->val - Vitr->val;

				}
			}
			lowVal = D[lowJ];
		}

		S[lowI] -= lowVal;
		D[lowJ] -= lowVal;

		isBasic[lowI][lowJ] = 1; 
		basicsEnd->val = lowVal;
		basicsEnd->i = lowI;
		basicsEnd->j = lowJ;
	
		basicsEnd->nextSnk = snkBasics[lowJ];
		if (snkBasics[lowJ] != NULL) snkBasics[lowJ]->prevSnk = basicsEnd;
		basicsEnd->nextSrc = srcBasics[lowI];
		if (srcBasics[lowI] != NULL) srcBasics[lowI]->prevSrc = basicsEnd;

		srcBasics[lowI] = basicsEnd;
		basicsEnd->prevSnk = NULL;
		snkBasics[lowJ] = basicsEnd;
		basicsEnd->prevSrc = NULL;

		basicsEnd++;
		
	}
	
	delete[] U;
	delete[] V;
	for (i = 0; i < n1; i++)
		delete[] Delta[i];
	delete[] Delta;
}
}
#endif

