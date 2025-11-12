#ifndef lint
static const char RCSid[] = "$Id$";
#endif
/*
 *  RcontribSimulManager.cpp
 *
 *	Rcontrib simulation manager implementation
 *
 *  Created by Greg Ward on 10/17/2024.
 */

#include <unistd.h>
#include <ctype.h>
#include "platform.h"
#include "selcall.h"
#include "RcontribSimulManager.h"
#include "func.h"
#include "resolu.h"
#include "source.h"

char			RCCONTEXT[] = "RC.";

extern const char	HDRSTR[];
extern const char	BIGEND[];
extern const char	FMTSTR[];

// new/exclusive, overwrite if exists, or recover data
int	RSDOflags[] = {RDSwrite|RDSexcl|RDSextend, RDSwrite|RDSextend,
				RDSread|RDSwrite};

static const char	ROWZEROSTR[] = "NROWS=0000000000000000\n";
#define LNROWSTR	6

// Modifier channel for recording contributions (no constructor/destructor)
struct RcontribMod {
	RcontribOutput *	opl;		// pointer to first output channel
	char *			params;		// parameters string
	EPNODE *		binv;		// bin expression (NULL if 1 bin)
	int			nbins;		// bin count this modifier
	int			coffset;	// column offset in bytes
	DCOLORV			cbin[1];	// bin accumulator (extends struct)
				/// Access specific (spectral) color bin
	DCOLORV *		operator[](int n) {
					if ((n < 0) | (n >= nbins)) return NULL;
					return cbin + n*NCSAMP;
				}
};

// Struct used to assign record calculation to child
struct RowAssignment {
	uint32			row;		// row to do
	uint32			ac;		// accumulation count
};

static const char	ROW_DONE[] = "ROW FINISHED\n";

// allocate rcontrib accumulator
static RcontribMod *
NewRcMod(const char *prms, const char *binexpr, int ncbins)
{
	if (binexpr && !*binexpr) binexpr = NULL;
	if (!prms) prms = "";
	if ((ncbins > 1) & !binexpr) {
		error(USER, "missing bin expression");
		return NULL;
	}
	if (ncbins < 1) ncbins = 1;
	
	RcontribMod *	mp = (RcontribMod *)ecalloc(1, sizeof(RcontribMod) +
						sizeof(DCOLORV)*(NCSAMP*ncbins-1) +
						strlen(prms)+1);

	if (binexpr) {				// get/check bin expression
		mp->binv = eparse(const_cast<char *>(binexpr));
		if (mp->binv->type == NUM) {	// constant expression (0)?
			if ((int)evalue(mp->binv) != 0) {
				sprintf(errmsg, "illegal non-zero constant for bin (%s)",
						binexpr);
				error(USER, errmsg);
			}
			if (ncbins > 1) {
				sprintf(errmsg, "bad bin count (%d should be 1)", ncbins);
				error(USER, errmsg);
			}
			epfree(mp->binv, true);
			mp->binv = NULL;
			prms = "";
			ncbins = 1;
		}
	}
	mp->params = strcpy((char *)(mp->cbin + ncbins*NCSAMP), prms);
	mp->nbins = ncbins;
	return mp;
}

// Free an RcontribMod (public for RcontribSimulManager constructor)
void
FreeRcMod(void *p)
{
	if (!p) return;
	EPNODE *	bep = (*(RcontribMod *)p).binv;
	if (bep) epfree(bep, true);
	efree(p);
}

// Get format identifier
const char *
formstr(int f)
{
	switch (f) {
	case 'a': return("ascii");
	case 'f': return("float");
	case 'd': return("double");
	case 'c': return(NCSAMP==3 ? COLRFMT : SPECFMT);
	}
	return("unknown");
}

// Standard file data share function
RdataShare *
fileDataShare(const char *name, RCOutputOp op, size_t siz)
{
	if (op == RCOrecover && access(name, R_OK|W_OK) < 0) {
		sprintf(errmsg, "cannot recover from '%s'", name);
		error(SYSTEM, errmsg);
		return NULL;
	}
	RdataShare *	rds = new RdataShareFile(name, RSDOflags[op],
						 siz*(op != RCOforce));
						 
	if (!rds || (op == RCOforce && rds->Resize(siz) < siz)) {
		delete rds;
		sprintf(errmsg, "cannot create %lu byte output file '%s'",
					(unsigned long)siz, name);
		error(SYSTEM, errmsg);
		return NULL;
	}
	return rds;
}

// Memory-mapped data share function
RdataShare *
mapDataShare(const char *name, RCOutputOp op, size_t siz)
{
	if (op == RCOrecover && access(name, R_OK|W_OK) < 0) {
		sprintf(errmsg, "cannot recover from '%s'", name);
		error(SYSTEM, errmsg);
		return NULL;
	}
	RdataShare *	rds = new RdataShareMap(name, RSDOflags[op],
						 siz*(op != RCOforce));

	if (!rds || (op == RCOforce && rds->Resize(siz) < siz)) {
		delete rds;
		sprintf(errmsg, "cannot create %lu byte output map '%s'",
					(unsigned long)siz, name);
		error(SYSTEM, errmsg);
		return NULL;
	}
	return rds;
}

// Set output format ('f', 'd', or 'c')
bool
RcontribSimulManager::SetDataFormat(int ty)
{
	if (outList) {
		error(INTERNAL, "cannot call SetDataFormat() after AddModifier()");
		return false;
	}
	switch (ty) {
	case 'f':
		dsiz = sizeof(float)*NCSAMP;
		break;
	case 'd':
		dsiz = sizeof(double)*NCSAMP;
		break;
	case 'c':
		dsiz = LSCOLR;
		break;
	default:
		sprintf(errmsg, "unsupported output format '%c'", ty);
		error(INTERNAL, errmsg);
		return false;
	}
	dtyp = ty;
	return true;
}

// static call-back for rcontrib ray-tracing
int
RcontribSimulManager::RctCall(RAY *r, void *cd)
{
	int	i;

	if (!r->ro || (i = r->ro->omod) == OVOID)
		return 0;		// hit nothing

	if (r->rsrc >= 0 && source[r->rsrc].so != r->ro)
		return 0;		// shadow ray not on source

	const char *		mname = objptr(i)->oname;
	RcontribSimulManager *	rcp = (RcontribSimulManager *)cd;
	RcontribMod *		mp = (RcontribMod *)lu_find(&rcp->modLUT,mname)->data;
	if (!mp)
		return 0;		// not in our modifier list

	if (rcp->HasFlag(RCcontrib)) {	// pre-emptive check for zero
		for (i = NCSAMP; i--; )
			if (r->rcoef[i]*r->rcol[i] > FTINY)
				break;
		if (i < 0)
			return 0;	// zero contribution
	} else if (sintens(r->rcoef) <= FTINY)
		return 0;		// zero coefficient

	int	bi = 0;			// get bin index
	if (mp->binv) {
		worldfunc(RCCONTEXT, r);
		set_eparams(mp->params);
		double		bval = evalue(mp->binv);
		if (bval <= -.5)
			return 0;	// silently ignore negative bin index
		bi = int(bval + .5);
	}
	DCOLORV *	dvp = (*mp)[bi];
	if (!dvp) {
		sprintf(errmsg, "bad bin number for '%s' (%d ignored)", mname, bi);
		error(WARNING, errmsg);
		return 0;
	}
	SCOLOR		contr;
	raycontrib(contr, r, PRIMARY);		// compute coefficient
	if (rcp->HasFlag(RCcontrib))
		smultscolor(contr, r->rcol);	// -> value contribution

	for (i = 0; i < NCSAMP; i++)
		*dvp++ += contr[i];		// accumulate color/spectrum
	return 1;
}

// check for given format type in string, return last char position
static int
hasFormat(const char *tst, const char *match)
{
	static const char	allPoss[] = "%diouxXfFeEgGaAcsb";
	const char *		s = tst;
	while (*s) {
		if (*s++ != '%')
			continue;
		while (!strchr(allPoss, *s))
			s++;
		if (strchr(match, *s))
			break;
		s++;
	}
	return (*s != '\0')*(s - tst);
}

// Add a modifier and arguments, opening output file(s)
bool
RcontribSimulManager::AddModifier(const char *modn, const char *outspec,
				const char *prms, const char *binval, int bincnt)
{
	if (!modn | !outspec || !*modn | !*outspec) {
		error(WARNING, "ignoring bad call to AddModifier()");
		return false;
	}
	if (*outspec == '!') {
		error(USER, "command output not supported by RcontribSimulManager");
		return false;
	}
	if (!nChan) {				// initial call?
		if ((xres < 0) | (yres <= 0)) {
			error(USER, "xres, yres must be set before first modifier");
			return false;
		}
		if (!SetDataFormat(dtyp))
			return false;
		nChan = NCSAMP;
	} else if (nChan != NCSAMP) {
		error(USER, "# spectral channels must be fixed in AddModifier()");
		return false;
	}
	if (Ready()) {
		error(INTERNAL, "call to AddModifier() after PrepOutput()");
		return false;
	}
	LUENT *		lp = lu_find(&modLUT, modn);
	if (lp->data) {
		sprintf(errmsg, "duplicate modifier '%s'", modn);
		error(USER, errmsg);
		return false;
	}
	if (!lp->key)			// create new entry
		lp->key = strcpy((char *)emalloc(strlen(modn)+1), modn);

	RcontribMod *	mp = NewRcMod(prms, binval, bincnt);
	if (!mp) return false;
	lp->data = (char *)mp;
	const int	modndx = hasFormat(outspec, "sb");
	const int	binndx = hasFormat(outspec, "diouxX");
	int		bin0 = 0;
	char		fnbuf[512];
	if (!modndx || (binndx > 0) & (modndx > binndx))
		sprintf(fnbuf, outspec, bin0, modn);
	else
		sprintf(fnbuf, outspec, modn, bin0);
	RcontribOutput *	olast = NULL;
	RcontribOutput *	op;
	for (op = outList; op; op = (olast=op)->next) {
		if (strcmp(op->GetName(), fnbuf))
			continue;
		if (modndx) {		// this ain't right
			sprintf(errmsg, "output name collision for '%s'", fnbuf);
			error(USER, errmsg);
			return false;
		}
		if ((binndx > 0) & (xres > 0)) {
			sprintf(fnbuf, outspec, ++bin0);
			continue;	// each bin goes to another image
		}
		mp->coffset = op->rowBytes;	// else add to what's there
		break;
	}
	if (!op) {			// create new output channel?
		op = new RcontribOutput(fnbuf);
		if (olast) olast->next = op;
		else outList = op;
	}
	mp->opl = op;			// first (maybe only) output channel
	if (modndx)			// remember modifier if part of name
		op->omod = lp->key;
	if (binndx) {			// append output image/bin list
		op->rowBytes += dsiz;
		op->obin = bin0;
		for (bincnt = 1; bincnt < mp->nbins; bincnt++) {
			if (!modndx || (binndx > 0) &  (modndx > binndx))
				sprintf(fnbuf, outspec, bin0+bincnt, modn);
			else
				sprintf(fnbuf, outspec, modn, bin0+bincnt);
			if (!op->next) {
				olast = op;
				olast->next = op = new RcontribOutput(fnbuf);
				if (modndx) op->omod = lp->key;
				op->obin = bin0+bincnt;
			} else {
				op = op->next;
				CHECK(op->obin != bin0+bincnt, CONSISTENCY,
					"bin number mismatch in AddModifier()");
			}
			CHECK(op->rowBytes != mp->coffset, CONSISTENCY,
					"row offset mismatch in AddModifier()");
			op->rowBytes += dsiz;
		}
	} else				// else send all results to this channel
		op->rowBytes += mp->nbins*dsiz;
	return true;
}

// Add a file of modifiers with associated arguments
bool
RcontribSimulManager::AddModFile(const char *modfn, const char *outspec,
						const char *prms,
						const char *binval, int bincnt)
{
	char *		path = getpath(const_cast<char *>(modfn),
					getrlibpath(), R_OK);
	FILE *		fp;
	if (!path || !(fp = fopen(path, "r"))) {
		sprintf(errmsg, "cannot %s modifier file '%s'",
				path ? "open" : "find", modfn);
		error(SYSTEM, errmsg);
		return false;
	}
	char		mod[MAXSTR];
	while (fgetword(mod, sizeof(mod), fp))
		if (!AddModifier(mod, outspec, prms, binval, bincnt)) {
			fclose(fp);
			return false;
		}
	fclose(fp);
	return true;
}

// call-back to check if modifier has been loaded
static int
checkModExists(const LUENT *lp, void *p)
{
	OBJECT	mod = modifier(lp->key);

	if ((mod != OVOID) & (mod < nsceneobjs))
		return 1;

	sprintf(errmsg, "tracked modifier '%s' not found in main scene", lp->key);
	error(WARNING, errmsg);
	return 0;
}

// Prepare output channels and return # completed rows
int
RcontribSimulManager::PrepOutput()
{
	if (!outList || !RtraceSimulManager::Ready()) {
		error(INTERNAL, "PrepOutput() called before octree & modifiers set");
		return -1;
	}
	if (!cdsF) {
		error(INTERNAL, "missing RdataShare constructor call (cdsF)");
		return -1;
	}
	if (lu_doall(&modLUT, checkModExists, NULL) < 0)
		return -1;

	outList->nRows = yres * (xres + !xres);	// all outputs have same #rows
	int	remWarnings = 20;
	for (RcontribOutput *op = outList; op; op = op->next) {
		if (op->rData) {
			error(INTERNAL, "output channel already open in PrepOutput()");
			return -1;
		}
		op->nRows = outList->nRows;
		op->rData = (*cdsF)(op->ofname, outOp,
					GetHeadLen()+1024 + op->nRows*op->rowBytes);
		freeqstr(op->ofname); op->ofname = NULL;
		if (outOp == RCOrecover) {
			int	rd = op->CheckHeader(this);
			if (rd < 0)
				return -1;
			if (rd >= op->nRows) {
				if (remWarnings > 0) {
					sprintf(errmsg, "recovered output '%s' is complete",
							op->GetName());
					error(WARNING, --remWarnings ? errmsg : "etc...");
				}
				rd = op->nRows;
			}
			if (!rInPos | (rInPos > rd)) rInPos = rd;
		} else if (!op->NewHeader(this))
			return -1;
						// make sure there's room
		if (op->rData->GetSize() < op->begData + op->nRows*op->rowBytes &&
				!op->rData->Resize(op->begData + op->nRows*op->rowBytes))
			return -1;		// calls error() for us
	}
	rowsDone.NewBitMap(outList->nRows);	// create row completion map
	rowsDone.ClearBits(0, rInPos, true);
	return nrDone = rInPos;
}

// Create header in open write-only channel
bool
RcontribOutput::NewHeader(const RcontribSimulManager *rcp)
{
	int			headlim = rcp->GetHeadLen() + 1024;
	char *			hdr = (char *)rData->GetMemory(0, headlim, 0);
	int			esiz;
	const int		etyp = rcp->GetFormat(&esiz);

	strcpy(hdr, HDRSTR);
	strcat(hdr, "RADIANCE\n");
	begData = strlen(hdr);
	strcpy(hdr+begData, rcp->GetHeadStr());
	begData += rcp->GetHeadLen();
	if (omod) {
		sprintf(hdr+begData, "MODIFIER=%s\n", omod);
		begData += strlen(hdr+begData);
	}
	if (obin >= 0) {
		sprintf(hdr+begData, "BIN=%d\n", obin);
		begData += strlen(hdr+begData);
	}
	strcpy(hdr+begData, ROWZEROSTR);
	rowCountPos = begData+LNROWSTR;
	begData += sizeof(ROWZEROSTR)-1;
	if (!rcp->xres | (rowBytes > esiz)) {
		sprintf(hdr+begData, "NCOLS=%d\n", int(rowBytes/esiz));
		begData += strlen(hdr+begData);
	}
	sprintf(hdr+begData, "%s%d\n", NCOMPSTR, NCSAMP);
	begData += strlen(hdr+begData);
	if (NCSAMP > 3) {
		sprintf(hdr+begData, "%s %g %g %g %g\n", WLSPLTSTR,
				WLPART[0], WLPART[1], WLPART[2], WLPART[3]);
		begData += strlen(hdr+begData);
	}
	if (etyp != 'c') {
		sprintf(hdr+begData, "%s%d\n", BIGEND, nativebigendian());
		begData += strlen(hdr+begData);
	}
	int	align = 0;
	switch (etyp) {
	case 'f':
		align = sizeof(float);
		break;
	case 'd':
		align = sizeof(double);
		break;
	case 'c':
		break;
	default:
		error(INTERNAL, "unsupported data type in NewHeader()");
		return false;
	}
	strcpy(hdr+begData, FMTSTR);	// write format string
	begData += strlen(hdr+begData);
	strcpy(hdr+begData, formstr(etyp));
	begData += strlen(hdr+begData);
	if (align)			// align data at end of header
		while ((begData+2) % align)
			hdr[begData++] = ' ';
	hdr[begData++] = '\n';		// EOL for data format
	hdr[begData++] = '\n';		// end of nominal header
	if ((rcp->xres > 0) & (rowBytes == esiz)) {	// tack on resolution string?
		sprintf(hdr+begData, PIXSTDFMT, rcp->yres, rcp->xres);
		begData += strlen(hdr+begData);
	}
	return rData->ReleaseMemory(hdr, RDSwrite);
}

// find string argument in header for named variable
static const char *
findArgs(const char *hdr, const char *vnm, int len)
{
	const char *	npos = strnstr(hdr, vnm, len);

	if (!npos) return NULL;

	npos += strlen(vnm);		// find start of (first) argument
	len -= npos - hdr;
	while (len-- > 0 && (*npos == '=') | isspace(*npos))
		if (*npos++ == '\n')
			return NULL;

	return (len >= 0) ? npos : NULL;
}

// Load and check header in read/write channel
int
RcontribOutput::CheckHeader(const RcontribSimulManager *rcp)
{
	int		esiz;
	const int	etyp = rcp->GetFormat(&esiz);
	const int	maxlen = rcp->GetHeadLen() + 1024;
	char *		hdr = (char *)rData->GetMemory(0, maxlen, RDSread);
	const char *	cp;
						// find end of header
	if (!hdr || !(cp = strnstr(hdr, "\n\n", maxlen))) {
		sprintf(errmsg, "cannot find end of header in '%s'", GetName());
		error(USER, errmsg);
		return -1;
	}
	begData = cp - hdr + 1;			// increment again at end
						// check # components
	if (((cp = findArgs(hdr, NCOMPSTR, begData)) ? atoi(cp) : 3) != NCSAMP) {
		sprintf(errmsg, "expected %s%d in '%s'", NCOMPSTR, NCSAMP, GetName());
		error(USER, errmsg);
		return -1;
	}
						// check format
	if (!(cp = findArgs(hdr, FMTSTR, begData)) ||
				strncmp(cp, formstr(etyp), strlen(formstr(etyp)))) {
		sprintf(errmsg, "expected %s%s in '%s'", FMTSTR, formstr(etyp), GetName());
		error(USER, errmsg);
		return -1;
	}
						// check #columns
	if ((cp = findArgs(hdr, "NCOLS=", begData)) ? (atoi(cp)*esiz != rowBytes)
						    : (!rcp->xres | (rowBytes > esiz))) {
		sprintf(errmsg, "expected NCOLS=%d in '%s'", int(rowBytes/esiz), GetName());
		error(USER, errmsg);
		return -1;
	}
						// find row count
	if (!(cp = findArgs(hdr, "NROWS=", begData))) {
		sprintf(errmsg, "missing NROWS in '%s'", GetName());
		error(USER, errmsg);
		return -1;
	}
	rowCountPos = cp - hdr;
	int	rlast = atoi(cp);
	begData++;				// advance past closing EOL
	if ((rcp->xres > 0) & (rowBytes == esiz)) {	// check/skip resolution string?
		char	rbuf[64];
		sprintf(rbuf, PIXSTDFMT, rcp->yres, rcp->xres);
		int	rlen = strlen(rbuf);
		if (strncmp(rbuf, hdr+begData, rlen)) {
			sprintf(errmsg, "bad resolution string in '%s'", GetName());
			error(USER, errmsg);
			return -1;
		}
		begData += rlen;
	}
	// XXX assume the rest is OK: endianness, wavelength splits, modifier, bin
	return rData->ReleaseMemory(hdr, 0) ? rlast : -1;
}

// Rewind calculation (previous results unchanged)
bool
RcontribSimulManager::ResetRow(int r)
{
	if (!rowsDone.Length() | (0 > r) | (r >= rInPos)) {
		error(WARNING, "ignoring bad call to ResetRow()");
		return (r == rInPos);
	}
	FlushQueue();			// finish current and reset
	for (RcontribOutput *op = outList; op; op = op->next)
		if (!op->SetRowsDone(r))
			return false;

	rowsDone.ClearBits(r, rInPos-r, false);
	nrDone = rInPos = r;
	return true;
}

// call-back for averaging each modifier's contributions to assigned channel(s)
static int
putModContrib(const LUENT *lp, void *p)
{
	RcontribMod *			mp = (RcontribMod *)lp->data;
	const RcontribSimulManager *	rcp = (const RcontribSimulManager *)p;
	const double			sca = 1./rcp->accum;
	const DCOLORV *			dvp = mp->cbin;
	RcontribOutput *		op = mp->opl;
	int				i, n;

	switch (rcp->GetFormat()) {	// conversion based on output type
	case 'd': {
		double *	dvo = (double *)op->InsertionP(mp->coffset);
		for (n = mp->nbins; n--; ) {
			for (i = 0; i < NCSAMP; i++)
				*dvo++ = *dvp++ * sca;
			if ((op->obin >= 0) & (n > 0))
				dvo = (double *)(op = op->Next())->InsertionP(mp->coffset);
		}
		} break;
	case 'f': {
		float *		fvo = (float *)op->InsertionP(mp->coffset);
		for (n = mp->nbins; n--; ) {
			for (i = 0; i < NCSAMP; i++)
				*fvo++ = float(*dvp++ * sca);
			if ((op->obin >= 0) & (n > 0))
				fvo = (float *)(op = op->Next())->InsertionP(mp->coffset);
		}
		} break;
	case 'c': {
		COLRV *		cvo = (COLRV *)op->InsertionP(mp->coffset);
		for (n = mp->nbins; n--; ) {
			SCOLOR	scol;
			for (i = 0; i < NCSAMP; i++)
				scol[i] = COLORV(*dvp++ * sca);
			scolor_scolr(cvo, scol);
			cvo += LSCOLR;
			if ((op->obin >= 0) & (n > 0))
				cvo = (COLRV *)(op = op->Next())->InsertionP(mp->coffset);
		}
		} break;
	default:
		error(CONSISTENCY, "unsupported output type in putModContrib()");
		return -1;
	}
						// clear for next tally
	memset(mp->cbin, 0, sizeof(DCOLORV)*NCSAMP*mp->nbins);
	return 1;
}

// Add a ray/bundle to compute next record
int
RcontribSimulManager::ComputeRecord(const FVECT orig_direc[])
{
	if (!Ready())
		return 0;
	if (rInPos >= outList->nRows) {
		error(WARNING, "ComputeRecord() called after last record");
		return 0;
	}
	if (nkids > 0) {			// in parent process?
		int	k = GetChild(false);	// updates output rows
		if (k < 0) return -1;		// someone died?
		RowAssignment	rass;
		rass.row = kidRow[k] = rInPos++;
		rass.ac = accum;
		if (write(kid[k].w, &rass, sizeof(rass)) != sizeof(rass) ||
				writebuf(kid[k].w, orig_direc, sizeof(FVECT)*2*accum)
						!= sizeof(FVECT)*2*accum) {
			error(SYSTEM, "cannot write to child; dead process?");
			return -1;
		}
		return accum;			// tracing/output happens in child
	}
	if (NCSAMP != nChan) {
		error(INTERNAL, "number of color channels must remain fixed");
		return -1;
	}
						// actual work is done here...
	if (EnqueueBundle(orig_direc, accum) < 0)
		return -1;

	RcontribOutput *	op;		// prepare output buffers
	for (op = outList; op; op = op->next)
		if (!op->GetRow(rInPos))
			return -1;
						// convert averages & clear
	if (lu_doall(&modLUT, putModContrib, this) < 0)
		return -1;
						// write buffers
	for (op = outList; op; op = op->next)
		op->DoneRow();

	if (!nkids)				// update row if solo process
		UpdateRowsDone(rInPos++);	// XXX increment here is critical

	return accum;
}

// Get next available child, returning index or -1 if forceWait & idling
int
RcontribSimulManager::GetChild(bool forceWait)
{
	if (nkids <= 0)
		return -1;
						// take inventory
	int	pn, n = 0;
	fd_set	readset, errset;
	FD_ZERO(&readset); FD_ZERO(&errset);
	for (pn = nkids; pn--; ) {
		if (kidRow[pn] < 0) {		// child already ready?
			if (forceWait) continue;
			return pn;		// good enough
		}
		FD_SET(kid[pn].r, &readset);	// will check on this one
		FD_SET(kid[pn].r, &errset);
		if (kid[pn].r >= n)
			n = kid[pn].r + 1;
	}
	if (!n)					// every child is idle?
		return -1;
						// wait on "busy" child(ren)
	while ((n = select(n, &readset, NULL, &errset, NULL)) <= 0)
		if (errno != EINTR) {
			error(SYSTEM, "select call failed in GetChild()");
			return -1;
		}
	char	buf[sizeof(ROW_DONE)] = "X";
	pn = -1;				// get flags set by select
	for (n = nkids; n--; )
		if (kidRow[n] >= 0 &&
				FD_ISSET(kid[n].r, &readset) |
				FD_ISSET(kid[n].r, &errset)) {
						// check for error
			if (FD_ISSET(kid[n].r, &errset) ||
					read(kid[n].r, buf, sizeof(ROW_DONE)) <= 0 ||
					memcmp(buf, ROW_DONE, sizeof(ROW_DONE)))
				return -1;
						// update output row counts
			UpdateRowsDone(kidRow[n]);
			kidRow[n] = -1;		// flag child available
			pn = n;
		}
	return pn;
}

// Update row completion bitmap and update outputs (does not change rInPos)
bool
RcontribSimulManager::UpdateRowsDone(int r)
{
	if (!rowsDone.TestAndSet(r)) {
		error(WARNING, "redundant call to UpdateRowsDone()");
		return false;
	}
	GetRowFinished();
	if (nrDone <= r)
		return true;			// nothing to update, yet
	for (RcontribOutput *op = outList; op; op = op->next)
		if (!op->SetRowsDone(nrDone))
			return false;
	return true;				// up-to-date
}

// Run rcontrib child process (never returns)
void
RcontribSimulManager::RunChild()
{
	FVECT *		vecList = NULL;
	RowAssignment	rass;
	ssize_t		nr;

	accum = 0;
	errno = 0;
	while ((nr = read(0, &rass, sizeof(rass))) == sizeof(rass)) {
		if (!rass.ac) {
			error(CONSISTENCY, "bad accumulator count in child");
			exit(1);
		}
		if (rass.ac > accum) {
			efree(vecList);
			vecList = (FVECT *)emalloc(sizeof(FVECT)*2*rass.ac);
		}
		accum = rass.ac;
		rInPos = rass.row;

		if (readbuf(0, vecList, sizeof(FVECT)*2*accum) !=
					sizeof(FVECT)*2*accum)
			break;

		if (ComputeRecord(vecList) <= 0)
			exit(1);
						// signal this row is done
		if (write(1, ROW_DONE, sizeof(ROW_DONE)) != sizeof(ROW_DONE))
			exit(1);
	}
	if (nr) {
		error(SYSTEM, "read error in child process");
		exit(1);
	}
	exit(0);
}

// Add child processes as indicated
bool
RcontribSimulManager::StartKids(int n2go)
{
	if ((n2go <= 0) | (nkids + n2go <= 1))
		return false;

	if (!nkids)
		cow_memshare();			// preload objects

	n2go += nkids;		// => desired brood size
	kid = (SUBPROC *)erealloc(kid, sizeof(SUBPROC)*n2go);
	kidRow = (int32 *)erealloc(kidRow, sizeof(int32)*n2go);

	fflush(stdout);				// shouldn't use, anyway
	while (nkids < n2go) {
		kid[nkids] = sp_inactive;
		int	rv = open_process(&kid[nkids], NULL);
		if (!rv) {			// in child process?
			while (nkids-- > 0) {
				close(kid[nkids].r);
				close(kid[nkids].w);
			}
			free(kid); free(kidRow);
			kid = NULL; kidRow = NULL;
			RunChild();		// should never return
			_exit(1);
		}
		if (rv < 0) {
			error(SYSTEM, "cannot fork worker process");
			return false;
		}
		kidRow[nkids++] = -1;		// newborn is ready
	}
	return true;
}

// Reap the indicated number of children (all if 0)
int
RcontribSimulManager::StopKids(int n2end)
{
	if (nkids <= 0)
		return 0;
	FlushQueue();
	int	status = 0;
	if (!n2end | (n2end >= nkids)) {	// end all subprocesses?
		status = close_processes(kid, nkids);
		free(kid); free(kidRow);
		kid = NULL; kidRow = NULL;
		nkids = 0;
	} else {				// else shrink family
		int	st;
		while (n2end-- > 0)
			if ((st = close_process(&kid[--nkids])) > 0)
				status = st;
		// could call realloc(), but it hardly matters
	}
	if (status) {
		sprintf(errmsg, "non-zero (%d) status from child", status);
		error(WARNING, errmsg);
	}
	return status;
}

// Set number of computation threads (0 => #cores)
int
RcontribSimulManager::SetThreadCount(int nt)
{
	if (!Ready()) {
		error(INTERNAL, "must call PrepOutput() before SetThreadCount()");
		return 0;
	}
	if (nt < 0)
		return NThreads();
	if (!nt) nt = GetNCores();
	int	status = 0;
	if (nt == 1)
		status = StopKids();
	else if (nt < nkids)
		status = StopKids(nkids-nt);
	else if (nt > nkids)
		StartKids(nt-nkids);
	if (status) {
		sprintf(errmsg, "non-zero (%d) status from child", status);
		error(WARNING, errmsg);
	}
	return NThreads();
}
