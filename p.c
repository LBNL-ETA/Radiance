#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>

main( int argc, char *argv[])
{
	int s,mb,q,a;
	unsigned char *mem,b;

	if( argc != 2 ) {
	  fprintf(stderr,"usage: p MB\n");
	  return(1);
	}
	s= 1024*1024* (mb=atoi(argv[1]));
	if (!(mem=malloc(s))) {
	  printf("failed\n");
	  return(1);
	}
	for (;;) {
	  for (q=0;q<mb;++q) {
	    a=q*1024*1024+ 1+(int)(1024*1024.0*rand()/(RAND_MAX+1.0)) ;
	    b=255*rand()/(RAND_MAX+1.0);
	    mem[a] = b;
	    if( mem[a] != b ) 
		fprintf(stderr,"hopla %x\n", a);
	  }
	}
	  
	return(0);
	
}
