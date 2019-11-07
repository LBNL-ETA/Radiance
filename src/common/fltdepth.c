#ifndef lint
static const char RCSid[] = "$Id: fltdepth.c,v 3.2 2019/11/07 23:33:45 greg Exp $";
#endif
/*
 * Function to open floating-point depth file, making sure it's
 * the correct length, and converting from encoded 16-bit
 * depth into a temporary file if needed
 */

#include "rtio.h"
#include "paths.h"

#if defined(_WIN32) || defined(_WIN64)

#include "platform.h"

/* Return file descriptor for open raw float file, or -1 on error */
/* Set expected_length to width*height, or 0 if unknown */
/* We don't translate encoded depth under Windows because of temp file issue */
int
open_float_depth(const char *fname, long expected_length)
{
	int	fd = open(fname, O_RDONLY|O_BINARY);

	if (fd < 0)
		return(-1);
	if (expected_length > 0) {
		off_t	flen = lseek(fd, 0, SEEK_END);
		if (flen != expected_length*sizeof(float)) {
			fprintf(stderr, "%s: expected length is %ld, actual length is %ld\n",
					fname, expected_length, (long)(flen/sizeof(float)));
			close(fd);
			return(-1);
		}
		lseek(fd, 0, SEEK_SET);
	}
	return(fd);
}

#else	/* ! defined(_WIN32) || defined(_WIN64) */

#include <ctype.h>
#include "fvect.h"
#include "depthcodec.h"

#define	FBUFLEN		1024		/* depth conversion batch size */

/* Return file descriptor for open raw float file, or -1 on error */
/* Set expected_length to width*height, or 0 if unknown */
int
open_float_depth(const char *fname, long expected_length)
{
	char		buf[FBUFLEN*sizeof(float)];
	DEPTHCODEC	dc;
	ssize_t		n, nleft;
	int		fd = open(fname, O_RDONLY);

	if (fd < 0)
		return(-1);
	dc.finp = NULL;
	if (expected_length <= 0) {	/* need to sniff file? */
		extern const char	HDRSTR[];
		const int		len = strlen(HDRSTR);
		n = read(fd, buf, len+1);
		if (n < len+1)
			goto gotEOF;
		if (lseek(fd, 0, SEEK_SET) != 0)
			goto seek_error;
		for (n = 0; n < len; n++)
			if (buf[n] != HDRSTR[n])
				break;
		if (n < len || !isprint(buf[len]))	
			return(fd);	/* unknown length raw float... */
	} else {
		off_t	flen = lseek(fd, 0, SEEK_END);
		if (flen < 0 || lseek(fd, 0, SEEK_SET) != 0)
			goto seek_error;
		if (flen == expected_length*sizeof(float))
			return(fd);	/* it's the right length for raw */
	}
	/* We're thinking it's an encoded depth file at this point... */
	set_dc_defaults(&dc);		/* load & convert if we can */
	dc.finp = fdopen(fd, "r");
	if (!dc.finp) {			/* out of memory?? */
		perror("fdopen");
		close(fd);
		return(-1);
	}
	dc.inpname = fname;		/* check encoded depth header */
	dc.hdrflags = HF_HEADIN|HF_RESIN|HF_STDERR;
	if (!process_dc_header(&dc, 0, NULL) || !check_decode_depths(&dc)) {
		fclose(dc.finp);	/* already reported error */
		return(-1);
	}
	if (expected_length > 0 && (long)dc.res.xr*dc.res.yr != expected_length) {
		fprintf(stderr, "%s: expected length is %ld, actual length is %ld\n",
				fname, (long)expected_length, (long)dc.res.xr*dc.res.yr);
		fclose(dc.finp);
		return(-1);
	}
	strcpy(buf, TEMPLATE);		/* create temporary file to hold raw */
	fd = mkstemp(buf);
	if (fd < 0) {
		fputs(buf, stderr);
		fputs(": cannot create temporary file\n", stderr);
		fclose(dc.finp);
		return(-1);
	}
	unlink(buf);			/* unlink it now (Windows forbids) */
	for (nleft = (ssize_t)dc.res.xr*dc.res.yr; nleft > 0; nleft -= FBUFLEN) {
		for (n = 0; n < FBUFLEN; n++) {
			double	d = decode_depth_next(&dc);
			if (d < 0) {
				if (n < nleft)
					goto gotEOF;
				break;
			}
			((float *)buf)[n] = (float)d;
		}
		n *= sizeof(float);
		if (write(fd, buf, n) != n) {
			perror("write");
			fclose(dc.finp);
			close(fd);
			return(-1);
		}
	}
	fclose(dc.finp);		/* all done -- clean up */
	if (lseek(fd, 0, SEEK_SET) != 0)
		goto seek_error;
	return(fd);
gotEOF:
	fputs(fname, stderr);
	fputs(": unexpected end-of-file\n", stderr);
	if (dc.finp) fclose(dc.finp);
	close(fd);
	return(-1);
seek_error:
	perror("lseek");
	close(fd);
	return(-1);
}

#endif
