#ifndef lint
static const char RCSid[] = "$Id$";
#endif
/*
 *  RdataShareFile.cpp
 *
 *	Shared data using a regular file
 *
 *  Created by Greg Ward on 10/14/2024
 */

#include "rtio.h"
#include "rterror.h"
#include "RdataShare.h"
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>

const char	RDSnoname[] = "<unnamed_channel>";

RdataShare::~RdataShare()
{
	freeqstr(chName);
}

// Struct needed for tracking allocated buffers
struct RDSbuffer {
	RDSbuffer *	next;			// next in list
	char *		buf;			// allocated buffer
	size_t		len;			// buffer length
	off_t		pos;			// offset in file

			RDSbuffer(size_t nb, RDSbuffer *nxt = NULL);
			~RDSbuffer() {
				delete next;
				if (len > 0) free(buf);
			}
			/// Find buffer matching pointer
	RDSbuffer *	Find(void *p) {
				RDSbuffer *	bp = this;
				while (p != (void *)bp->buf &&
						(bp = bp->next) != NULL)
					;
				return bp;
			}
};

// Allocate an i/o buffer
RDSbuffer::RDSbuffer(size_t nb, RDSbuffer *nxt)
{
	next = nxt;
	buf = NULL;
	pos = -1;
	if ((len = nb) > 0) {
		buf = (char *)malloc(nb);
		if (!buf) {
			sprintf(errmsg, "cannot allocate %lu-byte buffer",
					(unsigned long)nb);
			error(SYSTEM, errmsg);
			len = 0;
			return;
		}
	}
}

// Create memory-mapped file object
RdataShareFile::RdataShareFile(const char *name, int flags, size_t siz)
{
	fd = -1;
	blist = NULL;

	if (!name || !*name) {
		error(CONSISTENCY, "missing file name in RdataShareFile()");
		return;
	}
	if (!(flags & (RDSread|RDSwrite))) {
		error(CONSISTENCY, "RdataShareFile() flags must include RDSread or RDSwrite");
		return;
	}
	if ((flags & (RDSextend|RDSwrite)) == RDSextend) {
		error(CONSISTENCY, "bad RDSextend in RdataShareFile()");
		return;
	}
	int	oflags = O_CLOEXEC;
	switch (flags & (RDSread|RDSwrite)) {
	case RDSread|RDSwrite:
		oflags |= O_RDWR|O_CREAT;
		break;
	case RDSwrite:
		oflags |= O_WRONLY|O_CREAT;
		break;
	case RDSread:
		oflags |= O_RDONLY;
		break;
	}
	if (flags & RDSexcl) oflags |= O_EXCL;
	else if (flags & RDSextend && !siz) oflags |= O_TRUNC;
	fd = open(name, oflags, 0666);
	if (fd < 0) {
		sprintf(errmsg, "cannot open '%s'", name);
		error(SYSTEM, errmsg);
		return;
	}
	if (!(flags & RDSextend)) {
		struct stat	sbuf;
		if (flags & RDSexcl)
			siz = 0;
		else if (fstat(fd, &sbuf) >= 0)
			siz = sbuf.st_size;
		else  {
			sprintf(errmsg, "cannot stat '%s'", name);
			error(SYSTEM, errmsg);
			close(fd); fd = -1;
			return;
		}
	} else if (siz && ftruncate(fd, siz) < 0) {
		sprintf(errmsg, "cannot resize '%s'", name);
		error(SYSTEM, errmsg);
		close(fd); fd = -1;
		return;
	}
	osiz = siz;
	chName = savqstr(name);
	mode = flags;
}

RdataShareFile::~RdataShareFile()
{
	if (fd >= 0) close(fd);
	delete blist;
}

// Attempt to extend or shrink object (adjust if 0)
size_t
RdataShareFile::Resize(size_t new_siz)
{
	if (fd < 0)
		return 0;
	if (new_siz > 0) {
		if (new_siz == osiz)
			return osiz;
		if (!(mode & RDSwrite)) {
			error(CONSISTENCY, "cannot resize read-only file");
			return 0;
		}
	} else {			// sync to current file length
		struct stat	sbuf;
		if (fstat(fd, &sbuf) < 0) {
			sprintf(errmsg, "cannot stat '%s'", chName);
			error(SYSTEM, errmsg);
			return 0;
		}
		return osiz = sbuf.st_size;
	}				// else attempt to resize file
	if (ftruncate(fd, new_siz) < 0) {
		sprintf(errmsg, "cannot resize '%s'", chName);
		return 0;
	}
	return osiz = new_siz;
}

// Get data buffer
void *
RdataShareFile::GetMemory(size_t offs, size_t len, int fl)
{
	if (fd < 0)
		return NULL;
	if (fl & RDSextend && !Resize())	// sync to file length, first?
		return NULL;
	if (fl & mode & RDSread && offs + len > osiz) {
		if (fl & RDSextend) {	// resize if requested
			if (!Resize(offs + len))
				return NULL;
		} else {
			sprintf(errmsg, "requested block of %lu bytes past EOF in '%s'",
						(unsigned long)len, chName);
			error(CONSISTENCY, errmsg);
			return NULL;
		}
	}				// XXX should check for buffer overlap?
	blist = new RDSbuffer(len, blist);
	blist->pos = offs;
	if (fl & mode & RDSread &&	// reading from file?
				pread(fd, blist->buf, len, offs) != len) {
		sprintf(errmsg, "read error on '%s'", chName);
		error(SYSTEM, errmsg);
	}
	return blist->buf;
}

// Return data buffer
bool
RdataShareFile::ReleaseMemory(void *dp, int fl)
{
	if (!blist) return false;
	RDSbuffer *	bp = blist->Find(dp);
	if (!bp) {
		sprintf(errmsg, "return of unallocated block for '%s'", chName);
		error(CONSISTENCY, errmsg);
		return false;
	}
	if ((fl & (RDSwrite|RDSextend)) == RDSwrite && bp->pos + bp->len > osiz) {
		sprintf(errmsg, "write request past EOF in '%s'", chName);
		error(CONSISTENCY, errmsg);
		return false;
	}
	if (fl & mode & RDSwrite) {	// writing to file?
		if (pwrite(fd, bp->buf, bp->len, bp->pos) != bp->len) {
			sprintf(errmsg, "write error on '%s'", chName);
			error(SYSTEM, errmsg);
			return false;
		}
		if (bp->pos + bp->len > osiz)
			osiz = bp->pos + bp->len;
	}
	RDSbuffer	rbuf(0, blist);		// remove from buffer list
	RDSbuffer *	bp2 = &rbuf;
	while (bp2->next != bp)
		bp2 = bp2->next;
	bp2->next = bp->next;
	bp->next = NULL;
	delete bp;				// frees buffer memory
	blist = rbuf.next;
	return true;
}
