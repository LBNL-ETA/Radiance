#ifndef lint
static const char RCSid[] = "$Id: RdataShareMap.cpp,v 2.1 2024/10/29 00:36:54 greg Exp $";
#endif
/*
 *  RdataShareMap.cpp
 *
 *	Shared data using memory-mapped file
 *
 *  Created by Greg Ward on 10/14/2024
 */

#include "rtio.h"
#include "rterror.h"
#include "RdataShare.h"
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>

// Create memory-mapped file object
RdataShareMap::RdataShareMap(const char *name, int flags, size_t siz)
{
	mmorg = NULL;
	bufCount = 0;

	if (!(flags & (RDSread|RDSwrite))) {
		error(CONSISTENCY, "RdataShareMap() flags must include RDSread or RDSwrite");
		return;
	}
	if (name && !*name) name = NULL;
	if (!name && !siz | ((flags & (RDSextend|RDSread|RDSwrite)) !=
					(RDSextend|RDSread|RDSwrite))) {
		error(CONSISTENCY, "anonymous memory map must be read/write");
		return;
	}
	if ((flags & (RDSextend|RDSwrite)) == RDSextend) {
		error(CONSISTENCY, "bad RDSextend in RdataShareMap()");
		return;
	}
	int	mmprot = PROT_NONE;
	int	oflags = 0;
	switch (flags & (RDSread|RDSwrite)) {
	case RDSread|RDSwrite:
		mmprot |= PROT_READ|PROT_WRITE;
		oflags |= O_RDWR;
		break;
	case RDSread:
		mmprot |= PROT_READ;
		oflags |= O_RDONLY;
		break;
	case RDSwrite:
		mmprot |= PROT_WRITE;
		oflags |= O_WRONLY;
		break;
	}
	int	fd = -1;
	if (name) {			// opening a shared file
		if (flags & RDSexcl) oflags |= O_CREAT|O_EXCL;
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
			else {
				sprintf(errmsg, "cannot stat '%s'", chName);
				error(SYSTEM, errmsg);
				close(fd);
				return;
			}
		} else if (siz && ftruncate(fd, siz) < 0) {
			sprintf(errmsg, "cannot resize '%s'", name);
			error(SYSTEM, errmsg);
			close(fd);
			return;
		}
	}
	mmorg = (void *)mmap(NULL, siz, mmprot,
				MAP_SHARED|(name ? MAP_FILE : MAP_ANON), fd, 0);
	close(fd);
	if (mmorg == MAP_FAILED) {
		if (name)
			sprintf(errmsg, "cannot map '%s' to memory", name);
		else
			sprintf(errmsg, "cannot map anonymous map of %ld KBytes",
					long(siz/1024));
		error(SYSTEM, errmsg);
		mmorg = NULL;
		return;
	}
	osiz = siz;
	if (name) chName = savqstr(name);
	mode = flags;
}

RdataShareMap::~RdataShareMap()
{
	if (mmorg) munmap(mmorg, osiz);
}

// Attempt to extend or shrink object (adjust if 0)
size_t
RdataShareMap::Resize(size_t new_siz)
{
	if (!mmorg)
		return 0;
	if (new_siz > 0) {
		if (new_siz == osiz)
			return osiz;
		if (!(mode & RDSwrite)) {
			error(CONSISTENCY, "cannot resize read-only map");
			return 0;
		}
	}
	if (bufCount > 0) {
		error(INTERNAL, "cannot resize while memory is checked out");
		return 0;
	}
	if (!chName) {
		if (!new_siz)		// XXX should issue warning?
			return osiz;
		if (new_siz > osiz) {
			error(INTERNAL, "cannot grow anonymous map");
			return 0;
		}
		return osiz = new_siz;	// just pretend we shrank
	}
	if (!new_siz) {			// sync to current file length?
		struct stat	sbuf;
		if (stat(chName, &sbuf) < 0) {
			sprintf(errmsg, "cannot stat '%s'", chName);
			error(SYSTEM, errmsg);
			return 0;
		}
		if (sbuf.st_size <= osiz)
			return osiz = sbuf.st_size;

		new_siz = sbuf.st_size;
	}
	if (new_siz > osiz) {		// need to extend & remap
		int	fd = open(chName, mode&RDSread ? O_RDWR : O_WRONLY);
		if (fd < 0) {
			sprintf(errmsg, "cannot reopen '%s'", chName);
			error(SYSTEM, errmsg);
			return 0;
		}
		if (ftruncate(fd, new_siz) < 0) {
			sprintf(errmsg, "cannot grow '%s'", chName);
			error(SYSTEM, errmsg);
			close(fd);
			return 0;
		}
		munmap(mmorg, osiz);
		mmorg = mmap(NULL, new_siz,
				mode&RDSread ? PROT_READ|PROT_WRITE : PROT_WRITE,
				MAP_SHARED|MAP_FILE, fd, 0);
		close(fd);
		if (mmorg == MAP_FAILED) {
			sprintf(errmsg, "mmap() failed on '%s'", chName);
			error(SYSTEM, errmsg);
			mmorg = NULL;
			return osiz = 0;
		}
	} else if (truncate(chName, new_siz) < 0) {
		sprintf(errmsg, "cannot truncate '%s'", chName);
		return 0;
	}
	return osiz = new_siz;
}

// Get data buffer
void *
RdataShareMap::GetMemory(size_t offs, size_t len, int fl)
{
	if (chName && fl & RDSextend) {	// resize/extend?
		if (offs + len > osiz ? !Resize(offs + len) : !Resize())
			return NULL;
	}
	if (offs + len > osiz) {
		if (chName)
			sprintf(errmsg, "requested block of %lu bytes outside map '%s'",
					(unsigned long)len, chName);
		else
			sprintf(errmsg, "requested block of %lu bytes outside %lu byte map",
					(unsigned long)len, (unsigned long)osiz);
		error(CONSISTENCY, errmsg);
		return NULL;
	}
	++bufCount;
	return (char *)mmorg + offs;
}

// Return data buffer
bool
RdataShareMap::ReleaseMemory(void *dp, int fl)
{
	if ((dp < mmorg) | ((char *)dp >= (char *)mmorg + osiz)) {
		if (chName)
			sprintf(errmsg, "returned block outside map '%s'", chName);
		else
			sprintf(errmsg, "returned block outside %lu byte map",
					(unsigned long)osiz);
		error(CONSISTENCY, errmsg);
		return false;
	}
	bufCount -= (bufCount > 0);
	return true;
}
