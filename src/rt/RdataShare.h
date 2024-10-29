/* RCSid $Id: RdataShare.h,v 2.1 2024/10/29 00:36:54 greg Exp $ */
/*
 *  RdataShare.h
 *
 *	Radiance classes for interprocess data sharing
 *
 *	Errors end up calling error() from rterror.h, so
 *	checking return values is not normally needed.
 *
 *  Created by Greg Ward on 10/14/2024.
 */

#ifndef RdataShare_h
#define RdataShare_h

#include <sys/types.h>

/// Channel creation and i/o flags
enum {RDSexcl=1, RDSextend=2, RDSread=4, RDSwrite=8};

/*
 * When data sharing object is created, flag meanings are:
 *	RDSexcl =	Create a new file, or fail if it exists
 *	RDSextend =	Set data extent to match size parameter
 *	RDSread =	Data may be read
 *	RDSwrite =	Data may be written
 *
 * The RDSread flag should be given to GetMemory() if a data
 * read operation is desired.  Similarly, the RDSwrite flag
 * should be given when writing data in ReleaseMemory().
 * Memory-mapped files always read *and* write data if opened
 * in the appropriate mode.  (Beware that writing to read-only data
 * or reading write-only data will generate a system fault.)
 * The RDSextend flag may be given with RDSread in GetMemory()
 * and RDSwrite in ReleaseMemory() to allow the file size to
 * increase as needed for a given operation.
 * In some implementations, the RDSexcl flag may guarantee
 * data consistency by placing locks on the requested memory areas.
 * In such cases, it may make sense to include the RDSwrite flag
 * in GetMemory() calls if the intention is to write new results
 * in ReleaseMemory(), so an exclusive lock will be held.
 */

/// Data share object types returned by GetType()
enum RDSType {RDSTanonMap=1, RDSTfileMap, RDSTfile,
		RDSTcust1, RDSTcust2, RDSTcust3, RDSTcust4};

extern const char	RDSnoname[];

/// Abstract base class for shared memory object
class RdataShare {
protected:
	char *		chName;			// channel name
	size_t		osiz;			// current object size
	int		mode;			// open mode
public:
			RdataShare() {
				chName = NULL; osiz = 0; mode = 0;
			}
	virtual		~RdataShare();
			/// Get channel name, or RDSnoname if anonymous
	const char *	GetName() const {
				return chName ? chName : RDSnoname;
			}
			/// Get R/W flags
	int		GetMode() const {
				return mode & (RDSread|RDSwrite);
			}
			/// Get current object size
	size_t		GetSize() const {
				return osiz;
			}
			/// Return data sharing type
	virtual RDSType	GetType() const = 0;
			/// Attempt to extend or shrink object (adjust if 0)
	virtual size_t	Resize(size_t new_siz = 0) = 0;
			/// Get data buffer
	virtual void *	GetMemory(size_t offs, size_t len, int fl = RDSread) = 0;
			/// Return data buffer
	virtual bool	ReleaseMemory(void *bp, int fl = RDSwrite) = 0;
};

/// Memory sharing implemented with mmap'ed file (or anonymous area)
class RdataShareMap : public RdataShare {
	void *		mmorg;			// memory-mapped origin
	int		bufCount;		// count of allocated buffers
public:
			RdataShareMap(const char *name, int flags, size_t siz = 0);
	virtual		~RdataShareMap();
	RDSType		GetType() const {
				return chName ? RDSTfileMap : RDSTanonMap;
			}
	size_t		Resize(size_t new_siz = 0);
	void *		GetMemory(size_t offs, size_t len, int fl);
	bool		ReleaseMemory(void *dp, int fl);
};

struct RDSbuffer;	// Private struct for buffer list

/// Memory sharing implemented with simple file
class RdataShareFile : public RdataShare {
	int		fd;			// open file descriptor
	RDSbuffer *	blist;			// allocated buffer list
public:
			RdataShareFile(const char *name, int flags, size_t siz = 0);
	virtual		~RdataShareFile();
	RDSType		GetType() const {
				return RDSTfile;
			}
	size_t		Resize(size_t new_siz = 0);
	void *		GetMemory(size_t offs, size_t len, int fl);
	bool		ReleaseMemory(void *dp, int fl);
};

#endif		// RdataShare_h
