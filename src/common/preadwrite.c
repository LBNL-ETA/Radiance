#ifndef lint
static const char	RCSid[] = "$Id: preadwrite.c,v 3.1 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 * Substitute routines for pread(2) and pwrite(2)
 */

extern off_t	lseek();


int
pread(fd, buf, siz, offs)	/* read buffer from an open file */
int	fd;
char	*buf;
unsigned int	siz;
long	pos;
{
	if (lseek(fd, (off_t)offs, 0) != offs)
		return(-1);
	return(read(fd, buf, siz));
				/* technically, we should reset pointer here */
}


int
pwrite(fd, buf, siz, offs)	/* write buffer to an open file */
int	fd;
char	*buf;
unsigned int	siz;
long	pos;
{
	if (lseek(fd, (off_t)offs, 0) != offs)
		return(-1);
	return(write(fd, buf, siz));
				/* technically, we should reset pointer here */
}
