#ifndef lint
static const char	RCSid[] = "$Id: rholo4.c,v 3.32 2003/06/30 14:59:12 schorsch Exp $";
#endif
/*
 * Holodeck display process communication
 */

#include <sys/uio.h>
#include <string.h>

#include "rholo.h"
#include "rhdisp.h"
#include "rtprocess.h"

#ifndef HDSUF
#define HDSUF		".hdi"
#endif
#ifndef SLAVENAME
#define SLAVENAME	"slave"
#endif

#ifndef FNONBLK
#define FNONBLK		O_NONBLOCK
#endif

static int	inp_flags;
static SUBPROC	dpd;
static FILE	*dpout;


disp_open(dname)		/* open the named display driver */
char	*dname;
{
	char	buf[sizeof(HDGRID)+512], fd0[8], fd1[8], *cmd[5], *sfn;
	int	i, n, len;

	if (!strcmp(dname, SLAVENAME)) {
		dpd.r = 0;		/* read from stdin */
		dpout = stdout;		/* write to stdout */
		dpd.running = 0; /* we're the slave procees */
	} else {
					/* get full display program name */
#ifdef DEVPATH
		sprintf(buf, "%s/%s%s", DEVPATH, dname, HDSUF);
#else
		sprintf(buf, "dev/%s%s", dname, HDSUF);
#endif
					/* dup stdin and stdout */
		if (readinp)
			sprintf(fd0, "%d", dup(0));
		else
			strcpy(fd0, "-1");
		sprintf(fd1, "%d", dup(1));
					/* start the display process */
		cmd[0] = buf;
		cmd[1] = froot; cmd[2] = fd1; cmd[3] = fd0;
		cmd[4] = NULL;
		i = open_process(&dpd, cmd);
		if (i <= 0)
			error(USER, "cannot start display process");
		if ((dpout = fdopen(dpd.w, "w")) == NULL)
			error(SYSTEM, "problem opening display pipe");
					/* close dup'ed stdin and stdout */
		if (readinp)
			close(atoi(fd0));
		close(atoi(fd1));
	}
	dpd.w = -1;		/* causes ignored error in close_process() */
	inp_flags = 0;
				/* check if outside */
	if (vdef(OBSTRUCTIONS) && vbool(OBSTRUCTIONS))
		disp_result(DS_OUTSECT, 0, NULL);
				/* send eye separation if specified */
	if (vdef(EYESEP)) {
		sprintf(buf, "%.9e", vflt(EYESEP));
		disp_result(DS_EYESEP, strlen(buf)+1, buf);
	}
				/* write out hologram grids & octrees */
	for (i = 0; hdlist[i] != NULL; i++) {
		memcpy(buf, (void *)hdlist[i], sizeof(HDGRID));
		len = sizeof(HDGRID);
		n = vdef(GEOMETRY);
		sfn = i<n ? nvalue(GEOMETRY,i) :
				n ? nvalue(GEOMETRY,n-1) : vval(OCTREE);
		strcpy(buf+len, sfn);
		len += strlen(sfn) + 1;
		n = vdef(PORTS);
		sfn = i<n ? nvalue(PORTS,i) : n ? nvalue(PORTS,n-1) : "";
		strcpy(buf+len, sfn);
		len += strlen(sfn) + 1;
		disp_result(DS_ADDHOLO, len, buf);
	}
	disp_flush();
}


disp_packet(p)			/* display a packet */
register PACKHEAD	*p;
{
	disp_result(DS_BUNDLE, packsiz(p->nr), (char *)p);
}


disp_check(block)		/* check display process */
int	block;
{
	MSGHEAD	msg;
	int	n;
	char	*buf = NULL;

	if (dpout == NULL)
		return(-1);
					/* flush display output */
	disp_flush();
					/* check read blocking */
	if (block != (inp_flags == 0)) {
		inp_flags = block ? 0 : FNONBLK;
		if (fcntl(dpd.r, F_SETFL, inp_flags) < 0)
			goto fcntlerr;
	}
					/* read message header */
	n = read(dpd.r, (char *)&msg, sizeof(MSGHEAD));
	if (n != sizeof(MSGHEAD)) {
		if (n >= 0) {
			dpout = NULL;
			error(USER, "display process died");
		}
		if (errno != EAGAIN & errno != EINTR)
			goto readerr;
		return(2);		/* acceptable failure */
	}
	if (msg.nbytes) {		/* get the message body */
		if (msg.nbytes < 0)
			error(INTERNAL, "anti-message from display process");
		buf = (char *)malloc(msg.nbytes);
		if (buf == NULL)
			error(SYSTEM, "out of memory in disp_check");
		if (inp_flags != 0 && fcntl(dpd.r, F_SETFL, inp_flags=0) < 0)
			goto fcntlerr;
		if (readbuf(dpd.r, buf, msg.nbytes) != msg.nbytes)
			goto readerr;
	}
	switch (msg.type) {		/* take appropriate action */
	case DR_BUNDLE:		/* new bundle to calculate */
		if (msg.nbytes != sizeof(PACKHEAD))
			error(INTERNAL, "bad DR_BUNDLE from display process");
		bundle_set(BS_ADD, (PACKHEAD *)buf, 1);
		break;
	case DR_NEWSET:		/* new calculation set */
		if (msg.nbytes % sizeof(PACKHEAD))
			error(INTERNAL, "bad DR_NEWSET from display process");
		if (msg.nbytes)
			disp_result(DS_STARTIMM, 0, NULL);
		bundle_set(BS_NEW, (PACKHEAD *)buf, msg.nbytes/sizeof(PACKHEAD));
		if (msg.nbytes) {
			disp_result(DS_ENDIMM, 0, NULL);
			disp_flush();
		}
		break;
	case DR_ADDSET:		/* add to calculation set */
		if (!msg.nbytes)
			break;
		if (msg.nbytes % sizeof(PACKHEAD))
			error(INTERNAL, "bad DR_ADDSET from display process");
		disp_result(DS_STARTIMM, 0, NULL);
		bundle_set(BS_ADD, (PACKHEAD *)buf, msg.nbytes/sizeof(PACKHEAD));
		disp_result(DS_ENDIMM, 0, NULL);
		disp_flush();
		break;
	case DR_ADJSET:		/* adjust calculation set members */
		if (!msg.nbytes)
			break;
		if (msg.nbytes % sizeof(PACKHEAD))
			error(INTERNAL, "bad DR_ADJSET from display process");
		disp_result(DS_STARTIMM, 0, NULL);
		bundle_set(BS_ADJ, (PACKHEAD *)buf, msg.nbytes/sizeof(PACKHEAD));
		disp_result(DS_ENDIMM, 0, NULL);
		disp_flush();
		break;
	case DR_DELSET:		/* delete from calculation set */
		if (!msg.nbytes)
			break;
		if (msg.nbytes % sizeof(PACKHEAD))
			error(INTERNAL, "bad DR_DELSET from display process");
		bundle_set(BS_DEL, (PACKHEAD *)buf, msg.nbytes/sizeof(PACKHEAD));
		break;
	case DR_VIEWPOINT:	/* set target eye position */
		if (msg.nbytes != sizeof(VIEWPOINT))
			error(INTERNAL, "bad DR_VIEWPOINT from display process");
		copystruct(&myeye, (VIEWPOINT *)buf);
		break;
	case DR_ATTEN:		/* block for priority request */
		if (msg.nbytes)
			error(INTERNAL, "bad DR_ATTEN from display process");
					/* send acknowledgement */
		disp_result(DS_ACKNOW, 0, NULL);
		return(disp_check(1));	/* block on following request */
	case DR_KILL:		/* kill computation process(es) */
		if (msg.nbytes)
			error(INTERNAL, "bad DR_KILL from display process");
		if (nprocs > 0)
			done_rtrace();
		else
			error(WARNING, "no rtrace process to kill");
		break;
	case DR_RESTART:	/* restart computation process(es) */
		if (msg.nbytes)
			error(INTERNAL, "bad DR_RESTART from display process");
		if (ncprocs > nprocs)
			new_rtrace();
		else if (nprocs > 0)
			error(WARNING, "rtrace already runnning");
		else
			error(WARNING, "holodeck not open for writing");
		break;
	case DR_CLOBBER:	/* clobber holodeck */
		if (msg.nbytes)
			error(INTERNAL, "bad DR_CLOBBER from display process");
		if (force <= 0 | ncprocs <= 0)
			error(WARNING, "request to clobber holodeck denied");
		else {
			error(WARNING, "clobbering holodeck contents");
			hdclobber(NULL);
		}
		break;
	case DR_SHUTDOWN:	/* shut down program */
		if (msg.nbytes)
			error(INTERNAL, "bad DR_SHUTDOWN from display process");
		return(0);		/* zero return signals shutdown */
	case DR_NOOP:		/* do nothing */
		break;
	default:
		error(INTERNAL, "unrecognized request from display process");
	}
	if (msg.nbytes)			/* clean up */
		free(buf);
	return(1);			/* normal return value */
fcntlerr:
	error(SYSTEM, "cannot change display blocking mode");
readerr:
	error(SYSTEM, "error reading from display process");
}


int
disp_close()			/* close our display process */
{
	int	rval;

	if (dpout == NULL)
		return(-1);
	myeye.rng = 0;
	disp_result(DS_SHUTDOWN, 0, NULL);
	fclose(dpout);
	dpout = NULL;
	return(dpd.running ? close_process(&dpd) : 0);

}


disp_result(type, nbytes, p)	/* queue result message to display process */
int	type, nbytes;
char	*p;
{
	MSGHEAD	msg;
					/* consistency checks */
#ifdef DEBUG
	if (nbytes < 0 || nbytes > 0 & p == NULL)
		error(CONSISTENCY, "bad buffer handed to disp_result");
#endif
	if (dpout == NULL)
		return;
	msg.type = type;
	msg.nbytes = nbytes;
	fwrite((char *)&msg, sizeof(MSGHEAD), 1, dpout);
	if (nbytes > 0)
		fwrite(p, 1, nbytes, dpout);
}


disp_flush()			/* flush output to display */
{
	if (fflush(dpout) < 0)
		error(SYSTEM, "error writing to the display process");
}
