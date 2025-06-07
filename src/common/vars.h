/* RCSid $Id: vars.h,v 2.12 2025/06/07 05:09:45 greg Exp $ */
/*
 *  Header for programs that load variable files.
 */
#ifndef _RAD_VARS_H_
#define _RAD_VARS_H_
#ifdef __cplusplus
extern "C" {
#endif

typedef struct variable_s {
	const char	*name;	/* variable name */
	short	nick;		/* # characters required for nickname */
	short	nass;		/* # assignments made */
	char	*value;		/* assigned value(s) */
	void	(*fixval)(struct variable_s *);	/* assignment checking function */
} VARIABLE;		/* a variable-value pair */

/**** The following variables should be declared by calling program ****/

extern int	NVARS;		/* total number of variables */

extern VARIABLE	vv[];		/* variable-value pairs */

extern int	nowarn;		/* global boolean to turn warnings off */

/**** The rest is declared in loadvars.c ****/

#define UPPER(c)	((c)&~0x20)	/* ASCII trick */

#define vnam(vc)	(vv[vc].name)
#define vdef(vc)	(vv[vc].nass)
#define vval(vc)	(vv[vc].value)
#define vint(vc)	atoi(vval(vc))
#define vflt(vc)	atof(vval(vc))
#define vlet(vc)	UPPER(vval(vc)[0])
#define vscale		vlet
#define vbool(vc)	(vlet(vc)=='T')

#define HIGH		'H'
#define MEDIUM		'M'
#define LOW		'L'


extern void	loadvars(const char *rfname);
extern int	setvariable(const char *ass, VARIABLE *(*mv)(const char*));
extern VARIABLE	*matchvar(const char *nam);
extern char	*nvalue(int vn, int n);
extern void	checkvalues(void);
extern void	onevalue(VARIABLE *vp);
extern void	catvalues(VARIABLE *vp);
extern int	badmatch(char *tv, char *cv);
extern void	boolvalue(VARIABLE *vp);
extern void	qualvalue(VARIABLE *vp);
extern void	strvalue(VARIABLE *vp);
extern void	intvalue(VARIABLE *vp);
extern void	fltvalue(VARIABLE *vp);
extern int	singlevar(VARIABLE *vp);
extern void	printvars(FILE *fp);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_VARS_H_ */
