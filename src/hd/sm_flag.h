/* RCSid: $Id: sm_flag.h,v 3.3 2003/06/20 00:25:49 greg Exp $ */
/* sm_flag.h */

/* 32 bit FLAGS */
#define F_OFFSET(t)		((t)>>5)
#define F_BIT(t)		((t)&0x1f)
#define F_OP(f,t,op)	        ((f)[F_OFFSET(t)] op (0x1<<F_BIT(t)))
#define IS_FLAG(f,t)		F_OP(f,t,&)
#define SET_FLAG(f,t)		F_OP(f,t,|=)
#define CLR_FLAG(f,t)		F_OP(f,t,&=~)
#define FLAG_BYTES(n)            ((((n)+31) >>5)*sizeof(int32))
