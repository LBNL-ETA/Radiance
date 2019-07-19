#ifndef lint
static const char RCSid[] = "$Id: pictool.c,v 2.6 2019/07/19 17:37:56 greg Exp $";
#endif
#include "pictool.h"
#include "g3sphere.h"
#include <math.h>
#include "rtio.h"

#define        index                  strchr



int		pict_init(pict* p)
{
	p->pic = NULL;
	p->use_lut = 0;
	p->lut = NULL;
	p->pjitter = 0;
	p->comment[0] = '\0';
#ifdef PICT_GLARE
	p->pinfo = NULL;
	p->glareinfo = g3fl_create(PICT_GLSIZE);
#endif
	p->valid_view = 1;
	p->view = stdview;
	if (!pict_update_view(p))
		return 0;
	p->resol.rt = PIXSTANDARD;
	return pict_resize(p,512,512);
}

void	pict_new_gli(pict* p)
{
	int i;

	g3Float* gli;
	gli = g3fl_append_new(p->glareinfo);
	for(i=0;i<PICT_GLSIZE;i++)
		gli[i] = 0.0;
}

void	pict_set_comment(pict* p,const char* c)
{
	strcpy(p->comment,c);
}

pict*	pict_get_copy(pict* p)
{
	pict* res;

	if (p->use_lut) {
		fprintf(stderr,"pict_get_copy: can't copy in lut mode");
		return NULL;
	}
	res = pict_create();
	pict_resize(res,pict_get_xsize(p),pict_get_ysize(p));
	res->view = p->view;
	res->resol = p->resol;
	res->pjitter = p->pjitter;
	return res;
}

pict*	pict_create()
{
	pict* p = (pict*) malloc(sizeof(pict));
	
	if (!pict_init(p)) {
		pict_free(p);
		return NULL;
	}
	return p;
}

void	pict_free(pict* p)
{
	if (p->pic)
		free(p->pic);
#ifdef PICT_GLARE
	if (p->pinfo)
		free(p->pinfo);
	g3fl_free(p->glareinfo);
	p->pinfo = NULL;
#endif
	p->pic = NULL;
	p->lut = NULL;
	free(p);
}
static int	exp_headline(char* s,void* exparg) {
	double* exposure;
        if(strstr(s, EXPOSSTR) != NULL )  {
                fprintf(stderr,"EXP AND tab\n");
                }

	exposure = (double*) exparg;
	if (isexpos(s)) {
		(*exposure) *= exposval(s);
	}
	return 0;
}

static	double read_exposure(FILE* fp)
{
	double exposure = 1;


	if ((getheader(fp, exp_headline, &exposure)) < 0) {
		fprintf(stderr,"error reading header\n");
		return 0;
	}
	return exposure;
}

int	pict_update_view(pict* p)
{
	char* ret;
	if ((ret = setview(&(p->view)))) {
		fprintf(stderr,"pict_update_view: %s\n",ret);
		return 0;
	}
	return 1;
}

int	pict_set_pj(pict* p,double pj)
{
	if ((pj < 0) || (pj > 1))
		fprintf(stderr,"pict_set_pj: warning jitter factor out of range\n");
	p->pjitter = pj;
	return 1;
}

int	pict_setvup(pict* p,FVECT vup)
{
	VCOPY(p->view.vup,vup);
	return pict_update_view(p);
}

int		pict_trans_cam(pict* p,g3ATrans t)
{
	g3Vec o = g3v_create();
	g3Vec d = g3v_create();
	g3Vec u = g3v_create();
	
fprintf(stderr,"vdpre %f %f %f\n",p->view.vdir[0],p->view.vdir[1],p->view.vdir[2]);
fprintf(stderr,"vupre %f %f %f\n",p->view.vup[0],p->view.vup[1],p->view.vup[2]);
	g3v_fromrad(o,p->view.vp);
	g3v_fromrad(d,p->view.vdir);
	g3v_fromrad(u,p->view.vup);

	g3v_add(d,d,o);
	g3v_add(u,u,o);
	g3at_apply(t,o);
	g3at_apply(t,u);
	g3at_apply(t,d);
	g3v_sub(d,d,o);
	g3v_sub(u,u,o);
	g3v_torad(p->view.vp,o);
	g3v_torad(p->view.vdir,d);
	g3v_torad(p->view.vup,u);
fprintf(stderr,"vd %f %f %f\n",p->view.vdir[0],p->view.vdir[1],p->view.vdir[2]);
fprintf(stderr,"vu %f %f %f\n",p->view.vup[0],p->view.vup[1],p->view.vup[2]);

	g3v_free(o);
	g3v_free(d);
	g3v_free(u);
	return pict_update_view(p);
}

int	pict_resize(pict* p,int x,int y)
{
	p->resol.xr = x;
	p->resol.yr = y;

	if (!(p->pic = (COLOR*) realloc(p->pic,(p->resol.xr)*(p->resol.yr)*sizeof(COLOR)))) {
		fprintf(stderr,"Can't alloc picture\n");
		return 0;
	}
	if (pict_lut_used(p)) {
		if (!(p->lut = (COLOR*) realloc(p->lut,(p->resol.xr)*(p->resol.yr)*sizeof(COLOR)))) {
			fprintf(stderr,"Can't alloc LUT\n");
			return 0;
		}
	}
#ifdef PICT_GLARE
	if (!(p->pinfo = (pixinfo*) 
			realloc(p->pinfo,(p->resol.xr)*p->resol.yr*sizeof(pixinfo)))) {
		fprintf(stderr,"Can't alloc picture info\n");
		return 0;
	}
/*	fprintf(stderr,"resize %d %d %d\n",p->resol.xr,p->resol.yr,sizeof(pixinfo));*/
#endif
	return 1;
}


int		pict_zero(pict* p)
{
	int x,y,i;

	for(x=0;x<pict_get_xsize(p);x++)
		for(y=0;y<pict_get_ysize(p);y++) 
			for(i=0;i<3;i++)
				pict_get_color(p,x,y)[i] = 0;
	return 1;
}

struct hinfo {
	VIEW	*hv;
	int	ok;
	double exposure;
};


static int
gethinfo(				/* get view from header */
	char  *s,
	void  *v
)
{
        if(strstr(s, EXPOSSTR) != NULL && strstr(s, "\t") != NULL)  {
                
                fprintf(stderr,"error: header contains invalid exposure entry!!!!\n");
                fprintf(stderr,"check exposure and correct header setting !\n");
                
                
                fprintf(stderr,"stopping !!!!\n");
                exit(1);
                }
	if (isview(s) && sscanview(((struct hinfo*)v)->hv, s) > 0) {
		((struct hinfo*)v)->ok++;
	} else if (isexpos(s)) {
/*fprintf(stderr,"readexp %f %f\n",((struct hinfo*)v)->exposure,exposval(s));*/
		((struct hinfo*)v)->exposure *= exposval(s);
	}
	return(0);
}

int		pict_read_fp(pict* p,FILE* fp)
{
	struct hinfo	hi;
	int x,y,yy;


	hi.hv = &(p->view);
	hi.ok = 0;
	hi.exposure = 1;
	getheader(fp, gethinfo, &hi);
/*fprintf(stderr,"expscale %f\n",hi.exposure);*/

	if (!(pict_update_view(p)) || !(hi.ok))
		p->valid_view = 0;
/*	printf("dir %f %f %f\n",p->view.vdir[0],p->view.vdir[1],p->view.vdir[2]); */
	if (fgetsresolu(&(p->resol),fp) < 0) {
		fprintf(stderr,"No resolution given\n");
		return 0;
	}
	if (!(p->resol.rt & YMAJOR)) {
		x = p->resol.xr;
		p->resol.xr = p->resol.yr;
		p->resol.yr = x;
		p->resol.rt |= YMAJOR;
	}
	if (!pict_resize(p,p->resol.xr,p->resol.yr))
		return 0;
	for(yy=0;yy<p->resol.yr;yy++) {
		y = (p->resol.rt & YDECR) ? p->resol.yr - 1 - yy : yy;
		if (freadscan((p->pic + y*p->resol.xr),p->resol.xr,fp) < 0) {
			fprintf(stderr,"error reading line %d\n",y);
			return 0;
		}
		for(x=0;x<p->resol.xr;x++) {
            scalecolor(pict_get_color(p,x,y),1.0/hi.exposure);
		}
	}
#ifdef PICT_GLARE
	pict_update_evalglare_caches(p);
#endif
	return 1;
}

#ifdef PICT_GLARE
void pict_update_evalglare_caches(pict* p)
{
	int x,y,yy;
	float lum;
	for(yy=0;yy<p->resol.yr;yy++) {
		y = (p->resol.rt & YDECR) ? p->resol.yr - 1 - yy : yy;
		for(x=0;x<p->resol.xr;x++) {
			pict_get_omega(p,x,y) = pict_get_sangle(p,x,y);
            pict_get_dir(p,x,y,pict_get_cached_dir(p,x,y));
            pict_get_gsn(p,x,y) = 0;
            pict_get_pgs(p,x,y) = 0;
/*    make picture grey     */
            lum=luminance(pict_get_color(p,x,y))/179.0;
            pict_get_color(p,x,y)[RED]=lum;
            pict_get_color(p,x,y)[GRN]=lum;
            pict_get_color(p,x,y)[BLU]=lum;
		}
	}
}

#endif

int		pict_read(pict* p,char* fn)
{
	FILE* fp;

	if (!(fp = fopen(fn,"rb"))) {
		fprintf(stderr,"Can't open %s\n",fn);
		return 0;
	}
	if (!(pict_read_fp(p,fp)))
		return 0;
	fclose(fp);
	return 1;
}

static	void	ch_pct(char* line)
{
	char* pct;

	while ((pct = index(line,',')))
		*pct = '.';
}

char*		nextline(FILE* fp)
{
	static char* line = NULL;
	static int lsize = 500;
	int done;
	long fpos;

	if (!line) {
		if (!(line = (char*) malloc((lsize + 1)*sizeof(char)))) {
			fprintf(stderr,"allocation problem in nextline\n");
			return NULL;
		}
	}
	done = 0;
	fpos = ftell(fp);
	while(!done) {
		if (!fgets(line,lsize,fp))
			return NULL;
		if (!index(line,'\n')) {
			lsize *= 2;
			if (!(line = (char*) realloc(line,(lsize + 1)*sizeof(char)))) {
				fprintf(stderr,"allocation problem in nextline\n");
				return NULL;
			}
			fseek(fp,fpos,SEEK_SET);
		} else {
			done = 1;
		}
	}
	return line;
}

static	int	read_tt_txt(pict* p,int bufsel,int channel,char* fn)
{
	FILE* fp;
	char* line = NULL;
	int x,y,vsize,hsize;
	float* col;
	float c;
	char* curr;
	char** endp;
	COLOR* picbuf;

	if (!(fp = fopen(fn,"r"))) {
		fprintf(stderr,"pict_read_tt_txt: Can't open file %s\n",fn);
		return 0;
	}
	if (!(line = nextline(fp))) {
		fprintf(stderr,"pict_read_tt_txt: reading from file failed\n");
		return 0;
	}
	if (strncmp(line,"float",strlen("float"))) {
		fprintf(stderr,"pict_read_tt_txt: Unexpected format\n");
		return 0;
	}
	if (!(line = nextline(fp))) {
		fprintf(stderr,"pict_read_tt_txt: reading from file failed\n");
		return 0;
	}
	sscanf(line,"%d %d %d %d",&x,&vsize,&y,&hsize);
	vsize = vsize - x;
	hsize = hsize - y;
	if (!pict_resize(p,vsize,hsize))
		return 0;
	endp = malloc(sizeof(char*));
	picbuf = (bufsel == PICT_LUT) ? p->lut : p->pic;
	for(x=0;x<vsize;x++) {
		if (!(line = nextline(fp))) {
			fprintf(stderr,"pict_read_tt_txt: too few lines %d of %d\n",x,vsize);
			free(endp);
			return 0;
		}
		ch_pct(line);
		curr = line;
		for(y=0;y<hsize;y++) {
			col = pict_colbuf(p,picbuf,x,y);
			c = strtod(curr,endp);
			if (channel & PICT_CH1)
				col[0] = c;
			if (channel & PICT_CH2)
				col[1] = c;
			if (channel & PICT_CH3)
				col[2] = c;
			if (curr == (*endp)) {
				fprintf(stderr,"pict_read_tt_txt: too few values in line %d\n",x);	
				free(endp);
				return 0;
			}
			curr = (*endp);
		}
	}
	free(endp);
	return 1;
}

int		pict_read_tt_txt(pict* p,char* fn)
{
	return read_tt_txt(p,PICT_VIEW,PICT_ALLCH,fn);
}

int		pict_setup_lut(pict* p,char* theta_fn,char* phi_fn,char* omega_fn)
{
	int xs,ys;

	p->use_lut = 1;
	if (!read_tt_txt(p,PICT_LUT,PICT_CH1,theta_fn)) {
		fprintf(stderr,"pict_setup_lut: error reading theta file\n");
		p->use_lut = 0;
		return 0;
	}
	xs = pict_get_xsize(p);
	ys = pict_get_ysize(p);
	if (!read_tt_txt(p,PICT_LUT,PICT_CH2,phi_fn)) {
		fprintf(stderr,"pict_setup_lut: error reading phi file\n");
		p->use_lut = 0;
		return 0;
	}
	if ((xs != pict_get_xsize(p)) || (ys != pict_get_ysize(p))) {
		fprintf(stderr,"pict_setup_lut: different resolutions for lut files\n");
		p->use_lut = 0;
		return 0;
	}
	if (!read_tt_txt(p,PICT_LUT,PICT_CH3,omega_fn)) {
		fprintf(stderr,"pict_setup_lut: error reading omega file\n");
		p->use_lut = 0;
		return 0;
	}
	if ((xs != pict_get_xsize(p)) || (ys != pict_get_ysize(p))) {
		fprintf(stderr,"pict_setup_lut: different resolutions for lut files\n");
		p->use_lut = 0;
		return 0;
	}
	return 1;
}
		
int		pict_update_lut(pict* p)
{
	g3ATrans transf;
	g3Vec d;
	g3Vec u;
	int x,y;
	
	d = g3v_fromrad(g3v_create(),p->view.vdir);
	u = g3v_fromrad(g3v_create(),p->view.vup);
	g3at_ctrans(transf,d,u);
	for(x=0;x<p->resol.xr;x++) {
		for(y=0;y<p->resol.yr;y++) {
			d[G3S_RAD] = 1.0;
			d[G3S_THETA] = DEG2RAD(pict_colbuf(p,p->lut,x,y)[0]);
			d[G3S_PHI] = DEG2RAD(pict_colbuf(p,p->lut,x,y)[1]);
			g3s_sphtocc(d,d);
			g3at_apply(transf,d);
			g3s_cctosph(d,d);
			pict_colbuf(p,p->lut,x,y)[0] = RAD2DEG(d[G3S_THETA]);
			pict_colbuf(p,p->lut,x,y)[1] = RAD2DEG(d[G3S_PHI]);
		}
	}
	g3v_free(d);
	g3v_free(u);
	return 1;
}
	

int		pict_write_dir(pict* p,FILE* fp)
{
	int x,y;
	FVECT dir,orig;

	for(y=0;y<p->resol.yr;y++) {
		for(x=0;x<p->resol.xr;x++) {
			if (!(pict_get_ray(p,x,y,orig,dir))) {
				fprintf(fp,"%f %f %f 0.0 0.0 0.0\n",orig[0],orig[1],orig[2]);
			} else {
				fprintf(fp,"%f %f %f %f %f %f\n",orig[0],orig[1],orig[2],dir[0],dir[1],dir[2]);
			}
		}
	}
	return 1;
}

int		pict_read_from_list(pict* p,FILE* fp)
{
	int x,y,sz;
	char ln[1024];
	double val;

	sz = 0;
	for(y=0;y<p->resol.yr;y++) {
		for(x=0;x<p->resol.xr;x++) {
			if (!(fgets(ln,1023,fp))) {
				fprintf(stderr,"pict_read_from_list: Read only %d values\n",sz);
				return 0;
			}
			if (sscanf(ln,"%lg",&val) != 1) {
				fprintf(stderr,"pict_read_from_list: errline %d\n",sz);
			} else {
				sz++;
				setcolor(pict_get_color(p,x,y),val,val,val);
			}
		}
	}
	return 1;
}

int	pict_write(pict* p,char* fn)
{
	FILE* fp;

	if (!(fp = fopen(fn,"wb"))) {
		fprintf(stderr,"Can't open %s for writing\n",fn);
		return 0;
	}
	if (!(pict_write_fp(p,fp))) {
		fclose(fp);
		return 0;
	}
	fclose(fp);
	return 1;
}

int	pict_write_fp(pict* p,FILE* fp)
{
	char res[100];
	int y,yy;
	newheader("RADIANCE",fp);
	fputnow(fp);
	if (strlen(p->comment) > 0)
		fputs(p->comment,fp);
	fputs(VIEWSTR, fp);
	fprintview(&(p->view), fp);
	fprintf(fp,"\n");
	fputformat(COLRFMT,fp);
	fprintf(fp,"\n");

	resolu2str(res,&(p->resol));
	fprintf(fp,"%s",res);
	for(y=0;y<p->resol.yr;y++) {
		yy = (p->resol.rt & YDECR) ? p->resol.yr - 1 - y : y;
		if (fwritescan((p->pic + yy*p->resol.xr),p->resol.xr,fp) < 0) {
			fprintf(stderr,"writing pic-file failed\n");
			return 0;
		}
	}
	return 1;
}

int		pict_get_dir(pict* p,int x,int y,FVECT dir)
{
	FVECT  orig;
	if (pict_lut_used(p)) {
		if (!(p->lut)) {
			fprintf(stderr,"pict_get_dir: no lut defined\n");
			return 0;
		}
		dir[G3S_RAD] = 1.0;
		dir[G3S_THETA] = DEG2RAD(pict_colbuf(p,p->lut,x,y)[0]);
		dir[G3S_PHI] = DEG2RAD(pict_colbuf(p,p->lut,x,y)[1]);
		g3s_sphtocc(dir,dir);
		return 1;
	}
	return pict_get_ray(p,x,y,orig,dir);
}

int		pict_get_dir_ic(pict* p,double x,double y,FVECT dir)
{
	FVECT  orig;
	if (pict_lut_used(p)) {
		int res[2];
		loc2pix(res,&(p->resol),x,y);
		return pict_get_dir(p,res[0],res[1],dir);
	}
	return pict_get_ray_ic(p,x,y,orig,dir);
}

int		pict_get_ray(pict* p,int x,int y,FVECT orig,FVECT dir)
{
	RREAL pc[2];

	pix2loc(pc,&(p->resol),x,p->resol.yr - 1 - y);
	if (pict_lut_used(p)) {
		if (viewray(orig,dir,&(p->view),pc[0],pc[1]) <= 0)
			return 0;
		return pict_get_dir(p,x,y,dir);
	}
	return pict_get_ray_ic(p,pc[0],pc[1],orig,dir);
}

int		pict_get_ray_ic(pict* p,double x,double y,FVECT orig,FVECT dir)
{
	if (pict_lut_used(p)) {
		int res[2];
		loc2pix(res,&(p->resol),x,y);
		return pict_get_ray(p,res[0],res[1],orig,dir);
	}
	x = (x + (0.5*p->pjitter)*(0.5-gb_random())/p->resol.xr);
	y = (y + (0.5*p->pjitter)*(0.5-gb_random())/p->resol.yr);
	if (viewray(orig,dir,&(p->view),x,y) < 0)
		return 0;
	return 1;
}
	

static	int splane_normal(pict* p,double e1x,double e1y,
							double e2x,double e2y,FVECT n)
{
	FVECT e1 = {0,0,0};
	FVECT e2 = {0,0,0};
	
	pict_get_dir_ic(p,e1x,e1y,e1);
	pict_get_dir_ic(p,e2x,e2y,n);
	VSUB(e2,n,e1);
	VCROSS(n,e1,e2);
	if (normalize(n) == 0.0) {
		return 0;
	}
	return 1;
}

double	pict_get_sangle(pict* p,int x,int y)
{
	int i;
	double ang,a,hpx,hpy;
	RREAL pc[2];
	FVECT n[4] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};

	if (pict_lut_used(p)) {
		return pict_colbuf(p,p->lut,x,y)[2];
	}
	pix2loc(pc,&(p->resol),x,y);
	hpx = (0.5/p->resol.xr);
	hpy = (0.5/p->resol.yr);
	pc[0] -= hpx;
	pc[1] -= hpy;
	
	i = splane_normal(p,pc[0],pc[1],pc[0],(pc[1]+2.0*hpy),n[0]);
	i &= splane_normal(p,pc[0],(pc[1]+2.0*hpy),(pc[0]+2.0*hpx),(pc[1]+2.0*hpy),n[1]);
	i &= splane_normal(p,(pc[0]+2.0*hpx),(pc[1]+2.0*hpy),(pc[0]+2.0*hpx),pc[1],n[2]);
	i &= splane_normal(p,(pc[0]+2.0*hpx),pc[1],pc[0],pc[1],n[3]);

	if (!i) {
		return 0;
	}
	ang = 0;
	for(i=0;i<4;i++) {
		a = DOT(n[i],n[(i+1)%4]);
		a = acos(a);
		a = fabs(a);
		ang += M_PI - a;
	}
	ang = ang - 2.0*M_PI;
	if ((ang > (2.0*M_PI)) || ang < 0) {
		fprintf(stderr,"Normal error in pict_get_sangle %f %d %d\n",ang,x,y);
		return -1.0;
	}
	return ang;
}

int		pict_locate(pict* p,FVECT pt,int* x,int* y)
{
	FVECT pp;
	int res[2];

	if (pict_lut_used(p)) {
		fprintf(stderr,"pict_locate: Not implemented for lut mode\n");
		return 0;
	}

	if (viewloc(pp,&(p->view),pt) <= 0)
		return 0;
	loc2pix(res,&(p->resol),pp[0],pp[1]);
	*x = res[0];
	*y = res[1];
	return 1;
}

int	pict_get_sigma(pict* p,double x,double y,FVECT vdir,FVECT vup,double* s)
{
	FVECT pvdir;

	if (!pict_get_dir(p,x,y,pvdir))
		return 0;
	if (normalize(vdir) == 0.0) {
		fprintf(stderr,"get_sigma: view dir is zero\n");
		return 0;
	}
	if (normalize(vup) == 0.0) {
		fprintf(stderr,"get_sigma: view up is zero\n");
		return 0;
	}
	*s = acos(DOT(vdir,pvdir));
	return 1;
}

int	pict_get_tau(pict* p,double x,double y,FVECT vdir,FVECT vup,double* t)
{
	FVECT hv,pvdir; 
	double s;
	int i;
	if (!pict_get_sigma(p,x,y,vdir,vup,&s))
		return 0;
	if (!pict_get_dir(p,x,y,pvdir))
		return 0;
	VCOPY(hv,pvdir);
	normalize(hv);
	for(i=0;i<3;i++){
		hv[i] /= cos(s);
	}
	VSUB(hv,hv,vdir);
	normalize(hv);
	*t = acos(fdot(vup,hv));
	return 1;
}

static	int	ortho_coord(FVECT vdir,FVECT vup,FVECT vv,FVECT hv)
{
	if (normalize(vdir) == 0.0) {
		fprintf(stderr,"ortho_coord: view direction is zero\n");
		return 0;
	}
	if (normalize(vup) == 0.0) {
		fprintf(stderr,"ortho_coord: view direction is zero\n");
		return 0;
	}
	fcross(hv,vdir,vup);
	fcross(vv,vdir,hv);
	return 1;
}

int	pict_get_vangle(pict* p,double x,double y,FVECT vdir,FVECT vup,double* a)
{
	FVECT hv,vv,pvdir; 
	if (!ortho_coord(vdir,vup,vv,hv)) {
		return 0;
	}
	if (!pict_get_dir(p,x,y,pvdir))
		return 0;
	*a = acos(fdot(vv,pvdir)) - (M_PI/2.0) ; 
	return 1;
}

int	pict_get_hangle(pict* p,double x,double y,FVECT vdir,FVECT vup,double* a)
{
	FVECT hv,vv,pvdir; 
	if (!ortho_coord(vdir,vup,vv,hv)) {
		return 0;
	}
	if (!pict_get_dir(p,x,y,pvdir)) {
		return 0;
	}
	*a = ((M_PI/2.0) - acos(fdot(hv,pvdir)));
	return 1;
}
	
int		pict_setorigin(pict* p,FVECT orig)
{
	VCOPY(p->view.vp,orig);
	return pict_update_view(p);
}

int		pict_setdir(pict* p,FVECT dir)
{
	VCOPY(p->view.vdir,dir);
	return pict_update_view(p);
}

int		pict_sethview(pict* p,double ha)
{
	p->view.horiz = ha;
	return pict_update_view(p);
}

int		pict_setvview(pict* p,double va)
{
	p->view.vert = va;
	return pict_update_view(p);
}

int		pict_setvtype(pict* p,char vt)
{
	p->view.type = vt;
	return pict_update_view(p);
}

#ifdef	TEST_PICTTOOL

char*	progname = "pictool";

int	main(int argc,char** argv)
{
	pict* p = pict_create();
	FVECT dir = {0,-1,0};
	FVECT vdir,edir;
	FVECT up = {0,0,1};
	int x,y,i;
	double sang,h,v;

	if (argc < 2) {
		fprintf(stderr,"usage: %s <pic-file>\n",argv[0]);
		return EXIT_FAILURE;
	}
	fprintf(stderr,"readdd\n");
	pict_read(p,argv[1]);
	pict_set_comment(p,"pictool test output\n");
	/*
	pict_get_dir(p,100,100,edir); // YDECR
	for(y=0;y<p->resol.yr;y++) {
		for(x=0;x<p->resol.xr;x++) {
			pict_get_hangle(p,x,y,dir,up,&h);
			pict_get_vangle(p,x,y,dir,up,&v);
			pict_get_dir(p,x,y,vdir);
			if (acos(DOT(vdir,edir)) < DEG2RAD(3))
				setcolor(pict_get_color(p,x,y),h+M_PI/2,v+M_PI/2,1);
			else
				setcolor(pict_get_color(p,x,y),h+M_PI/2,v+M_PI/2,0);
		}
	}
	*/
	printf("resolution: %dx%d\n",pict_get_xsize(p),pict_get_ysize(p));
	sang = 0.0;

	pict_new_gli(p);
	for(x=0;x<pict_get_xsize(p);x++)
		for(y=0;y<pict_get_ysize(p);y++) {
			pict_get_dir(p,x,y,dir);
			if (DOT(dir,p->view.vdir) >= 0.0) {
				pict_get_av_omega(p,0) += pict_get_omega(p,x,y);
			}
			sang = pict_get_sangle(p,x,y);
			for(i=0;i<3;i++)
				if (sang <= 0 && pict_is_validpixel(p,x,y)) {
				pict_get_color(p,x,y)[i] = 1;
				} else {
				pict_get_color(p,x,y)[i] =0;
				}
		}
	printf("s %f\n",pict_get_av_omega(p,0));
	pict_write(p,"test.pic");
	while(1) {
		printf("pos\n");
		if (scanf("%d %d",&x,&y) != 2)
			continue;
		pict_get_dir(p,x,y,dir);
		printf("dir: \t\t %f %f %f\n",dir[0],dir[1],dir[2]);
		printf("solid angle: \t %f \n",pict_get_sangle(p,x,y));
	}
	return EXIT_SUCCESS;
}

#endif
