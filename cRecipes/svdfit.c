#define NRANSI
#include "nrutil.h"
#include "cRecipes.h"
#define TOL 1.0e-5

void svdfit(void *x, double *y, double *sig, int ndata, double *a, int ma,
	double **u, double **v, double *w, double *chisq,
	void (*funcs)(void *, int i,double *, int) )
{
	int j,i;
	double wmax,tmp,thresh,sum,*b,*afunc;

	b=dvector(1,ndata);
	afunc=dvector(1,ma);
	for (i=1;i<=ndata;i++) {
		(*funcs)(x,i,afunc,ma);
		tmp=1.0/sig[i];
		for (j=1;j<=ma;j++) u[i][j]=afunc[j]*tmp;
		b[i]=y[i]*tmp;
	}
	svdcmp(u,ndata,ma,w,v);
	wmax=0.0;
	for (j=1;j<=ma;j++)
		if (w[j] > wmax) wmax=w[j];
	thresh=TOL*wmax;
	for (j=1;j<=ma;j++)
		if (w[j] < thresh) w[j]=0.0;
	svbksb(u,w,v,ndata,ma,b,a);
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x,i,afunc,ma);
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
		*chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
	}
	free_dvector(afunc,1,ma);
	free_dvector(b,1,ndata);
}
#undef TOL
#undef NRANSI
