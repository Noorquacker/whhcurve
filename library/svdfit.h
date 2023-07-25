#ifndef _SVDFIT_H
#define _SVDFIT_H

#include "tnt/tnt.h"
#include "jama/jama_svd.h"
using namespace TNT;

// npmax and nfmax SHOULD BE MAXIMUM ALLOWED VALUES. FIGURE OUT HOW TO ALLOW SMALLER ONES WHEN NECESSARY

extern const char* SVDFIT_INDEX_OUT_OF_BOUNDS;
extern const char* SVDFIT_TOO_MANY_FUNCTIONS;
extern const char* SVDFIT_TOO_MANY_POINTS;
extern const char* SVDFIT_DIM_MISMATCH;


template <class T>
class SVDFit
{
    public:
    SVDFit(unsigned npmax, unsigned nfmax, T (*phifn)(const unsigned, const T) ) { init(npmax,nfmax,phifn); }

    // pointer to function that determines our fit type
    T (*phi)(const unsigned, const T);

    void fit(const Array1D<T>& x, const Array1D<T>& y, const Array1D<T>& sigma, unsigned nf);
    void fit(const Array1D<T>& x, const Array1D<T>& y, unsigned nf) { fit(x,y,Array1D<T>(x.dim(),1.0),nf); }
    void fit(const Array1D<T>& x, const Array1D<T>& y, const Array1D<T>& sigma) { fit(x,y,sigma,NFmax); }
    void fit(const Array1D<T>& x, const Array1D<T>& y) { fit(x,y,Array1D<T>(x.dim(),1.0),NFmax); }
    
    const T& operator[](unsigned i) const { if(i>=NFfit) throw SVDFIT_INDEX_OUT_OF_BOUNDS; return a[i]; }
    const Array1D<T>& params() const { return a; }
    const T chisquare() const { return chisq; }
    const Array2D<T> covariance() const { return cvm; }
    const Array1D<T> stddev() const
    {
        Array1D<T> std(cvm.dim1());
        for(unsigned i=0; i<std.dim(); i++)
            std[i] = sqrt(cvm[i][i]);
        return std;        
    }

    protected:
    SVDFit() { }    // default constructor
    void init(unsigned npmax, unsigned nfmax, T (*phifn)(const unsigned j, const T x));
    void svdvar();
    void svbksb();

    // number of points and functions given memory allocation (maximum)
    unsigned NPmax;
    unsigned NFmax;

    // number of points and functions in use this fit (fit <= max)
    unsigned NPfit;
    unsigned NFfit;

    Array1D<T> a;	// ending coefficients
    
    Array1D<T> b;	// some temporary storage for svdfit [NPmax]
    Array1D<T> afunc;	// some temporary storage for svdfit [NFmax]
    Array1D<T> tmp;	// some temporary storage for svbksb [NFmax]
    
    // output of svdcmp is the set of matrices u,v,w
    Array2D<T> u;
    Array2D<T> v;
    Array1D<T> w;
    
    // output of svdvar is the covariance matrix
    Array2D<T> cvm;
    
    T chisq;
    bool validCVM;

    T TOL;
};

// x and y are source data, sigma is uncertainty in y, nf is number of functions to use (up to NF)
template <class T>
void SVDFit<T>::fit(const Array1D<T>& x, const Array1D<T>& y, const Array1D<T>& sigma, unsigned nf)
{
    using namespace JAMA;

    #ifdef TNT_NR_DIM_CHECK
    // some checks: everything should have the same dim, and nf should be <= NFmax
    if(nf > (*this).NFmax)
        throw SVDFIT_TOO_MANY_FUNCTIONS;
    if(x.dim() != y.dim() || x.dim() != sigma.dim())
        throw SVDFIT_DIM_MISMATCH;
    if(x.dim() > (*this).NPmax)
        throw SVDFIT_TOO_MANY_POINTS;
    #endif
        
    // set our state
    (*this).validCVM = false;
        
    // everything seems to check out. do fit.
    (*this).NPfit = x.dim();
    (*this).NFfit = nf;
    
	int i,j;
	T wmax,tmp,thresh,sum;
	   	
	for(i=0; i<(*this).NPfit; i++)
	{
	    tmp=1.0/sigma[i];
		for (j=0;j<(*this).NFfit;j++)
        {
		    u[i][j]=(*(this->phi))(j,x[i]) * tmp;
		}
		(*this).b[i]=y[i]*tmp;
    }

    // We should do the decomposition on a subarray of u, not the
    // whole thing. figure out how to do this without destroying everything.

    // these are reference arrays that are fit-size rather than max size,
    // but point to the same memory space as u, v, and w.
    Array2D<T> smallU;
    Array2D<T> smallV;
    Array1D<T> smallW;
    
    smallU = u.subarray(0,NPfit-1,0,NFfit-1);
    smallV = v.subarray(0,NFfit-1,0,NFfit-1);
    smallW = w.subarray(0,NFfit-1);

    Array2D<T> newU, newV;
    Array1D<T> newW;

    // using JAMA implementation of SVD
    SVD<T> decomposer(smallU);
    decomposer.getU(newU);
    decomposer.getV(newV);
    decomposer.getSingularValues(newW);

    // copy into our own arrays
    for(i=0;i<(*this).NPfit;i++)
    {
        for(j=0;j<(*this).NFfit;j++)
        {
            if(j<newU.dim2())
                smallU[i][j] = newU[i][j];
            else
                smallU[i][j] = 0.0; // too many functions!
        }
    }
    for(i=0;i<(*this).NFfit;i++)
        for(j=0;j<(*this).NFfit;j++)
                smallV[i][j] = newV[i][j];
    for(i=0;i<(*this).NFfit;i++)
        smallW[i] = newW[i];                

	wmax=0.0;
	for (j=0; j<(*this).NFfit; j++)
	    if ((*this).w[j] > wmax)
	        wmax=(*this).w[j];

	thresh=(*this).TOL*wmax;
	for (j=0;j<(*this).NFfit;j++)
	    if ((*this).w[j] < thresh)
            (*this).w[j]=0.0;
	
	(*this).svbksb();
	
	(*this).chisq=0.0;
	for (i=0; i<(*this).NPfit; i++) 
	{
	    for(sum=0.0,j=0; j<(*this).NFfit; j++)
            sum += a[j] * (*(this->phi))(j,x[i]);
	    
	    (*this).chisq += (tmp = (y[i]-sum) / sigma[i], tmp*tmp);
	}
} // SVDFit::fit()

template <class T>
void SVDFit<T>::svdvar()
{
    int k,j,i;
    T sum;
    Array1D<T> wti(NFfit);

    for(i=0;i<NFfit;i++) 
    {
        wti[i]=0.0;
        if (w[i] != 0.0)
            wti[i]=1.0/(w[i]*w[i]);
    }
    for(i=0;i<NFfit;i++)
    {
        for(j=0;j<i;j++)
        {
            for(sum=0.0,k=0;k<(*this).NFfit;k++)
                sum += ((*this).v[i][k] * (*this).v[j][k] * wti[k]);
            
            (*this).cvm[j][i] = (*this).cvm[i][j] = sum;
        }
    }
    // set our state
    (*this).validCVM = true;
}

template <class T>
void SVDFit<T>::svbksb()
{
	int jj,j,i;
	double s;
	
	for (j=0; j<(*this).NFfit; j++)
    {
        s=0.0;
        if ((*this).w[j] != 0.0)
        {
            for (i=0;i<(*this).NPfit;i++)
                s += (*this).u[i][j] * (*this).b[i];
            s /= (*this).w[j];
        }
        (*this).tmp[j]=s;
    }
	for(j=0; j<(*this).NFfit; j++)
    {
        s=0.0;
        for(jj=0; jj<(*this).NFfit; jj++)
            s += (*this).v[j][jj] * (*this).tmp[jj];
        (*this).a[j] = s;
    }
}


template <class T>
void SVDFit<T>::init(unsigned npmax, unsigned nfmax, T (*phifn)(const unsigned j, const T x))
{
    // constants
    (*this).TOL = 1.0e-18;

    // allocate some memory for the biggest fit we are asked to process
    (*this).u = Array2D<T>(npmax,nfmax);
    (*this).v = Array2D<T>(nfmax,nfmax);
    (*this).w = Array1D<T>(nfmax);
    (*this).cvm = Array2D<T>(nfmax,nfmax);

    (*this).NFmax = nfmax;
    (*this).NPmax = npmax;    
    
    // set our state
    (*this).validCVM = false;

    (*this).a = Array1D<T>(nfmax);	// ending coefficients
    
    (*this).b = Array1D<T>(npmax);	// some temporary storage for svdfit [NPmax]
    (*this).afunc = Array1D<T>(nfmax);	// some temporary storage for svdfit [NFmax]
    (*this).tmp = Array1D<T>(nfmax);	// some temporary storage for svbksb [NFmax]

    // set fit type
    (*this).phi = phifn;
}

#endif
