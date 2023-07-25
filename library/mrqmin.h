#ifndef _MRQMIN_H
#define _MRQMIN_H

#include "tnt/tnt.h"
#include "tntnr.h"
#include "mrqcof.h"
#include "covsrt.h"
#include "gaussj.h"

using namespace TNT;


extern const char* MRQMIN_DIM_MISMATCH1;
extern const char* MRQMIN_DIM_MISMATCH2;
extern const char* MRQMIN_BAD_STATE_START;
extern const char* MRQMIN_BAD_STATE_STEP;
extern const char* MRQMIN_BAD_STATE_STOP;
extern const char* MRQMIN_BAD_STATE_READ;

const unsigned LM_STATE_INIT = 1;
const unsigned LM_STATE_START = 2;
const unsigned LM_STATE_STEP = 3;
const unsigned LM_STATE_STOP = 4;

template <class T>
class LevenbergMarquardtFit
{
    public:
    LevenbergMarquardtFit(unsigned xndatamax, unsigned xmamax) { init(xndatamax, xmamax); }

    const T chisquare() const { return chisq; }
    const Array2D<T> covariance() const { return covar; }
    const Array1D<T> stddev() const
    {
        Array1D<T> std(covar.dim1());
        for(unsigned i=0; i<std.dim(); i++)
            std[i] = sqrt(covar[i][i]);
        return std;        
    }
    T operator [](unsigned int n) const
    {
        #ifdef TNT_NR_DIM_CHECK
        if(state < LM_STATE_START)
            throw MRQMIN_BAD_STATE_READ;
        #endif
        return a[n];
    }
    const Array1D<T>& params() const
    {
        #ifdef TNT_NR_DIM_CHECK
        if(state < LM_STATE_START)
            throw MRQMIN_BAD_STATE_STEP;
        #endif
        return a;
    }
    

    const Array1D<T> fit(const Array1D<T> qx, const Array1D<T> qy, const Array1D<T> qsig, const Array1D<T> guessa, T (*xfuncs)(T, const Array1D<T>, Array1D<T>) )
    {
        start(qx,qy,qsig,guessa,xfuncs);
        T lastchisq;
        int stoppingcondition = 0;
        
        // conditions to continue stepping
        while(stoppingcondition < 2)
        {
            // do the step
            lastchisq=chisq;
            step();

            // assume first that we're prepared to get closer to finishing
            stoppingcondition++;
            
            // then if we're still not there, dash the loop's hopes of completion
            if(chisq-lastchisq < -0.01 && (chisq-lastchisq)/chisq < -1e-3)
                stoppingcondition = 0;
            if(chisq-lastchisq >= 0)
                stoppingcondition = 0;

        }
        stop();
        return a;
    }

    // gets called to start a fit, returns chisq
    T start(const Array1D<T> qx, const Array1D<T> qy, const Array1D<T> qsig, const Array1D<T> guessa, T (*xfuncs)(T, const Array1D<T>, Array1D<T>) )
    {
        x = qx;
        y = qy;
        sig = qsig;
        funcs = xfuncs;
        a = guessa.copy();

        #ifdef TNT_NR_DIM_CHECK
        if(state < LM_STATE_INIT)
            throw MRQMIN_BAD_STATE_START;

        // checks: x=y=sig=>ndata, a=covar1=covar2=alpha1=alpha2=>ma
        if(x.dim() != y.dim() || y.dim() != sig.dim())
            throw MRQMIN_DIM_MISMATCH1;

        if(a.dim() != covar.dim1() || covar.dim1() != covar.dim2() || covar.dim2() != alpha.dim1() || alpha.dim1() != alpha.dim2())
            throw MRQMIN_DIM_MISMATCH2;
        #endif

        ndata = x.dim();
        ma = a.dim();

        int j;
        for(mfit=0,j=0; j<ma; j++)
            if(ia[j]) mfit++;
        oneda = Array2D<T>(mfit,1);
        alamda = 0.001;
        mrqcof<T>(x, y, sig, a, ia, alpha, beta, chisq, funcs);
        ochisq = chisq;
        for(j=0; j<ma; j++) atry[j] = a[j];

        state = LM_STATE_START;
        nstep++;
        return chisq;        
    }

    
    T step() // returns chisq
    {
        #ifdef TNT_NR_DIM_CHECK
        if(state < LM_STATE_START)
            throw MRQMIN_BAD_STATE_STEP;
        #endif

        int j,k,l;

        for (j=0; j<mfit; j++)
        {
            for(k=0; k<mfit; k++) covar[j][k] = alpha[j][k];
            covar[j][j]=alpha[j][j]*(1.0+alamda);
            oneda[j][0]=beta[j];
        }
        gaussj<T>(covar.subarray(0,mfit-1,0,mfit-1),oneda);
        for(j=0; j<mfit; j++)
            da[j]=oneda[j][0];
        if(alamda == 0.0)   // if we have been asked to finish and report results
        {
            covsrt<T>(covar,ia,mfit);
            return chisq;
        }
        for(j=-1,l=0; l<ma; l++)
            if(ia[l]) atry[l]=a[l]+da[++j];
        mrqcof<T>(x, y, sig, atry, ia, covar, da, chisq, funcs);
        if(chisq < ochisq)
        {
            alamda *= 0.1;
            ochisq = chisq;

            for(j=0; j<mfit; j++)
            {
                for(k=0; k<mfit; k++) alpha[j][k] = covar[j][k];
                beta[j] = da[j];
            }
            for(l=0; l<ma; l++) a[l] = atry[l];
        }
        else
        {
            alamda *= 10.0;
            chisq = ochisq;
        }

        state = LM_STATE_STEP;
        nstep++;
        return chisq;       
    }
    
    T stop()
    {
        #ifdef TNT_NR_DIM_CHECK
        if(state < LM_STATE_START)
            throw MRQMIN_BAD_STATE_READ;
        #endif

        alamda = 0.0;
        T temp = step();
        state = LM_STATE_STOP;
        return temp;
    }
    
    void dump() { this->dump(std::cout); }
    void dump(ostream& out)
    {
        cout << 
            "state: " << state << endl <<
            "nstep: " << nstep << endl <<
            "ndatamax: " << ndatamax << endl <<
            "mamax: " << mamax << endl <<
            "x: " << x << endl <<
            "y: " << y << endl <<
            "sig: " << sig << endl <<
            "a: " << a << endl <<
            "ia: " << ia << endl <<
            "dyda: " << dyda << endl <<
            "covar: " << covar << endl <<
            "alpha: " << alpha << endl <<
            "beta: " << beta << endl <<
            "atry: " << atry << endl <<
            "da: " << da << endl <<
            "oneda: " << oneda << endl <<
            "chisq: " << chisq << endl <<
            "ochisq: " << ochisq << endl <<
            "alamda: " << alamda << endl <<
            "mfit: " << mfit << endl <<
            "ndata: " << ndata << endl <<
            "ma: " << ma << endl << endl << endl;
    }

    protected:

    // program state
    int state, nstep;

    // data members
    unsigned ndatamax, mamax;
    Array1D<T> x;
    Array1D<T> y;
    Array1D<T> sig;
    Array1D<T> a;
    Array1D<int> ia;
	Array1D<T> dyda;
    Array2D<T> covar;
    Array2D<T> alpha;
	Array1D<T> beta;
	Array1D<T> atry;
	Array1D<T> da;
	Array2D<T> oneda;
	T chisq, ochisq, alamda;
    int mfit,ndata,ma;
    T (*funcs)(T x, const Array1D<T> a, Array1D<T> dyda);


    void init(unsigned xndatamax, unsigned xmamax)
    {
        ndatamax = xndatamax;
        mamax = xmamax;
        
        atry = Array1D<T>(mamax);
        beta = Array1D<T>(mamax);
        da = Array1D<T>(mamax);
        dyda = Array1D<T>(mamax);
        ia = Array1D<int>(mamax,1);
        covar = Array2D<T>(mamax,mamax);
        alpha = Array2D<T>(mamax,mamax);
        
        state = LM_STATE_INIT;
        nstep = -1;
    }
};

#endif 
// _MRQMIN_H
