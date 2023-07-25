#ifndef _SAVITZKYGOLAY_H
#define _SAVITZKYGOLAY_H

#include <string>

#include "tnt/tnt.h"
#include "jama/jama_lu.h"

using namespace TNT;
using namespace JAMA;

#include "twave1d.h"
#include "fwave1d.h"
#include "wave1d_operators.h"

#define SAVITZKYGOLAY_INIT  1
#define SAVITZKYGOLAY_LU    2
#define SAVITZKYGOLAY_COEFFICIENTS  3

template <class T>
class SavitzkyGolay
{
    public:
    // constructors
        SavitzkyGolay(const SavitzkyGolay& sg)
        {
            MODE = SAVITZKYGOLAY_INIT;
            m = sg.m;
            nl = sg.nl;
            nr = sg.nr;
            ld = sg.ld;
            b = Array1D<T>();            
            test_params();
            
            lu = 0;
        }
        SavitzkyGolay(unsigned ord, unsigned pts) { params(ord,pts,0); lu = 0; MODE = SAVITZKYGOLAY_INIT; }
        SavitzkyGolay(unsigned ord, unsigned pts, unsigned deriv) { params(ord,pts,deriv); lu = 0; MODE = SAVITZKYGOLAY_INIT; }
        SavitzkyGolay() { lu = 0; MODE = 0;}
        ~SavitzkyGolay() { delete lu; }
        
        /** get/set filter order */
        unsigned order(unsigned _m) { m = _m; test_params(); MODE = SAVITZKYGOLAY_INIT; return m; }
        unsigned order() const { return m; }
        /** get/set number of points */
        unsigned points(unsigned pts) { nl = (pts-1)/2; nr =(pts-1)/2; test_params(); MODE = SAVITZKYGOLAY_INIT; return points(); }
        unsigned points() const { return nl+nr+1; }
        /** set number of points before and after seperately */
        void points(unsigned b, unsigned a) { nl = b; nr = a; test_params(); MODE = SAVITZKYGOLAY_INIT; }
        /** get/set number of points before */
        unsigned pointsbefore(unsigned pts) { nl = pts; test_params(); MODE = SAVITZKYGOLAY_INIT; return nl; }
        unsigned pointsbefore() const { return nl; }
        /** get/set number of points after */
        unsigned pointsafter(unsigned pts) { nr = pts; test_params(); MODE = SAVITZKYGOLAY_INIT; return nr; }
        unsigned pointsafter() const { return nr; }
        /** get/set derivative order */
        unsigned derivative_order(unsigned o) { ld = o; test_params(); MODE = SAVITZKYGOLAY_INIT; return ld; }
        unsigned derivative_order() const { return ld; }
        /** set all params */
        void params(unsigned ord, unsigned pts) { m = ord; points(pts); if(MODE > SAVITZKYGOLAY_LU) MODE = SAVITZKYGOLAY_INIT; else MODE = SAVITZKYGOLAY_INIT; }
        void params(unsigned ord, unsigned pts, unsigned deriv) { m = ord; nl = (pts-1)/2; nr =(pts-1)/2; ld = deriv; test_params(); MODE = SAVITZKYGOLAY_INIT; }

        TWave1D<T,T> filter(const TWave1D<T,T> twave)
        {
            test_params();
            if(MODE < SAVITZKYGOLAY_COEFFICIENTS)
                make_coeff();
//cerr << "running convolution" << endl;            
            FWave1D<T,T> fwave;
            fwave.isPeriodic(false);
            fwave = twave;
            TWave1D<T,T> tkernel(fwave.dim());
            
            int kk;
            int k, mm;
            T sum, fac, derivfac;
            
			derivfac = pow(twave.dt(),-ld);
            // zero out kernel to begin
            for(kk=0; kk<tkernel.dim(); kk++)
                tkernel[kk] = 0.0;
            for(k = -nl; k <= nr; k++)
            {
                sum = b[0];
                fac = 1.0;
                for(mm=0; mm < m; mm++)
                    sum += b[mm+1]*(fac *= k);
                kk = (tkernel.dim()-k) % tkernel.dim();   // wrap it around the end
                tkernel[kk] = sum * derivfac;
//cerr << "kernel point " << k << ", position " << kk << " = " << sum << endl;
            }
//cerr << tkernel;
            FWave1D<T,T> fkernel;
            fkernel.isPeriodic(true);
            fkernel = tkernel;  // fft
            fkernel.copyTime(fwave);
            fkernel.isPeriodic(true);   // again, since copyTime resets it
            
            // still not sure if wave operations are debugged completely,
            // so let's do the convolution explicitly
            fkernel[0] *= fwave[0];
            fkernel[1] *= fwave[1];
            for(kk=2; kk<fkernel.dim(); kk += 2)
            {
                T real = fkernel[kk]*fwave[kk] - fkernel[kk+1]*fwave[kk+1];
                fkernel[kk+1] = fkernel[kk+1]*fwave[kk] + fkernel[kk]*fwave[kk+1];
                fkernel[kk] = real;
            }

//for(int i=0; i<fkernel.dim(); i++)
//{
//    cout << fkernel[i] << "\t" << fwave[i] << "\t" << endl;
//}

            
            fkernel.copyTime(fwave);
            tkernel = fkernel;  // inverse-fft
            return tkernel;
        }

    protected:
        void make_lu()    // find coefficient matrix
        {
            if(MODE < SAVITZKYGOLAY_INIT)
                mode_error();
//cerr << "making lu" << endl;
			// library/savitzkygolay.h:127: warning: unused variable ‘j’                
            // int imj, ipj, j, k, mm;
			int imj, ipj, k, mm;
            T sum;
            Array2D<T> a(this->m+1,this->m+1);
            
            // set up normal equations for least squares
            for(ipj = 0; ipj <= (m << 1); ipj++)
            {
                sum = (ipj ? 0.0 : 1.0);
                for(k=1; k <= nr; k++)
                    sum += pow( (double) k, (double) ipj);
                for(k=1; k <= nl; k++)
                    sum += pow( (double) -k, (double) ipj);
                mm = min(ipj, 2*m - ipj);
                for(imj = -mm; imj <= mm; imj += 2)
                {
                    a[(ipj+imj)/2][(ipj-imj)/2] = sum;
//cerr << "setting a[" << (ipj+imj)/2 << "][" << (ipj-imj)/2 << "] = " << sum << ", ipj=" << ipj << ", imj=" << imj << endl;
                }
            }
            
            // do the LU decomposition
            if(lu) delete lu;
            lu = new LU<T>(a);
            MODE = SAVITZKYGOLAY_LU;
        }
        void make_coeff()
        {
//cerr << "making coefficients" << endl;
            test_params();
            if(MODE < SAVITZKYGOLAY_LU)
                make_lu();
            
            // set up b
            b = Array1D<T>(m+1);
            for(int j=0; j < m+1; j++)
                b[j] = 0.0;
            b[ld] = 1.0;

            b = lu->solve(b);

//cerr << "b = " << b;            
            MODE = SAVITZKYGOLAY_COEFFICIENTS;
        }
        
        void mode_error()
        {
            throw "SavitzkyGolay not in correct mode for operation.";
        }
        void param_error()
        {
            throw "SavitzkyGolay parameters invalid.";
        }
        void test_params()
        {
//cerr << "testing params: MODE=" << MODE << ", m=" << m << ", nl=" << nl << ", nr=" << nr << ", ld=" << ld << endl;
            if(ld > m || nl+nr < m)
                param_error();
        }
        int min(int a, int b) /// utility function for minima
        {
            return (a < b) ? a : b;            
        }
    // data
        int MODE;       /// modes: INIT, MATRIX, COEFFICIENTS
        int m;     /// filter order
        int nl;    /// points before
        int nr;    /// points after
        int ld;    /// derivative order
        
        LU<T>* lu;
        
        Array1D<T> b;
};

#endif
