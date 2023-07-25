#ifndef _POLYFIT_H
#define _POLYFIT_H

#include "svdfit.h"
/**
Polynomial fitting:
<pre>
int main(int,char**)
{
    Array1D<double> x(20);
    Array1D<double> y(20);
    for(unsigned i=0; i<20; i++)
    {
        x[i] = 1.0*i;
        y[i] = 7.0*sin(i/30.0*2*3.14159);

        cout << x[i] << "\t" << y[i] << endl;        
    }
    
    PolyFit<double> f(20,3);
    try{
    f.fit(x,y);

    // parameters are f[0], f[1], f[2], etc...
    
    }catch(const char* foo) { cout << foo << endl; }
    
    return 0;
}
</pre>
*/
template <class T>
class PolyFit : public SVDFit<T>
{
    public:
    PolyFit(unsigned npmax, unsigned nfmax) { SVDFit<T>::init(npmax,nfmax,PolyFit<T>::polyphi); }
    
    static T polyphi(const unsigned j, const T x);

};

template <class T>
T PolyFit<T>::polyphi(const unsigned j, const T x)
{
    T q = 1.0;
    for(unsigned i=0; i<j; i++)
    {
        q *= x;        
    }
    
    return q;
}

#endif
