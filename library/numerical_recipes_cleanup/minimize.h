#ifndef _MINIMIZE_H
#define _MINIMIZE_H

#include <iostream>

using namespace std;


// brent wants some things defined; we also use SHFT() to aid caching in fourier_mag()
template <class T>
void SHFT(T& a, T& b, T& c, const T& d)
    { a=b; b=c; c=d; }


#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define ITMAX 100
#define CGOLD 0.3819660
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define ZEPS 1.0e-10
/* Here ITMAX is the maximum allowed number of iterations; CGOLD is the golden ratio; ZEP is a small number
that protects against trying to achieve fractional accuracy for a minimum that happens to be exactly zero. */


// brent's method for finding the peak of a fourier transform. last three params are just pass-through
double dbrent_four(double ax, double bx, double cx, double tol, double *xmin, const WORD* data, const unsigned long num_points, const double dt)
{
//cout << "Starting Brent with ax=" << ax << ", bx=" << bx << ", cx=" << cx << "\n";

    int iter;
    double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    double e=0.0;	// This will be the distance moved on the step before last.
    
    a=(ax < cx ? ax : cx);	// a and b must be in ascending order, but input abscissas need not be.
    b=(ax > cx ? ax : cx);
    x=w=v=bx;	// Initializations...
    
    // negative sign b/c we want to find a maximum, and brent calculates a minimum.
    fw=fv=fx=-1.0*fourier_mag(data,num_points,dt,x);
    
//cout << x << "\t" << fx << "\n";        
    
    
    for (iter=1;iter<=ITMAX;iter++) {	// main program loop
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
        if (fabs(x-xm) <= (tol2-0.5*(b-a))) {	// test for done here.
//cout << "Apparently, Brent thinks we're done already. tol1=" << tol1 << ", tol2=" << tol2 << "\n";            
            *xmin=x;
            return fx;
        }
        if (fabs(e) > tol1) {	// construct a trial parabolic fit.
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if (q > 0.0) p = -p;
            q=fabs(q);
            etemp=e;
            e=d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                d=CGOLD*(e=(x >= xm ? a-x : b-x));
            /* The above conditions determine the acceptability of the parabolic fit. Here we take the 
            golden section into the larger of the two segments. */
            else {
                d=p/q;	// Take the parabolic step.
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                    d=SIGN(tol1,xm-x);
            }
        } else {
            d = CGOLD*(e=(x >= xm ? a-x : b-x));
        }
        u=(fabs(d) >= tol1 ? x+d : x+ SIGN(tol1,d));
        // This is the one function evaluation per iteration
        fu=-1.0*fourier_mag(data,num_points,dt,u);
        
//cout << "Brent: " << u << "\t" << fu << "\n";        
        
        if (fu <= fx) {	// Now decide what to do with our function evaluation.
            if (u >= x) a=x; else b=x;
            SHFT(v,w,x,u)
            SHFT(fv,fw,fx,fu)
        } else {
            if (u < x) a=u; else b=u;
            if (fu <= fw || w == x) {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            } else if (fu <= fv || v== x || v == w) {
                v=u;
                fv=fu;
            }
        }
    }	// done with housekeeping. back for another iteration.
    spit_error("Too many iterations in brent.");
    *xmin = x;
    return fx;
}

// Bracket our Fourier peak with lower and upper limits. Again, last three params pass through to fourier_mag()
// From Numerical Recipes: (partly out of date considering modifications) Given a function func, and given distinct initial points ax and bx, this routine searches in the downhill direction (defined by the function as evaluated at the initial points) and returns new points ax, bx, cx that bracket a minimum of the function. Also returned are the function values at the three points fa, fb, and fc.
void dpeakbracket_four(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, const WORD* data, const unsigned long num_points, const double dt)
{
    double ulim,u,r,q,fu,dum;
    
    *fa=-1.0*fourier_mag(data,num_points,dt,*ax);
    *fb=-1.0*fourier_mag(data,num_points,dt,*bx);
    if (*fb > *fa) {	// Switch roles of a and b so that we can go downhill in the direction from a to b
        SHFT(dum,*ax,*bx,dum)
        SHFT(dum,*fb,*fa,dum)
    }
    *cx=(*bx)+GOLD*(*bx-*ax);	// First guess for c
    *fc=-1.0*fourier_mag(data,num_points,dt,*cx);
    while (*fb > *fc) {	// Keep returning here until we bracket
        r=(*bx-*ax)*(*fb-*fc);	// Compute u by parabolic extrapolation from a,b,c.
        q=(*bx-*cx)*(*fb-*fa);
        u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(DMAX(fabs(q-r),TINY),q-r));
        ulim=(*bx)+GLIMIT*(*cx-*bx);
        
        // We won't go farther than this. Test various possibilities:
        if((*bx-u)*(u-*cx) > 0.0) {	// Parabolic u is between b and c: try it.
            fu=-1.0*fourier_mag(data,num_points,dt,u);
            if(fu<*fc) {	// Got a minimum between b and c.
                *ax=(*bx);
                *bx=u;
                *fa=(*fb);
                *fb=fu;
                return;
            } else if (fu > *fb) { // Got a minimum between a and u
                *cx=u;
                *fc=fu;
                return;
            }
            u=(*cx)+GOLD*(*cx-*bx);	// Parabolic fit was no use. Use default magnification.
            fu=-1.0*fourier_mag(data,num_points,dt,u);
        } else if ((*cx-u)*(u-ulim) > 0.0) { // parabolic fit is between c and its allowed limit
            fu=-1.0*fourier_mag(data,num_points,dt,u);
            if (fu < *fc) {
                SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
                SHFT(*fb,*fc,fu,-1.0*fourier_mag(data,num_points,dt,u))
            }
        } else if ((u-ulim)*(ulim-*cx) >= 0.0) { // limit parabolic u to maximum allowed value.
            u=ulim;
            fu=-1.0*fourier_mag(data,num_points,dt,u);
        } else { // reject parabolic u, use default magnification
            u=(*cx)+GOLD*(*cx-*bx);
            fu=-1.0*fourier_mag(data,num_points,dt,u);
        }
        SHFT(*ax,*bx,*cx,u)	// eliminate oldest point and continue
        SHFT(*fa,*fb,*fc,fu)
    }
}



#endif // _MINIMIZE_H