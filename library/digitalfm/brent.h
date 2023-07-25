#ifndef _BRENT_H
#define _BRENT_H

/** Sometimes we have to shift some variables around. */
template<class T>
void BrentSHFT(T& a, T& b, T& c, const T& d)
{
    a = b;
    b = c;
    c = d;
}

/** Sometimes we have to exchange variables. Just like SHFT, but
    has a built-in allocation call. */
template<class T>
void BrentEXCH(T& a, T& b)
{
    T foo = a;
    a = b;
    b = foo;    
}

template<class T>
T BrentSIGN(T a, T b)
{
    return b >= 0.0 ? BrentABS<T>(a) : -1.0*BrentABS<T>(a);
}

template<class T>
T BrentABS(T a)
{
    return a >= 0.0 ? a : -1.0*a;
}

template<class T>
T BrentM(this->ax)(T a, T b)
{
    return a > b ? a : b;
}


/** Virtual base class to perform a Brent's Method search for the peak of a function.
Tcoord is a type defining the inputs to the function. Tres defines the output type
of the function.

Tres MUST be a scalar: i.e., a Tres * Tcoord must be a Tcoord;

*/
template <class Tcoord, class Tres>
class Brent
{
    public:
    Brent() { is_valid = false; }
    /** I don't know why you would copy this, but... */
    Brent(const Brent& br)
        { ax=br.ax; bx=br.bx; cx=br.cx; xmin=br.xmin; fa=br.fa; fb=br.fb; fc=br.fc; fmin=br.fmin; is_valid=br.is_valid; }
        
    /** Use find_peak() to do actual Brent calculation */
    bool find_peak(Tcoord ina, Tcoord inb, double tol)
        { bracket_peak(ina,inb); return is_valid = find_bracketed_peak(tol); }
    /** Find out if data is valid */
    bool isvalid() { return is_valid; }
    void reset() { is_valid = false; }
    
    protected:
    void bracket_peak(Tcoord ina, Tcoord inb);
    bool find_bracketed_peak(double tol);
    Tres fn(Tcoord x) = 0;
    
    // constants
    static double CGOLD;
    static double GOLD;
    static double GLIMIT;
    static double TINY;
    static double ZEPS;
    
    /** bracketed coordinates */
    Tcoord ax, bx, cx, xmin;
    /** function evaluations at bracketed coordinates */
    Tres fa, fb, fc, fmin;

    bool is_valid;    
};

template <class Tcoord, class Tres>
Brent<Tcoord,Tres>::CGOLD = 0.3819660;

template <class Tcoord, class Tres>
Brent<Tcoord,Tres>::GOLD = 1.618034;

template <class Tcoord, class Tres>
Brent<Tcoord,Tres>::GLIMIT = 100.0;

template <class Tcoord, class Tres>
Brent<Tcoord,Tres>::TINY = 1e-20;

template <class Tcoord, class Tres>
Brent<Tcoord,Tres>::ZEPS = 1e-10;


// Bracket our peak with lower and upper limits. Again, last three params pass through to fourier_mag()
// From Numerical Recipes: (partly out of date considering modifications) Given a function func, and given distinct initial points (this->ax) and bx, this routine searches in the downhill direction (defined by the function as evaluated at the initial points) and returns new points (this->ax), bx, cx that bracket a minimum of the function. Also returned are the function values at the three points fa, fb, and fc.
template <class Tcoord, class Tres>
void Brent<Tcoord,Tres>::bracket_peak(T ina, T inb)
{
    Tcoord ulim,u,r,q,
    Tres fu;
    
    // start over
    if(is_valid)
        this->reset();
    
    this->ax = ina;
    this->bx = inb;
    
    (this->fa) = this->fn(this->ax);
    (this->fb) = this->fn(this->bx);
    if ((this->fb) > (this->fa)) {	// Switch roles of a and b so that we can go downhill in the direction from a to b
        BrentEXCH<Tcoord>((this->ax), this->bx);
        BrentEXCH<Tres>(this->fa, this->fb);
    }
    (this->cx) = (this->bx) + Brent::GOLD * ((this->bx)-((this->ax)));	// First guess for c
    (this->fc) = this->fn(this->cx);
    while ((this->fb) > (this->fc))
    {	// Keep returning here until we bracket
        r = ((this->bx)-(this->ax)) * ((this->fb)-(this->fc));	// Compute u by parabolic extrapolation from a,b,c.
        q = ((this->bx)-(this->cx)) * ((this->fb)-(this->fa));
        u = (this->bx)-(((this->bx) - (this->cx)) * q - ((this->bx)-(this->ax))*r)/(2.0*BrentSIGN<Tcoord>(BrentM(this->ax)<Tcoord>(BrentABS<Tcoord>(q-r),Brent::TINY),q-r));
        ulim = ((this->bx))+ Brent::GLIMIT * ((this->cx) - (this->bx));
        
        // We won't go farther than this. Test various possibilities:
        if(((this->bx)-u)*(u-(this->cx)) > 0.0) {	// Parabolic u is between b and c: try it.
            fu = this->fn(u);
            if(fu<(this->fc)) {	// Got a minimum between b and c.
                ((this->ax))=((this->bx));
                (this->bx)=u;
                (this->fa)=((this->fb));
                (this->fb)=fu;
                return;
            } else if (fu > (this->fb)) { // Got a minimum between a and u
                (this->cx)=u;
                (this->fc)=fu;
                return;
            }
            u=((this->cx))+Brent::GOLD*((this->cx)-(this->bx));	// Parabolic fit was no use. Use default magnification.
            fu=this->fn(u);
        } else if (((this->cx)-u)*(u-ulim) > 0.0) { // parabolic fit is between c and its allowed limit
            fu=this->fn(u);
            if (fu < (this->fc)) {
                BrentSHFT<Tcoord>(this->bx,this->cx,u,(this->cx)+GOLD*((this->cx)-(this->bx)));
                BrentSHFT<Tres>((this->fb),(this->fc),fu,this->fn(u))
            }
        } else if ((u-ulim)*(ulim-(this->cx)) >= 0.0) { // limit parabolic u to m(this->ax)imum allowed value.
            u=ulim;
            fu=this->fn(u);
        } else { // reject parabolic u, use default magnification
            u=((this->cx))+Brent::GOLD*((this->cx)-(this->bx));
            fu=this->fn(u);
        }
        BrentSHFT<Tcoord>(((this->ax)),(this->bx),(this->cx),u);	// eliminate oldest point and continue
        BrentSHFT<Tcoord>((this->fa),(this->fb),(this->fc),fu);
    }
}

template <class Tcoord, class Tres>
bool Brent<Tcoord,Tres>::find_bracketed_peak(double tol);
{
//cout << "Starting Brent with (this->ax)=" << (this->ax) << ", (this->bx)=" << (this->bx) << ", (this->cx)=" << (this->cx) << "\n";

    int iter;
    Tcoord a,b,d,p,q,r,u,x,w,v,xm;
    Tres fu,fv,fw,fx;
    Tcoord etemp;
    Tcoord e=0.0;	// This will be the distance moved on the step before last.
    double tol1,tol2;
    
    a=((this->ax) < (this->cx) ? (this->ax) : (this->cx));	// a and b must be in ascending order, but input abscissas need not be.
    b=((this->ax) > (this->cx) ? (this->ax) : (this->cx));
    x=w=v=(this->bx);	// Initializations...
    
    // negative sign b/c we want to find a m(this->ax)imum, and brent calculates a minimum.
    fw=fv=fx=this->fn(x);
    
//cout << x << "\t" << fx << "\n";        
    
    
    for (iter=1;iter<=ITMAX;iter++) {	// main program loop
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*BrentABS<Tcoord>(x)+Brent::ZEPS);
        if (BrentABS<Tcoord>(x-xm) <= (tol2-0.5*(b-a))) {	// test for done here.
//cout << "Apparently, Brent thinks we're done already. tol1=" << tol1 << ", tol2=" << tol2 << "\n";            
            this->xmin=x;
            this->fmin=fx;
            return true;
        }
        if (BrentABS<Tcoord>(e) > tol1) {	// construct a trial parabolic fit.
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if (q > 0.0) p = -p;
            q=BrentABS<Tcoord>(q);
            etemp=e;
            e=d;
            if (BrentABS<Tcoord>(p) >= BrentABS<Tcoord>(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                d=Brent::CGOLD*(e=(x >= xm ? a-x : b-x));
            /* The above conditions determine the acceptability of the parabolic fit. Here we take the 
            golden section into the larger of the two segments. */
            else {
                d=p/q;	// Take the parabolic step.
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                    d=BrentSIGN<Tcoord>(tol1,xm-x);
            }
        } else {
            d = Brent::CGOLD*(e=(x >= xm ? a-x : b-x));
        }
        u=(BrentABS<Tcoord>(d) >= tol1 ? x+d : x+ BrentSIGN<Tcoord>(tol1,d));
        // This is the one function evaluation per iteration
        fu=this->fn(u);
        
//cout << "Brent: " << u << "\t" << fu << "\n";        
        
        if (fu <= fx) {	// Now decide what to do with our function evaluation.
            if (u >= x) a=x; else b=x;
            BrentSHFT<Tcoord>(v,w,x,u)
            BrentSHFT<Tres>(fv,fw,fx,fu)
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
//    spit_error("Too many iterations in brent.");
    this->xmin = x;
    this->fmin = fx;
    return false
}

typedef Brent<double,double> DBrent;

#endif
