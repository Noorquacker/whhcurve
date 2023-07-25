#ifndef _SMOOTHSORT_H
#define _SMOOTHSORT_H

#include "tnt/tnt.h"
#include "leonardo_number.h"

using namespace TNT;


/**	SmoothSort
	
	Implements the Smoothsort algorithm using keys on TNT arrays.
	Usage:
	
	Array1D<T> primary = ...
	Array1D<U> secondary = ...
	
	SmoothSort<T> s(primary);
	s.reorder(secondary);
	
	With the first line, primary will be properly ordered, and s will contain the information to
	perform the same ordering transformation on any array of length primary.dim(). The second line
	uses this transformation information to reorder secondary. Note that primary and secondary
	need not be arrays of the same type, although they must be the same length. A max->min sort
	may be accomplished using the template parameter rsort:
	
	SmoothSort<T,true> s(primary);
	
	Since it is much faster to sort in the direction that is most closely sorted already, SmoothSort
	provides the ability to reverse the keys after sorting. If the array starts out nearly reverse-
	sorted, but ultimately must end up forward-sorted, the proper treatment is to sort it with
	rsort=true, then reverse it afterward:
	
	SmoothSort<T> s(primary);
	s.reverse();
	s.reorder(secondary);
	
*/
template <class T, bool rsort = false>
class SmoothSort
{
	public:
	SmoothSort(const Array1D<T>& comp_)
	{
		if(comp_.dim() == 0) throw "SmoothSort: Key array has zero size.";
		doSort(comp_);
	}
	
	void reverse()
	{
		for(int i=0; i<comp.dim()/2; i++)
		{
			swap(i,comp.dim()-1-i);
		}
	}
	
	
	// reorder the elements of a according to the sorting key and return
	// the result.
	template <class U>
	Array1D<U> reorder(const Array1D<U>& a) const
	{
		Array1D<U> b(a.dim());
		if(a.dim() != indices.dim()) throw "SmoothSort: Array to reorder not the same dimension as index array";
		for(int i=0; i<indices.dim(); i++)
		{
			b[i] = a[indices[i]];
		}
		return b;
	}

	// in-place reorder the elements of a and b according to the sorting key. function is
	// not provided for more than two arrays because the speed advantages are negligible and memory
	// usage higher.
	template <class U>
	void reorder(Array1D<U>& a, Array1D<U>& b) const
	{
		Array1D<U> aa = a.copy();
		Array1D<U> bb = b.copy();
		if(a.dim() != indices.dim()) throw "SmoothSort: Array to reorder (a) not the same dimension as index array";
		if(b.dim() != indices.dim()) throw "SmoothSort: Array to reorder (b) not the same dimension as index array";
		for(int i=0; i<indices.dim(); i++)
		{
			a[i] = aa[indices[i]];
			b[i] = bb[indices[i]];
		}
	}
	
	protected:
	Array1D<int> indices;
	Array1D<T> comp;
		
	void doSort(const Array1D<T>& comp_)
	{
		comp = comp_;
		
		indices = Array1D<int>(comp.dim());
		for(int i=0; i<comp.dim(); i++)
			indices[i] = i;
		
		smoothsort();
	}
	
	void swap(int i1, int i2)
	{
		T comptemp;
		int itemp;
		comptemp = comp[i1];
		itemp = indices[i1];
		comp[i1] = comp[i2];
		indices[i1] = indices[i2];
		comp[i2] = comptemp;
		indices[i2] = itemp;
	}
	
	bool ge(int i1, int i2)
	{
		return (rsort ? comp[i1] <= comp[i2] : comp[i1] >= comp[i2]);
	}
	
	
	// This section from http://en.wikibooks.org/wiki/Algorithm_Implementation/Sorting/Smoothsort
	
	
	/**
	 **  SmoothSort function template + helper functions.
	 **
	 **    Formal type T should have a comparison operator >= with prototype:
	 **
	 **      bool T::operator >= (const T &) const throw ();
	 **
	 **    which should compare its arguments by the given relation
	 **     (possibly taking advantage of the type itself).
	 **
	 **
	 **/
	
	void smoothsort()
	
	/**
	 **  Sorts the given array in ascending order.
	 **
	 **    Usage: smoothsort (<array>, <size>)
	 **
	 **    Where: <array> pointer to the first element of the array in question.
	 **            <size> length of the array to be sorted.
	 **
	 **
	 **/
	
    {
		int _n = comp.dim();
		
		int p = 1;
		LeonardoNumber b;
		
		for (int q = 0; ++q < _n ; ++p)
			if (p % 8 == 3)
			{
				sift(q - 1, b);
				
				++++b; p >>= 2;
			}
		
			else if (p % 4 == 1)
			{
				if (q + ~b < _n)  sift(q - 1, b);
				else  trinkle(q - 1, p, b);
				
				for (p <<= 1; --b > 1; p <<= 1)  ;
			}
		
		trinkle(_n - 1, p, b);
		
		for (--p; _n-- > 1; --p)
			if (b == 1)
				for ( ; !(p % 2); p >>= 1)  ++b;
		
			else if (b >= 3)
			{
				if (p)  semitrinkle(_n - b.gap(), p, b);
				
				--b; p <<= 1; ++p;
				semitrinkle(_n - 1, p, b);
				--b; p <<= 1; ++p;
			}
    }
	
		void sift (int _r, LeonardoNumber _b)
		
		/**
		 **  Sifts up the root of the stretch in question.
		 **
		 **    Usage: sift (<array>, <root>, <number>)
		 **
		 **    Where:     <array> Pointer to the first element of the array in
		 **                       question.
		 **                <root> Index of the root of the array in question.
		 **              <number> Current Leonardo number.
		 **
		 **
		 **/
		
        {
			int r2;
			
			while (_b >= 3)
            {
				if (ge(_r - _b.gap(), _r - 1))
					r2 = _r - _b.gap();
				else
                { r2 = _r - 1; --_b; }
				
				if (ge(_r, r2))  break;
				else
                { swap(_r, r2); _r = r2; --_b; }
            }
			
			
			return;
        }
		
		
		void semitrinkle (int _r, int long long _p, LeonardoNumber _b)
		
		/**
		 **  Trinkles the roots of the stretches of a given array and root when the
		 **  adjacent stretches are trusty.
		 **
		 **    Usage: semitrinkle (<array>, <root>, <standard_concat>, <number>)
		 **
		 **    Where:           <array> Pointer to the first element of the array in
		 **                             question.
		 **                      <root> Index of the root of the array in question.
		 **           <standard_concat> Standard concatenation's codification.
		 **                    <number> Current Leonardo number.
		 **
		 **
		 **/
		
        {
			if ( ge(_r - ~_b, _r))
            {
				swap(_r, _r - ~_b);
				trinkle(_r - ~_b, _p, _b);
            }
        }
		
		
		void trinkle(int _r, int _p, LeonardoNumber _b)
		
		/**
		 **  Trinkles the roots of the stretches of a given array and root.
		 **
		 **    Usage: trinkle (<array>, <root>, <standard_concat>, <number>)
		 **
		 **    Where:           <array> Pointer to the first element of the array in
		 **                             question.
		 **                      <root> Index of the root of the array in question.
		 **           <standard_concat> Standard concatenation's codification.
		 **                    <number> Current Leonardo number.
		 **
		 **
		 **/
		
        {
			while (_p)
            {
				for ( ; !(_p % 2); _p >>= 1)  ++_b;
				
				if (!--_p || (ge(_r, _r - _b)))  break;
				else
					if (_b == 1)
					{ swap(_r, _r - _b); _r -= _b; }
				
					else if (_b >= 3)
					{
						int r2 = _r - _b.gap (), r3 = _r - _b;
						
						if( ge(_r - 1, r2))
						{ r2 = _r - 1; _p <<= 1; --_b; }
						
						if ( ge(r3, r2) )
						{ swap(_r, r3); _r = r3; }
						
						else
						{ swap(_r, r2); _r = r2; --_b; break; }
					}
            }
			
			sift(_r, _b);
        }
	
};

#endif
