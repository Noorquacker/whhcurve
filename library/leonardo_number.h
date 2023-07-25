#ifndef _LEONARDO_NUMBER_H
#define _LEONARDO_NUMBER_H

// Leonardo numbers are used by the sorting algorithm Smoothsort, which is
// one of the best if the data are already nearly in order. In the case of 
// splining data, our arrays are usually already sorted either forwards
// or backwards, and Smoothsort can be used optimally. This implementation
// of Smoothsort is most immediately credited to WikiBooks:
// http://en.wikibooks.org/wiki/Algorithm_Implementation/Sorting/Smoothsort

class LeonardoNumber

/**
 **  Helper class for manipulation of Leonardo numbers
 **
 **/

{
public:
	/**  Default ctor.  **/
	LeonardoNumber (void) throw () : b (1), c (1)
	{ return; }
	
	/**  Copy ctor.  **/
	LeonardoNumber (const LeonardoNumber & _l) throw () : b (_l.b), c (_l.c)
	{ return; }
	
	/**  
	 **  Return the "gap" between the actual Leonardo number and the
	 **  preceding one.
	 **/
	int gap (void) const throw ()
	{ return b - c; }
	
	
	/**  Perform an "up" operation on the actual number.  **/
	LeonardoNumber & operator ++ (void) throw ()
	{ int s = b; b = b + c + 1; c = s; return * this; }
	
	/**  Perform a "down" operation on the actual number.  **/
	LeonardoNumber & operator -- (void) throw ()
	{ int s = c; c = b - c - 1; b = s; return * this; }
	
	/**  Return "companion" value.  **/
	int operator ~ (void) const throw ()
	{ return c; }
	
	/**  Return "actual" value.  **/
	operator int (void) const throw ()
	{ return b; }
	
	
private:
	int b;   /**  Actual number.  **/
	int c;   /**  Companion number.  **/
};


#endif