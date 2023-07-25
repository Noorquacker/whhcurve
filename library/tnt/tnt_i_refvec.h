/*
*
* Template Numerical Toolkit (TNT)
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*/


#ifndef TNT_I_REFVEC_H
#define TNT_I_REFVEC_H

#include <cstdlib>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef TNT_BOUNDS_CHECK
#include <assert.h>
#endif

#ifndef NULL
#define NULL 0
#endif

namespace TNT
{
/*
	Internal representation of ref-counted array.  The TNT
	arrays all use this building block_.

	<p>
	If an array block_ is created by TNT, then every time 
	an assignment is made, the left-hand-side reference 
	is decreased by one, and the right-hand-side refernce
	count is increased by one.  If the array block_ was
	external to TNT, the refernce count is a NULL pointer
	regardless of how many references are made, since the 
	memory is not freed by TNT.


	
*/
template <class T>
class i_refvec
{


  private:
    T* data_;                  
    int *ref_count_;
	
	#ifdef _OPENMP
		omp_lock_t* lock_;
	#endif
	
	inline  void destroy();


  public:

			 i_refvec();
	explicit i_refvec(int n);
	inline	 i_refvec(T* data);
	inline	 i_refvec(const i_refvec &v);
	inline   T*		 begin();
	inline const T* begin() const;
	inline  T& operator[](int i);
	inline const T& operator[](int i) const;
	inline  i_refvec<T> & operator=(const i_refvec<T> &V);
		    void copy_(T* p, const T* q, const T* e); 
		    void set_(T* p, const T* b, const T* e); 
	inline 	int	 ref_count() const;
	inline  int is_null() const;
			 ~i_refvec();
			
};

template <class T>
void i_refvec<T>::copy_(T* p, const T* q, const T* e)
{
	for (T* t=p; q<e; t++, q++)
		*t= *q;
}

template <class T>
i_refvec<T>::i_refvec() : data_(NULL), ref_count_(NULL)
#ifdef _OPENMP
	,lock_(NULL)
#endif
{}

/**
	In case n is 0 or negative, it does NOT call new. 
*/
template <class T>
i_refvec<T>::i_refvec(int n) : data_(NULL), ref_count_(NULL)
#ifdef _OPENMP
	,lock_(NULL)
#endif	
{
	if (n >= 1)
	{
#ifdef TNT_DEBUG
		std::cout  << "new data storage.\n";
#endif
		data_ = new T[n];
		ref_count_ = new int;
		*ref_count_ = 1;

		#ifdef _OPENMP
			lock_ = new omp_lock_t;
			omp_init_lock(lock_);
		#endif
		
	}
}

template <class T>
inline	 i_refvec<T>::i_refvec(const i_refvec<T> &V): data_(V.data_), ref_count_(V.ref_count_)
#ifdef _OPENMP
	,lock_(V.lock_)
#endif
{
	#ifdef _OPENMP
		// must recast as a non-const because locking really
		// doesn't affect the const-ness of the underlying class
		#pragma omp flush
		if(V.lock_) omp_set_lock( (omp_lock_t*) V.lock_ );
	#endif
	if (V.ref_count_ != NULL)
	{
		(*(V.ref_count_))++;
	}
	#ifdef _OPENMP
		#pragma omp flush
		if(V.lock_) omp_unset_lock( (omp_lock_t*) V.lock_ );
	#endif
	
}


template <class T>
i_refvec<T>::i_refvec(T* data) : data_(data), ref_count_(NULL)
#ifdef _OPENMP
	,lock_(NULL)
#endif
{}

template <class T>
inline T* i_refvec<T>::begin()
{
	return data_;
}

template <class T>
inline const T& i_refvec<T>::operator[](int i) const
{
	return data_[i];
}

template <class T>
inline T& i_refvec<T>::operator[](int i)
{
	return data_[i];
}


template <class T>
inline const T* i_refvec<T>::begin() const
{
	return data_;
}



template <class T>
i_refvec<T> & i_refvec<T>::operator=(const i_refvec<T> &V)
{
	if (this == &V)
		return *this;
	
	#ifdef _OPENMP
		if(lock_) omp_set_lock(lock_);
		#pragma omp flush
	#endif
	if (ref_count_ != NULL)
	{
		
		(*ref_count_) --;
		if ((*ref_count_) == 0)
		{
			destroy();
		}
	}
	#ifdef _OPENMP
		#pragma omp flush
		if(lock_) omp_unset_lock(lock_);
	#endif
	
	#ifdef _OPENMP	// with _OPENMP
		#pragma omp flush
		if(V.lock_) omp_set_lock( (omp_lock_t*) V.lock_ );
		if(V.ref_count_ != NULL)
		{
			data_ = V.data_;
			ref_count_ = V.ref_count_;
			lock_ = V.lock_;
			(*(V.ref_count_))++;
			#pragma omp flush
		}
		else
		{
			data_ = NULL;
			ref_count_ = NULL;
		}
		if(V.lock_) omp_unset_lock( (omp_lock_t*) V.lock_ );
	#else	// without _OPENMP (no changes)
		data_ = V.data_;
		ref_count_ = V.ref_count_;
	
		if (V.ref_count_ != NULL)
			(*(V.ref_count_))++;
	#endif	// _OPENMP detection
	
	return *this;
}

template <class T>
void i_refvec<T>::destroy()
{
	if (ref_count_ != NULL)
	{

#ifdef TNT_DEBUG
		std::cout << "destorying data... \n";
#endif
		delete ref_count_;
		ref_count_ = NULL;

#ifdef TNT_DEBUG
		std::cout << "deleted ref_count_ ...\n";
#endif
		if (data_ != NULL)
			delete []data_;
#ifdef TNT_DEBUG
		std::cout << "deleted data_[] ...\n";
#endif
		data_ = NULL;
		
		#ifdef _OPENMP
//			omp_unset_lock(lock_);
			if(lock_) omp_destroy_lock(lock_);
			delete lock_;
			lock_ = NULL;
			#pragma omp flush
		#endif

	}
}

/*
* return 1 is vector is empty, 0 otherwise
*
* if is_null() is false and ref_count() is 0, then
* 
*/
template<class T>
int i_refvec<T>::is_null() const
{
	#ifdef _OPENMP
		#pragma omp flush
	#endif
	return (data_ == NULL ? 1 : 0);
}

/*
*  returns -1 if data is external, 
*  returns 0 if a is NULL array,
*  otherwise returns the positive number of vectors sharing
*  		this data space.
*/
template <class T>
int i_refvec<T>::ref_count() const
{
	#ifdef _OPENMP
		#pragma omp flush
	#endif
	if (data_ == NULL)
		return 0;
	else
		return (ref_count_ != NULL ? *ref_count_ : -1) ; 
}

template <class T>
i_refvec<T>::~i_refvec()
{
	#ifdef _OPENMP
		if(lock_) omp_set_lock(lock_);
		#pragma omp flush
	#endif	
	if (ref_count_ != NULL)
	{
		
		(*ref_count_)--;

		if (*ref_count_ == 0)
		{			
			destroy();
			return;
		}
		
	}
	#ifdef _OPENMP
		if(lock_) omp_unset_lock(lock_);
	#endif
}


} /* namespace TNT */


#endif
/* TNT_I_REFVEC_H */
