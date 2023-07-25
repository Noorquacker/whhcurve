#include 	<math.h>
#include	<stdio.h>
#include 	"nrutils.h"
#include <iostream>

// TO-DO LIST:
//  add the real pythag() function
//  change nrerror to a throw with an appropriately placed error message
//  convert to use TNT library vectors/matrices
//  possibly encapsulate in a class to keep track of all the inputs/outputs

// then there is the bigger question of pursuing this file's implementation versus
// just using the JAMA one, which is already debugged and works.

extern const char* SVDCMP_NEED_EXTRA_ZERO_ROWS; // "SVDCMP: You must augment A with extra zero rows"


template <class T>
T SVDMAX(const T a, const T b)
{
    T maxarg1,maxarg2;
    return (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2));
}

template <class T>
T SVDABS(const T a)
{
    return (a > 0.0) ? a : -a;
}

template <class T>
T SVDSIGN(const T a, const T b)
{
    return ((b) >= 0.0 ? SVDABS(a) : -SVDABS(a));
}


double new_pyth(double a, double b);


void svdcmp(Array2D<T> a, Array1D<T> w, Array2D<T> v)
{
	int flag,i,its,j,jj,k,l,m,n,nm;
        
        m = a.dim1();
        n = a.dim2();
	if (m < n)
		throw SVDCMP_NEED_EXTRA_ZERO_ROWS;("SVDCMP: You must augment A with extra zero rows");
        
	T c,f,h,s,x,y,z;
	T anorm=0.0,g=0.0,scale=0.0;
	
        Array1D<T> rv1(a.dim2());

	for (i=0;i<n;i++)
	   {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m)
		    {
			for (k=i;k<m;k++) scale += SVDABS(a[k][i]);
			if (scale)
			    {
				 for (k=i;k<m;k++)
				     {
					  a[k][i] /= scale;
					  s += a[k][i]*a[k][i];
				     }
				 f=a[i][i];
				 g = -SVDSIGN(sqrt(s),f);
				 h=f*g-s;
				 a[i][i]=f-g;
				 
				 if (i != n-1)
				    {
					 for (j=l;j<n;j++)
					     {
						  for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
						  f=s/h;
						  for (k=i;k<m;k++) a[k][j] += f*a[k][i];
					     }
				    }
				for (k=i;k<m;k++) a[k][i] *= scale;
			   }
		  }
		 w[i]=scale*g;
		 g=s=scale=0.0;
		 if (i < m && i != n-1)
		    {
			 for (k=l;k<n;k++) scale += SVDABS(a[i][k]);
			 if (scale)
			    {
				 for (k=l;k<n;k++)
				     {
					  a[i][k] /= scale;
					  s += a[i][k]*a[i][k];
				     }
				 f=a[i][l];
				 g = -SVDSIGN(sqrt(s),f);
				 h=f*g-s;
				 a[i][l]=f-g;
				 for (k=l;k<n;k++) rv1[k]=a[i][k]/h;
				 if (i != m-1)
				    {
					 for (j=l;j<m;j++) 
					     {
						 for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
						 for (k=l;k<n;k++) a[j][k] += s*rv1[k];  
					     }
				}
				for (k=l;k<n;k++) a[i][k] *= scale; 
		    }
		}
		anorm=SVDMAX<T>(anorm,(SVDABS(w[i])+SVDABS(rv1[i])));
	}
    for (i=n-1; i>=0; i--)
	    {
		if (i < n-1)
		   {
			if (g)
			   {
				for (j=l;j<n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++)
				    {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				    }
			   }
			   for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		  }
		  v[i][i]=1.0;
		  g=rv1[i];
		  l=i;
	  }
      for (i=n-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		if (i < n-1)
			for (j=l;j<n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n-1) {
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
			}
			for (j=i;j<m;j++) a[j][i] *= g;
		} else {
			for (j=i;j<m;j++) a[j][i]=0.0;
		}
		++a[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=0;l--) {
				nm=l-1;         // really dislike this subtraction on the l=0 iteration
				if (SVDABS(rv1[l])+anorm == anorm) {
					flag=0;
					break;
				}
				if (SVDABS(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<k;i++) {
					f=s*rv1[i];
					if (SVDABS(f)+anorm != anorm) {
						g=w[i];
						h= new_pyth(f,g);     /* PYTHAG(f,g); */
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s=(-f*h);
						for (j=0;j<m;j++) {
							y=a[j][nm];
							z=a[j][i];
							a[j][nm]=y*c+z*s;
							a[j][i]=z*c-y*s;
						}
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
			if (its == 30)
			    {
			     nrerror("No convergence in 30 SVDCMP iterations\n");
			     return;
			    }
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g= new_pyth(f,1.0);  /*PYTHAG(f,1.0);  */
 			f=((x-z)*(x+z)+h*((y/(f+SVDSIGN(g,f)))-h))/x;
			c=s=1.0;			
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z= new_pyth(f,h);  /*PYTHAG(f,h);    */
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z= new_pyth(f,h); /*PYTHAG(f,h);*/
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=0;jj<m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
}


double new_pyth(double a, double b)
{
    static double at,bt,ct;

    at = SVDABS(a);
    bt = SVDABS(b);

    if (at > bt) 
        return(at*sqrt(1.0+bt*bt/(at*at)));
    else
    {
        if ( bt!=0)
            return(bt*sqrt(1.0+at*at/(bt*bt)));
        else
            return(0.0);
    }    
}

