#include <stdio.h>
#include <math.h>
#include <immintrin.h>



#ifndef _baibai_h_
#define _baibai_h_


struct dd_real {
  double hi;
  double lo;


  inline dd_real operator=(double b) {
	  this->hi = b;
	  this->lo = 0.0;
	  return *this;
  }
};



struct double_comp {
	double re;
	double im;

	inline double_comp operator=(double b) {
		this->re = b;
		this->im = 0.0;
		return *this;
	}


};


struct dd_comp {
  dd_real re;
  dd_real im;

  inline dd_comp operator=(double b) {
	  this->re.hi = b;
	  this->re.lo = 0.0;
	  this->im.hi = 0.0;
	  this->im.lo = 0.0;
	  return *this;
  }

  inline dd_comp operator=(double_comp b) {
	  this->re.hi = b.re;
	  this->re.lo = 0.0;
	  this->im.hi = b.im;
	  this->im.lo = 0.0;
	  return *this;
  }

  inline dd_comp operator=(dd_real b) {
	  this->re = b;
	  this->im.hi = 0.0;
	  this->im.lo = 0.0;
	  return *this;
  }

};



inline dd_comp zero_comp(void);

inline dd_real twoproduct(double a,double b){
  dd_real seki;
  seki.hi=a*b;
  seki.lo=fma(a,b,-seki.hi);
    return seki;
  
}
  
inline dd_real twosum(double a,double b){
  dd_real add;
  double z;
  add.hi=a+b;
  z=add.hi-a;
  add.lo=(a-(add.hi-z))+(b-z);

  return add;
}

inline dd_real twodiff(double a, double b) {
	dd_real dif;
	double z;
	dif.hi = a - b;
	z = dif.hi - a;
	dif.lo = (a - (dif.hi - z)) - (b + z);

	return dif;
}

inline dd_real fasttwosum(double a,double b){
  dd_real add;
  add.hi=a+b;
  add.lo=b-(add.hi-a);

  return add;
}

inline dd_comp fasttwosum_comp(dd_comp a) {

	a.re = fasttwosum(a.re.hi, a.re.lo);
	a.im = fasttwosum(a.im.hi, a.im.lo);

	return a;
}
inline dd_real cray_add(dd_real x,dd_real y){
  dd_real z;
  z=twosum(x.hi,y.hi);
  z.lo=z.lo+x.lo+y.lo;
  z=fasttwosum(z.hi,z.lo);
  
  return z;
}

inline dd_real cray_add_dd_d(dd_real x, double y) {
	dd_real z;
	z = twosum(x.hi, y);
	z.lo += x.lo;
//	z = fasttwosum(z.hi, z.lo);

	return z;
}

inline  dd_real cray_sub(dd_real x,dd_real y){
  dd_real z;
  z=twodiff(x.hi,y.hi);
  z.lo=z.lo+x.lo-y.lo;
  z=fasttwosum(z.hi,z.lo);
  
  return z;
}

inline  dd_real cray_sub_dd_d(dd_real x, double y) {
	dd_real z;
	z = twodiff(x.hi, y);
	z.lo += x.lo;
	z = fasttwosum(z.hi, z.lo);

	return z;
}

inline  dd_real cray_sub_d_dd(double x, dd_real y) {
	dd_real z;
	z = twodiff(x, y.hi);
	z.lo -= y.lo;
	z = fasttwosum(z.hi, z.lo);

	return z;
}

inline dd_real mul(dd_real x,dd_real y){
  dd_real z;
  
  z=twoproduct(x.hi,y.hi);
  z.lo=fma(x.hi,y.lo,z.lo);
  z.lo=fma(x.lo,y.hi,z.lo);
  z=fasttwosum(z.hi,z.lo);
  return z;
  
}


inline dd_real mul_dd_d(dd_real x,double y){
  dd_real z;

  z.hi=x.hi*y;
  z.lo=fma(x.hi,y,-z.hi);
  z.lo=fma(x.lo,y,z.lo);
  z=fasttwosum(z.hi,z.lo);
  
 return z;

}


extern inline dd_real fast_div(dd_real x,dd_real y){
  dd_real z,t;
  double y_r;
  if(x.hi==0.0 && x.lo==0.0)
    return x;


  y_r=1.0/y.hi;
  z.hi=x.hi*y_r;
  t=twoproduct(z.hi,y.hi);
  z.lo=((x.hi-t.hi)-t.lo)*y_r+z.hi*(x.lo/x.hi-y.lo*y_r);
  z=fasttwosum(z.hi,z.lo);

  return z;
}

extern inline double short_div(dd_real x, dd_real y) {
	dd_real z, t;
	double y_r;
	if (x.hi == 0.0 && x.lo == 0.0)
		return x.hi;


	y_r = 1.0 / y.hi;
	z.hi = x.hi*y_r;
	t = twoproduct(z.hi, y.hi);
	z.lo = ((x.hi - t.hi) - t.lo)*y_r + z.hi*(x.lo / x.hi - y.lo*y_r);
	z.hi += z.lo;

	return z.hi;
}


inline dd_real fast_div_dd_d(dd_real x,double y){
  dd_real z,t;
  double y_r;
  if(x.hi==0.0 && x.lo==0.0)
    return x;

  y_r=1.0/y;
  z.hi=x.hi*y_r;
  t=twoproduct(z.hi,y);
  z.lo=fma(z.hi,x.lo/x.hi,((x.hi-t.hi)-t.lo)*y_r);
  z=fasttwosum(z.hi,z.lo);

  return z;
}

extern inline dd_real fast_div_d_dd(double x, dd_real y) {
	dd_real z, t;
	double y_r;


	y_r = 1.0 / y.hi;
	z.hi = x * y_r;
	t = twoproduct(z.hi, y.hi);
	z.lo = ((x - t.hi) - t.lo)*y_r - z.hi * y.lo * y_r;
	z = fasttwosum(z.hi, z.lo);

	return z;
}

extern inline dd_comp my_plus(dd_comp a,dd_comp b){
  dd_comp c;
   
  c.re=cray_add(a.re,b.re);
  c.im=cray_add(a.im,b.im);

  return c;
}

inline dd_comp my_short_plus(dd_comp a,dd_comp b){
  dd_comp c;
   
  c.re.hi=a.re.hi+b.re.hi;
  c.im.hi=a.im.hi+b.im.hi;

  c.re.lo=0.0;
  c.im.lo=0.0;
  
  return c;
}

extern inline dd_comp my_plus_d_d(double_comp a, double_comp b) {
	dd_comp c;

	c.re = twosum(a.re, b.re);
	c.im = twosum(a.im, b.im);

	return c;
}

extern inline dd_comp my_plus_dd_d(dd_comp a, double_comp b) {
	dd_comp c;

	c.re = cray_add_dd_d(a.re, b.re);
	c.im = cray_add_dd_d(a.im, b.im);

	return c;
}

inline double_comp  my_d_plus(double_comp a, double_comp b) {
	double_comp c;
	c.re = a.re + b.re;
	c.im = a.im + b.im;

	return c;
}

extern inline dd_comp my_minu(dd_comp a,dd_comp b){
  dd_comp c;
   
  c.re=cray_sub(a.re,b.re);
  c.im=cray_sub(a.im,b.im);
  
  return c;
}

extern inline dd_comp my_minu_dd_d(dd_comp a, double_comp b) {
	dd_comp c;

	c.re = cray_sub_dd_d(a.re, b.re);
	c.im = cray_sub_dd_d(a.im, b.im);

	return c;
}

extern inline dd_comp my_minu_d_dd(double_comp a, dd_comp b) {
	dd_comp c;

	c.re = cray_sub_d_dd(a.re, b.re);
	c.im = cray_sub_d_dd(a.im, b.im);

	return c;
}

extern inline dd_comp my_minu_d_d(double_comp a, double_comp b) {
	dd_comp c;

	c.re = twodiff(a.re, b.re);
	c.im = twodiff(a.im, b.im);

	return c;
}

inline dd_comp  my_short_minu(dd_comp a,dd_comp b){
  dd_comp c;
  c.re.hi=a.re.hi-b.re.hi;
  c.im.hi=a.im.hi-b.im.hi;

  c.re.lo=0.0;
  c.im.lo=0.0;
  return c;
}


inline double_comp  my_d_minu(double_comp a,double_comp b){
  double_comp c;
  c.re=a.re-b.re;
  c.im=a.im-b.im;

  return c;
}


inline dd_real cray_add_for_dot(dd_real x, dd_real y) {
	dd_real z;
	z = twosum(x.hi, y.hi);
	z.lo = z.lo + x.lo + y.lo;


	return z;
}

extern inline dd_comp my_plus_for_dot(dd_comp a, dd_comp b) {
	dd_comp c;

	c.re = cray_add_for_dot(a.re, b.re);
	c.im = cray_add_for_dot(a.im, b.im);

	return c;
}

extern inline dd_comp  my_mult_for_dot(dd_comp a, dd_comp b) {
	dd_comp c;
	dd_real t1, t2;
	double d1, d2, d3, d4;
	//c.re=a.re*b.re - a.im*b.im;


	d1 = a.re.hi * b.re.hi;
	d2 = fma(a.re.hi, b.re.hi, -d1);
	d2 = fma(a.re.hi, b.re.lo, d2);
	d2 = fma(a.re.lo, b.re.hi, d2);




	d3 = a.im.hi*b.im.hi;
	d4 = fma(a.im.hi, b.im.hi, -d3);
	d4 = fma(a.im.hi, b.im.lo, d4);
	d4 = fma(a.im.lo, b.im.hi, d4);

	c.re = twosum(d1, -d3);
	c.re.lo += d2 - d4;

	c.re = fasttwosum(c.re.hi, c.re.lo);

	//  c.im=a.re*b.im + a.im*b.re;

	d1 = a.re.hi * b.im.hi;
	d2 = fma(a.re.hi, b.im.hi, -d1);
	d2 = fma(a.re.hi, b.im.lo, d2);
	d2 = fma(a.re.lo, b.im.hi, d2);




	d3 = a.im.hi*b.re.hi;
	d4 = fma(a.im.hi, b.re.hi, -d3);
	d4 = fma(a.im.hi, b.re.lo, d4);
	d4 = fma(a.im.lo, b.re.hi, d4);

	c.im = twosum(d1, d3);
	c.im.lo += d2 + d4;



	return c;
}


extern inline dd_comp  my_mult(dd_comp a, dd_comp b) {
	dd_comp c;
	dd_real t1, t2;
	double d1, d2, d3, d4;
	//c.re=a.re*b.re - a.im*b.im;


	d1 = a.re.hi * b.re.hi;
	d2 = fma(a.re.hi, b.re.hi, -d1);
	d2 = fma(a.re.hi, b.re.lo, d2);
	d2 = fma(a.re.lo, b.re.hi, d2);




	d3 = a.im.hi*b.im.hi;
	d4 = fma(a.im.hi, b.im.hi, -d3);
	d4 = fma(a.im.hi, b.im.lo, d4);
	d4 = fma(a.im.lo, b.im.hi, d4);

	c.re = twosum(d1, -d3);
	c.re.lo += d2 - d4;

	c.re = fasttwosum(c.re.hi, c.re.lo);

	//  c.im=a.re*b.im + a.im*b.re;

	d1 = a.re.hi * b.im.hi;
	d2 = fma(a.re.hi, b.im.hi, -d1);
	d2 = fma(a.re.hi, b.im.lo, d2);
	d2 = fma(a.re.lo, b.im.hi, d2);




	d3 = a.im.hi*b.re.hi;
	d4 = fma(a.im.hi, b.re.hi, -d3);
	d4 = fma(a.im.hi, b.re.lo, d4);
	d4 = fma(a.im.lo, b.re.hi, d4);

	c.im = twosum(d1, d3);
	c.im.lo += d2 + d4;

	c.im = fasttwosum(c.im.hi, c.im.lo);


	return c;
}

inline dd_comp  my_short_mult(dd_comp a,dd_comp b){
  dd_comp c;

  c.re.hi=a.re.hi*b.re.hi - a.im.hi*b.im.hi;
  c.im.hi=a.re.hi*b.im.hi + a.im.hi*b.re.hi;

  c.re.lo=0.0;
  c.im.lo=0.0;
  
  return c;
}

inline dd_comp  my_short_mult_dd_d(dd_comp a,double_comp b){
  dd_comp c;

  c.re.hi=a.re.hi*b.re - a.im.hi*b.im;
  c.im.hi=a.re.hi*b.im + a.im.hi*b.re;

  c.re.lo=0.0;
  c.im.lo=0.0;
  
  return c;
}

inline double_comp  my_d_mult_dd_d(dd_comp a,double_comp b){
  double_comp c;
 

  c.re=a.re.hi*b.re - a.im.hi*b.im;
  c.im=a.re.hi*b.im + a.im.hi*b.re;

  return c;
}

inline double_comp  my_d_mult_d_dd(double_comp a,dd_comp b){
  double_comp c;
 

  c.re=a.re*b.re.hi - a.im*b.im.hi;
  c.im=a.re*b.im.hi + a.im*b.re.hi;

  return c;
}

inline dd_comp  my_short_mult_d_dd(double_comp a,dd_comp b){
  dd_comp c;

  c.re.hi=a.re*b.re.hi - a.im*b.im.hi;
  c.im.hi=a.re*b.im.hi + a.im*b.re.hi;

  c.re.lo=0.0;
  c.im.lo=0.0;
  
  return c;
}


inline double_comp  my_d_mult(double_comp a,double_comp b){
  double_comp c;

  c.re=a.re*b.re - a.im*b.im;
  c.im=a.re*b.im + a.im*b.re;


  
  return c;
}



extern inline dd_comp  my_mult_dd_d(dd_comp a, double_comp b) {
	dd_comp c;
	dd_real t1, t2;
	double d1, d2, d3, d4;
	//c.re=a.re*b.re - a.im*b.im;


	d1 = a.re.hi * b.re;
	d2 = fma(a.re.hi, b.re, -d1);
	d2 = fma(a.re.lo, b.re, d2);




	d3 = a.im.hi*b.im;
	d4 = fma(a.im.hi, b.im, -d3);
	d4 = fma(a.im.lo, b.im, d4);

	c.re = twosum(d1, -d3);
	c.re.lo += d2 - d4;

	c.re = fasttwosum(c.re.hi, c.re.lo);

	//  c.im=a.re*b.im + a.im*b.re;

	d1 = a.re.hi * b.im;
	d2 = fma(a.re.hi, b.im, -d1);
	d2 = fma(a.re.lo, b.im, d2);




	d3 = a.im.hi*b.re;
	d4 = fma(a.im.hi, b.re, -d3);
	d4 = fma(a.im.lo, b.re, d4);

	c.im = twosum(d1, d3);
	c.im.lo += d2 + d4;

	c.im = fasttwosum(c.im.hi, c.im.lo);


	return c;
}


inline dd_comp  my_mult_d_dd(double_comp a, dd_comp b) {
	dd_comp c;
	dd_real t1, t2;
	double d1, d2, d3, d4;
	//c.re=a.re*b.re - a.im*b.im;


	d1 = a.re * b.re.hi;
	d2 = fma(a.re, b.re.hi, -d1);
	d2 = fma(a.re, b.re.lo, d2);




	d3 = a.im*b.im.hi;
	d4 = fma(a.im, b.im.hi, -d3);
	d4 = fma(a.im, b.im.lo, d4);

	c.re = twosum(d1, -d3);
	c.re.lo += d2 - d4;

	c.re = fasttwosum(c.re.hi, c.re.lo);

	//  c.im=a.re*b.im + a.im*b.re;

	d1 = a.re * b.im.hi;
	d2 = fma(a.re, b.im.hi, -d1);
	d2 = fma(a.re, b.im.lo, d2);




	d3 = a.im*b.re.hi;
	d4 = fma(a.im, b.re.hi, -d3);
	d4 = fma(a.im, b.re.lo, d4);

	c.im = twosum(d1, d3);
	c.im.lo += d2 + d4;

	c.im = fasttwosum(c.im.hi, c.im.lo);


	return c;
}

inline dd_comp  my_mult_d_d(double_comp a, double_comp b) {
	dd_comp c;
	dd_real t1, t2;
	double d1, d2, d3, d4;
	//c.re=a.re*b.re - a.im*b.im;


	d1 = a.re * b.re;
	d2 = fma(a.re, b.re, -d1);


	d3 = a.im*b.im;
	d4 = fma(a.im, b.im, -d3);

	c.re = twosum(d1, -d3);
	c.re.lo += d2 - d4;

	c.re = fasttwosum(c.re.hi, c.re.lo);

	//  c.im=a.re*b.im + a.im*b.re;

	d1 = a.re * b.im;
	d2 = fma(a.re, b.im, -d1);




	d3 = a.im*b.re;
	d4 = fma(a.im, b.re, -d3);

	c.im = twosum(d1, d3);
	c.im.lo += d2 + d4;

	c.im = fasttwosum(c.im.hi, c.im.lo);


	return c;
}

inline dd_comp  my_mult_dd_ddreal(dd_comp a, dd_real b) {
	dd_comp c;
	dd_real t1, t2;

	//c.re=a.re*b.re - a.im*b.im;
	c.re = mul(a.re, b);
	c.im = mul(a.im, b);

	return c;
}

inline dd_comp  my_mult_dd_double(dd_comp a,double b){
  dd_comp c;
  dd_real t1,t2;

  //c.re=a.re*b.re - a.im*b.im;
  c.re=mul_dd_d(a.re,b);  
  c.im=mul_dd_d(a.im,b);

   return c;
}


inline dd_comp  my_mult_d_ddreal(double_comp a, dd_real b) {
	dd_comp c;

	//c.re=a.re*b.re - a.im*b.im;
	c.re = mul_dd_d(b, a.re);
	c.im = mul_dd_d(b, a.im);

	return c;
}

inline dd_comp  my_mult_d_double(double_comp a, double b) {
	dd_comp c;

	//c.re=a.re*b.re - a.im*b.im;
	c.re = twoproduct(a.re, b);
	c.im = twoproduct(a.im, b);

	return c;
}

extern inline dd_comp my_divv(dd_comp a, dd_comp b) {
	dd_comp c;
	dd_real t1, t2, t3, t4;
	double d1, d2, d3, d4;
	//  c.re=(a.re*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im);

	if (b.im.hi == 0.0) {
		//  c.re=a.re/b.re;
		c.re = fast_div(a.re, b.re);

		//  c.im=a.im/b.re;
		c.im = fast_div(a.im, b.re);

		return c;

	}
	else {

		//  c.re=(a.re*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im);
		d1 = b.re.hi * b.re.hi;
		d2 = fma(b.re.hi, b.re.hi, -d1);
		d2 = fma(b.re.hi*2.0, b.re.lo, d2);


		d3 = b.im.hi * b.im.hi;
		d4 = fma(b.im.hi, b.im.hi, -d3);
		d4 = fma(b.im.hi*2.0, b.im.lo, d4);

		t4 = twosum(d1, d3);
		t4.lo += d2 + d4;


		d1 = a.re.hi * b.re.hi;
		d2 = fma(a.re.hi, b.re.hi, -d1);
		d2 = fma(a.re.hi, b.re.lo, d2);
		d2 = fma(a.re.lo, b.re.hi, d2);

		d3 = a.im.hi*b.im.hi;
		d4 = fma(a.im.hi, b.im.hi, -d3);
		d4 = fma(a.im.hi, b.im.lo, d4);
		d4 = fma(a.im.lo, b.im.hi, d4);


		c.re = twosum(d1, d3);
		c.re.lo += d2 + d4;


		c.re = fast_div(c.re, t4);


		//  c.im=(a.im*b.re - a.re*b.im)/(b.re*b.re+b.im*b.im);
		d1 = a.im.hi * b.re.hi;
		d2 = fma(a.im.hi, b.re.hi, -d1);
		d2 = fma(a.im.hi, b.re.lo, d2);
		d2 = fma(a.im.lo, b.re.hi, d2);

		d3 = a.re.hi*b.im.hi;
		d4 = fma(a.re.hi, b.im.hi, -d3);
		d4 = fma(a.re.hi, b.im.lo, d4);
		d4 = fma(a.re.lo, b.im.hi, d4);


		c.im = twosum(d1, -d3);
		c.im.lo += d2 - d4;
		c.im = fast_div(c.im, t4);

		return c;
	}
}

extern inline double_comp my_d_divv_dd_dd(dd_comp a, dd_comp b) {
	double_comp c;
	dd_real t1, t2, t3, t4;
	double d1, d2, d3, d4;
	//  c.re=(a.re*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im);



		//  c.re=(a.re*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im);
		d1 = b.re.hi * b.re.hi;
		d2 = fma(b.re.hi, b.re.hi, -d1);
		d2 = fma(b.re.hi*2.0, b.re.lo, d2);


		d3 = b.im.hi * b.im.hi;
		d4 = fma(b.im.hi, b.im.hi, -d3);
		d4 = fma(b.im.hi*2.0, b.im.lo, d4);

		t4 = twosum(d1, d3);
		t4.lo += d2 + d4;


		d1 = a.re.hi * b.re.hi;
		d2 = fma(a.re.hi, b.re.hi, -d1);
		d2 = fma(a.re.hi, b.re.lo, d2);
		d2 = fma(a.re.lo, b.re.hi, d2);

		d3 = a.im.hi*b.im.hi;
		d4 = fma(a.im.hi, b.im.hi, -d3);
		d4 = fma(a.im.hi, b.im.lo, d4);
		d4 = fma(a.im.lo, b.im.hi, d4);


		t1 = twosum(d1, d3);
		t1.lo += d2 + d4;


		c.re = short_div(t1, t4);


		//  c.im=(a.im*b.re - a.re*b.im)/(b.re*b.re+b.im*b.im);
		d1 = a.im.hi * b.re.hi;
		d2 = fma(a.im.hi, b.re.hi, -d1);
		d2 = fma(a.im.hi, b.re.lo, d2);
		d2 = fma(a.im.lo, b.re.hi, d2);

		d3 = a.re.hi*b.im.hi;
		d4 = fma(a.re.hi, b.im.hi, -d3);
		d4 = fma(a.re.hi, b.im.lo, d4);
		d4 = fma(a.re.lo, b.im.hi, d4);


		t2 = twosum(d1, -d3);
		t2.lo += d2 - d4;
		c.im = short_div(t2, t4);

		return c;
	
}


inline dd_comp my_short_divv(dd_comp a,dd_comp b){
  dd_comp c;
  double t;
  t=b.re.hi*b.re.hi+b.im.hi*b.im.hi;
  
  c.re.hi=(a.re.hi*b.re.hi + a.im.hi*b.im.hi)/t;
  c.im.hi=(a.im.hi*b.re.hi - a.re.hi*b.im.hi)/t;

  c.re.lo=0.0;
  c.im.lo=0.0;
  
  return c;
}


inline dd_comp my_short_divv_dd_d(dd_comp a, double_comp b) {
	dd_comp c;
	double t;
	t = b.re*b.re + b.im*b.im;

	c.re.hi = (a.re.hi*b.re + a.im.hi*b.im) / t;
	c.im.hi = (a.im.hi*b.re - a.re.hi*b.im) / t;

	c.re.lo = 0.0;
	c.im.lo = 0.0;

	return c;
}

inline double_comp my_d_divv(double_comp a, double_comp b) {
	double_comp c;
	double t;
	t = b.re*b.re + b.im*b.im;

	c.re = (a.re*b.re + a.im*b.im) / t;
	c.im = (a.im*b.re - a.re*b.im) / t;


	return c;
}

extern inline double_comp my_acu_d_divv(double_comp a, double_comp b) {

	double_comp c;
	dd_real t3, t4, t5, t6;
	double d1, d2;

	//  c.re=(a.re*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im);

	t3 = twoproduct(a.re, b.re);
	t4 = twoproduct(a.im, b.im);
	t5 = cray_add(t3, t4);
	t3 = twoproduct(b.re, b.re);
	t4 = twoproduct(b.im, b.im);
	t6 = cray_add(t3, t4);
	c.re = short_div(t5, t6);

	t3 = twoproduct(a.im, b.re);
	t4 = twoproduct(a.re, b.im);
	t5 = cray_sub(t3, t4);

	c.im = short_div(t5, t6);

	return c;
}



extern inline dd_comp my_divv_d_d(double_comp a,double_comp b){

  dd_comp c;
  dd_real t3,t4,t5,t6;
  double d1,d2;
  
  //  c.re=(a.re*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im);

  t3=twoproduct(a.re,b.re);
  t4=twoproduct(a.im,b.im); 
  t5=cray_add(t3,t4);
  t3=twoproduct(b.re,b.re);  
  t4=twoproduct(b.im,b.im);
  t6=cray_add(t3,t4);  
  c.re=fast_div(t5,t6); 

  t3=twoproduct(a.im,b.re);  
  t4=twoproduct(a.re,b.im);   
  t5=cray_sub(t3,t4);

  c.im=fast_div(t5,t6);
  
  return c;
}



extern inline dd_comp my_divv_dd_d(dd_comp a, double_comp b) {
	dd_comp c;
	dd_real t1, t2, t3, t4;
	double d1, d2, d3, d4;
	//  c.re=(a.re*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im);

	if (b.im == 0.0) {
		//  c.re=a.re/b.re;
		c.re = fast_div_dd_d(a.re, b.re);

		//  c.im=a.im/b.re;
		c.im = fast_div_dd_d(a.im, b.re);

		return c;

	}
	else {

		//  c.re=(a.re*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im);
		d1 = b.re * b.re;
		d2 = fma(b.re, b.re, -d1);


		d3 = b.im * b.im;
		d4 = fma(b.im, b.im, -d3);

		t4 = twosum(d1, d3);
		t4.lo += d2 + d4;


		d1 = a.re.hi * b.re;
		d2 = fma(a.re.hi, b.re, -d1);
		d2 = fma(a.re.lo, b.re, d2);

		d3 = a.im.hi*b.im;
		d4 = fma(a.im.hi, b.im, -d3);
		d4 = fma(a.im.lo, b.im, d4);


		c.re = twosum(d1, d3);
		c.re.lo += d2 + d4;

		c.re = fast_div(c.re, t4);


		//  c.im=(a.im*b.re - a.re*b.im)/(b.re*b.re+b.im*b.im);
		d1 = a.im.hi * b.re;
		d2 = fma(a.im.hi, b.re, -d1);
		d2 = fma(a.im.lo, b.re, d2);

		d3 = a.re.hi*b.im;
		d4 = fma(a.re.hi, b.im, -d3);
		d4 = fma(a.re.lo, b.im, d4);


		c.im = twosum(d1, -d3);
		c.im.lo += d2 - d4;

		c.im = fast_div(c.im, t4);

		return c;
	}
}


extern inline dd_comp my_divv_d_dd(double_comp a, dd_comp b) {
	dd_comp c;
	dd_real t1, t2, t3, t4;
	double d1, d2, d3, d4;
	//  c.re=(a.re*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im);

	if (b.im.hi == 0.0) {
		//  c.re=a.re/b.re;
		c.re = fast_div_d_dd(a.re, b.re);

		//  c.im=a.im/b.re;
		c.im = fast_div_d_dd(a.im, b.re);

		return c;

	}
	else {

		//  c.re=(a.re*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im);
		d1 = b.re.hi * b.re.hi;
		d2 = fma(b.re.hi, b.re.hi, -d1);
		d2 = fma(b.re.hi*2.0, b.re.lo, d2);


		d3 = b.im.hi * b.im.hi;
		d4 = fma(b.im.hi, b.im.hi, -d3);
		d4 = fma(b.im.hi*2.0, b.im.lo, d4);

		t4 = twosum(d1, d3);
		t4.lo += d2 + d4;


		d1 = a.re * b.re.hi;
		d2 = fma(a.re, b.re.hi, -d1);
		d2 = fma(a.re, b.re.lo, d2);

		d3 = a.im*b.im.hi;
		d4 = fma(a.im, b.im.hi, -d3);
		d4 = fma(a.im, b.im.lo, d4);


		c.re = twosum(d1, d3);
		c.re.lo += d2 + d4;


		c.re = fast_div(c.re, t4);


		//  c.im=(a.im*b.re - a.re*b.im)/(b.re*b.re+b.im*b.im);
		d1 = a.im * b.re.hi;
		d2 = fma(a.im, b.re.hi, -d1);
		d2 = fma(a.im, b.re.lo, d2);

		d3 = a.re*b.im.hi;
		d4 = fma(a.re, b.im.hi, -d3);
		d4 = fma(a.re, b.im.lo, d4);


		c.im = twosum(d1, -d3);
		c.im.lo += d2 - d4;
		c.im = fast_div(c.im, t4);

		return c;
	}
}


inline dd_comp my_divv_dd_double(dd_comp a,double b){

  dd_comp c;
  dd_real t3,t4,t5,t6;
  double d1,d2;
  
  //  c.re=(a.re)/(b);
  c.re=fast_div_dd_d(a.re,b); 
  
  //  c.im=(a.im)/(b);
  c.im=fast_div_dd_d(a.im,b);
  
  return c;
}


extern inline dd_comp my_divv_double_dd(double a, dd_comp b) {
	dd_comp c;
	dd_real t1, t2, t3, t4;
	double d1, d2, d3, d4;
	//  c.re=(a.re*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im);

	if (b.im.hi == 0.0) {
		//  c.re=a.re/b.re;
		c.re = fast_div_d_dd(a, b.re);

		//  c.im=a.im/b.re;
		c.im.hi=0.0;
		c.im.lo = 0.0;

		return c;

	}
	else {

		//  c.re=(a.re*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im);
		d1 = b.re.hi * b.re.hi;
		d2 = fma(b.re.hi, b.re.hi, -d1);
		d2 = fma(b.re.hi*2.0, b.re.lo, d2);


		d3 = b.im.hi * b.im.hi;
		d4 = fma(b.im.hi, b.im.hi, -d3);
		d4 = fma(b.im.hi*2.0, b.im.lo, d4);

		t4 = twosum(d1, d3);
		t4.lo += d2 + d4;


		c.re.hi = a * b.re.hi;
		c.re.lo = fma(a, b.re.hi, -d1);
		c.re.lo = fma(a, b.re.lo, d2);



		c.re = fasttwosum(c.re.hi, c.re.lo);

		c.re = fast_div(c.re, t4);


		//  c.im=(a.im*b.re - a.re*b.im)/(b.re*b.re+b.im*b.im);

		c.im.hi = a*b.im.hi;
		c.im.lo = fma(a, b.im.hi, -d3);
		c.im.lo = fma(a, b.im.lo, d4);


		c.im = fasttwosum(-c.im.hi, -c.im.lo);

		c.im = fast_div(c.im, t4);

		return c;
	}
}


inline dd_comp zero_comp(void){

  dd_comp a;

  a.re.hi=0.0;
  a.re.lo=0.0;
  a.im.hi=0.0;
  a.im.lo=0.0;
  
  return a;
}

inline double_comp zero_double_comp(void) {

	double_comp a;

	a.re = 0.0;
	a.im = 0.0;

	return a;
}

inline dd_comp one_comp(void){

  dd_comp a;

  a.re.hi=1.0;
  a.re.lo=0.0;
  a.im.hi=0.0;
  a.im.lo=0.0;
  
  return a;
}

inline double_comp one_double_comp(void) {

	double_comp a;

	a.re = 1.0;
	a.im = 0.0;

	return a;
}


inline dd_real zero_dd_real(void){
  dd_real a;

  a.hi=0.0;
  a.lo=0.0;
  
  return a;
}


extern inline dd_comp  my_dot(int n, dd_comp *a, dd_comp *b) {
	int i;
	dd_comp c;
	c = zero_comp();
	for (i = 0; i < n; i++) {
		c = my_plus_for_dot(c, my_mult_for_dot(a[i], b[i]));
		if (i % 10 == 0) {
			c.re = fasttwosum(c.re.hi, c.re.lo);
			c.im = fasttwosum(c.im.hi, c.im.lo);
		}
	}
	c.re = fasttwosum(c.re.hi, c.re.lo);
	c.im = fasttwosum(c.im.hi, c.im.lo);

	return c;
}



inline dd_real my_sqrt_dd_real(dd_real a){
  int i,n;
  dd_real x,two_dd;

  n=50;
  x=two_dd=zero_dd_real();


  two_dd.hi=2.0;

  x.hi=1.0;

  
  for(i=0;i<n;i++)
    x=cray_sub(x,fast_div(cray_sub(mul(x,x),a),mul(x,two_dd)));
  
  
  return x;
}


inline dd_comp my_sqrt(dd_comp a){
  int i,n;
  dd_comp x,two_dd;

  n=50;
  x=two_dd=zero_comp();


  two_dd.re.hi=2.0;

  x.im.hi=1.0;

  
  for(i=0;i<n;i++)
    x=my_minu(x,my_divv(my_minu(my_mult(x,x),a),my_mult(x,two_dd)));
  
  
  return x;
}



inline dd_comp my_short_sqrt(dd_comp a){
  int i,n;
  dd_comp x,two_dd;

  n=50;
  x=two_dd=zero_comp();


  two_dd.re.hi=2.0;

  x.im.hi=1.0;

  
  for(i=0;i<n;i++)
    x=my_short_minu(x,my_short_divv(my_short_minu(my_short_mult(x,x),a),my_short_mult(x,two_dd)));
  
  
  return x;
}




inline dd_real comp_abs(dd_comp a){

  dd_real x;

  x=cray_add(mul(a.re,a.re),mul(a.im,a.im));
  x=my_sqrt_dd_real(x);
  
  
  return x;
}


inline double shortd_abs(dd_comp a) {

	double x;

	x = a.re.hi* a.re.hi + a.im.hi * a.im.hi;
	x = sqrt(x);


	return x;
}

inline dd_comp my_short(dd_comp a){

  a.re.lo=0.0;
  a.im.lo=0.0;

  return a;
}




inline dd_comp my_conj(dd_comp a){

  dd_comp x;

  x.re=a.re;
  x.im.hi=-a.im.hi;
  x.im.lo=-a.im.lo;

  return x;
}


inline dd_real dtodd(double a) {

	dd_real x;

	x.hi = a;
	x.lo = 0.0;
	return x;

}


inline dd_comp dtodd_comp(double_comp a) {

	dd_comp x;

	x.re.hi = a.re;
	x.re.lo = 0.0;
	x.im.hi = a.im;
	x.im.lo = 0.0;
	return x;

}

inline double_comp ddtod_comp(dd_comp a) {

	double_comp x;

	x.re = a.re.hi;
	x.im = a.im.hi;
	return x;

}




//----------real addition---------------//
inline dd_real operator+(dd_real a, dd_real b) {
	return cray_add(a, b);
}

inline dd_real operator+(dd_real a, double b) {
	return cray_add_dd_d(a, b);
}

inline dd_real operator+(double a, dd_real b) {
	return cray_add_dd_d(b, a);
}


inline dd_real operator+=(dd_real &a, dd_real b) {

	a = cray_add(a, b);
}

inline dd_real operator+=(dd_real &a, double b) {

	a = cray_add_dd_d(a, b);
}

inline double operator+=(double &a, dd_real b) {
	a += b.hi;
}

//----------real addition---------------//



//----------real subtracition---------------//
inline dd_real operator-(dd_real &a, dd_real &b) {
	return cray_sub(a, b);
}

inline dd_real operator-(dd_real &a, double b) {
	return cray_sub_dd_d(a, b);
}

inline dd_real operator-(double a, dd_real &b) {
	return cray_sub_d_dd(a, b);
}

inline dd_real operator-=(dd_real &a, dd_real &b) {

	a = cray_sub(a, b);
}

inline dd_real operator-=(dd_real &a, double b) {

	a = cray_sub_dd_d(a, b);
}

//----------real subtracition---------------//


//----------real multiplication---------------//
inline dd_real operator*(dd_real a, dd_real b) {
	return mul(a, b);
}

inline dd_real operator*(dd_real a, double b) {
	return mul_dd_d(a, b);
}

inline dd_real operator*(double a, dd_real b) {
	return mul_dd_d(b, a);
}

inline dd_real operator*=(dd_real &a, dd_real b) {

	a = mul(a, b);
}

inline dd_real operator*=(dd_real &a, double b) {

	a = mul_dd_d(a, b);
}

//----------real multiplication---------------//


//----------real division---------------//
inline dd_real operator/(dd_real a, dd_real b) {
	return fast_div(a, b);
}

inline dd_real operator/(dd_real a, double b) {
	return fast_div_dd_d(a, b);
}

inline dd_real operator/(double a, dd_real b) {
	return fast_div_d_dd(a, b);
}

inline dd_real operator/=(dd_real &a, dd_real b) {

	a = fast_div(a, b);
}

inline dd_real operator/=(dd_real &a, double b) {

	a = fast_div_dd_d(a, b);
}

//----------real division---------------//



//----------complex addition---------------//

inline dd_comp operator+(dd_comp a, dd_comp b) {
	return my_plus(a, b);
}

inline dd_comp operator+(dd_comp a, double_comp b) {
	return my_plus_dd_d(a, b);
}

inline dd_comp operator+(double_comp a, dd_comp b) {
	return my_plus_dd_d(b, a);
}

inline dd_comp operator+(double_comp a, double_comp b) {
	return my_plus_d_d(a, b);
}

inline dd_comp operator+=(dd_comp &a, dd_comp b) {
	a = my_plus(a, b);
}

inline dd_comp operator+=(dd_comp &a, double_comp b) {
	a = my_plus_dd_d(a, b);
}



//----------complex addition---------------//

//----------complex subtraction---------------//

inline dd_comp operator-(dd_comp a, dd_comp b) {
	return my_minu(a, b);
}

inline dd_comp operator-(dd_comp a, double_comp b) {
	return my_minu_dd_d(a, b);
}

inline dd_comp operator-(double_comp a, dd_comp b) {
	return my_minu_d_dd(a, b);
}

inline dd_comp operator-(double_comp a, double_comp b) {
	return my_minu_d_d(a, b);
}

inline dd_comp operator-=(dd_comp &a, dd_comp b) {
	a = my_minu(a, b);
}

inline dd_comp operator-=(dd_comp &a, double_comp b) {
	a = my_minu_dd_d(a, b);
}

//----------complex subtraction---------------//

//----------complex multiplication---------------//(double, double_comp, dd_real, dd_comp)‚Ì4‚Â‚ª‚ ‚é‚©‚ç16’Ê‚è

inline dd_comp operator*(dd_comp a, dd_comp b) {
	return my_mult(a, b);
}

inline dd_comp operator*(dd_comp a, double_comp b) {
	return my_mult_dd_d(a, b);
}

inline dd_comp operator*(dd_comp a, dd_real b) {
	return my_mult_dd_ddreal(a, b);
}

inline dd_comp operator*(dd_comp a, double b) {
	return my_mult_dd_double(a, b);
}



inline dd_comp operator*(double_comp a, dd_comp b) {
	return my_mult_dd_d(b, a);
}

inline dd_comp operator*(double_comp a, double_comp b) {
	return my_mult_d_d(a, b);
}

inline dd_comp operator*(double_comp a, dd_real b) {
	return my_mult_d_ddreal(a, b);
}

inline dd_comp operator*(double_comp a, double b) {
	return my_mult_d_double(a, b);
}



inline dd_comp operator*(dd_real a, dd_comp b) {
	return my_mult_dd_ddreal(b, a);
}

inline dd_comp operator*(dd_real a, double_comp b) {
	return my_mult_d_ddreal(b, a);
}



inline dd_comp operator*(double a, dd_comp b) {
	return my_mult_dd_double(b, a);
}

inline dd_comp operator*(double a, double_comp b) {
	return my_mult_d_double(b, a);
}


inline dd_comp operator*=(dd_comp &a, dd_comp b) {
	a = my_mult(a, b);
}

inline dd_comp operator*=(dd_comp &a, double_comp b) {
	a = my_mult_dd_d(a, b);
}




//----------complex multiplication---------------//

//----------complex division---------------//

inline dd_comp operator/(dd_comp a, dd_comp b) {
	return my_divv(a, b);
}

inline dd_comp operator/(dd_comp a, double_comp b) {
	return my_divv_dd_d(a, b);
}

inline dd_comp operator/(double_comp a, dd_comp b) {
	return my_divv_d_dd(a, b);
}

inline dd_comp operator/(double_comp a, double_comp b) {
	return my_divv_d_d(a, b);
}

inline dd_comp operator/(double a, dd_comp b) {
	return my_divv_double_dd(a, b);
}

inline dd_comp operator/=(dd_comp &a, dd_comp b) {
	a = my_divv(a, b);
}

inline dd_comp operator/=(dd_comp &a, double_comp b) {
	a = my_divv_dd_d(a, b);
}

//----------complex division---------------//

#endif