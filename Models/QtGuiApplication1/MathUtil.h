//	MathUtil.h: Mathクラス
//
//  	Copyright(C) 2001,  T. Ishii
//////////////////////////////////////////////////////////////////////

#if !defined(__MATH_UTIL_H_INCLUDED__)
#define __MATH_UTIL_H_INCLUDED__ 
#include "varray.h"
#include"globalFunc.h"
namespace base {
//一般的な数学メソッドをまとめたクラス
#if !defined(FASTCALL) 
#define FASTCALL _fastcall
#endif

namespace Math
{
	const float  Pi = 3.14159265358979323846f;
	const float  Pi2       = 6.28318530717958623200f; // 2 * Pi
	const float  PiDiv2   = 1.57079632679489655800f; // Pi / 2
	const float  PiDiv4   = 0.78539816339744827900f; // Pi / 4
	const float  InvPi     = 0.31830988618379069122f; // 1 / Pi
	const float  DegToRad   = 0.01745329251994329547f; // Degrees to Radians
	const float  RadToDeg   = 57.29577951308232286465f; // Radians to Degrees
	const float  Huge       = 1.0e+38f;                // Huge number for FLOAT
	const float  Epsilon    = 1.0e-6f;                 // Tolerance for FLOATs
	const float  EPS3       = 1.0e-3f;
	const float  EPS2       = 1.0e-2f;
	const float  EPS1       = 1.0e-1f;

	  float FASTCALL  Acos(float x);
	  float FASTCALL  Asin(float x);
	  float FASTCALL  Cos(float x);
	  float FASTCALL Sin(float x);
	  float FASTCALL Tan(float x);
	  float FASTCALL Atan(float x);
	  float FASTCALL Atan2(float x, float y);
	  float FASTCALL Sqrt(float x);
	  float FASTCALL Fabs(float x);
	  float FASTCALL Pow(float x, float y);
	  float FASTCALL Log(float x);
	  float FASTCALL Log10(float x);
	 // float FASTCALL Exp(float x);

	inline float Deg2Rad(float val)	 {return (float(val)*DegToRad);}
	inline float Rad2Deg(float val)	 {return (float(val)*RadToDeg);}
	//inline float Random()	 {return ((float)rand()/(float)RAND_MAX );}	
}

//ﾒﾔﾏﾂｺｯﾊ?ﾎｪﾗﾔｱ滷ﾑ?ﾌ?ｺｯﾊ?｣ｬﾐｧﾂﾊﾆ貭?｣ｬﾉ?ﾓﾃ｡｣
 //ﾐｴｷｨﾊﾇｶﾔｵﾄ｣ｬｵｫﾊﾇﾐｧﾂﾊﾌｫｵﾍ｡｣
  float		BernsteinFun(int n, int i,float t); 
  double		BSplineFunRecursive(int n,int i,double t,const varray<double>& knots);  //ｴﾋｺｯﾊ?ﾗﾔｱ爛ｬﾐｧﾂﾊﾆ貭?
  float		BSplineFun(int n,int i,float t);    //
  void		BSplineFunArray(int n,float t,varray<float>& allRatio);
  double      TwoValueDivide(double vup,double vdown);

//ﾒﾔﾏﾂﾎｪｱ?ﾗｼｺｯﾊ?ｿ?
//ﾒﾔﾏﾂｼｸｸ?ｺｯﾊ?ﾊﾇBﾑ?ﾌ?ﾏ犹ﾘｺｯﾊ?｡｣
//order ｽﾗﾊ?｣ｬdegree ﾎｪｴﾎﾊ?｣ｬpｱ?ﾊｾdegree｣ｬ
//n+1ﾎｪｿﾘﾖﾆｵ飜?ﾄｿ｣ｬpﾎｪﾇ?ﾏﾟｽﾗﾊ?｣ｬｽﾚｵ飜?ﾄｿﾎｪm+1,U[]ﾎｪｽﾚｵ飜?ﾗ?
//ﾒ?ｴﾋU[]ﾊ?ﾗ魘?ﾐ｡ﾎｪm + 1｣ｬU[m]ｼｴｱ?ﾊｾﾗ?ｺ?ﾒｻｸ?ｽﾚｵ罍｣
//ﾒ?ｴﾋ m = n + p + 1
//NﾎｪN(i,p)ﾔﾚuﾖｵｴｦｲｻﾎｪﾁ羞ﾄｻ?ｺｯﾊ?ﾊ?ﾗ鬘｣
  int FindSpan(int n, int p, double u, const varray<double>& U);
  void BasisFuns( int i, double u, int p, const varray<double>& U, varray<double>& N);
  void dersBasisFuns(int i, double u, int p, int order, double knot[], double **ders);
  double OneBasisFun(int p, int m, const varray<double>& U, int i, double u);
  void dersOneBasisFuns(int p, int m, double U[], int i, double u, int n, double* ders);
  double** init2DArray(int x, int y);
  void free2Darray(double **array, int x);


}

#endif

