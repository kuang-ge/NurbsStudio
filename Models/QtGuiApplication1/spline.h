// Spline.h: Spline クラスのインターフェイス
//
//////////////////////////////////////////////////////////////////////

#ifndef _SPLINE_H__INCLUDED_
#define _SPLINE_H__INCLUDED_

#include <Windows.h>
#include "varray.h"
#include "XVec.h"
#include "MathUtil.h"

using namespace std;
using namespace base;
namespace base {
//////////////////////////////////
	class SplineBase
	{
	public:

		/*ｼﾆﾋ羶ﾚｵ耘ﾂｱ?
		x｣ｺｽﾚｵ聊ｵ
		degree｣ｺｴﾎﾊ?
		CtrlPtsNum｣ｺｿﾘﾖﾆｵ飜?ﾁｿ
		knots｣ｺｽﾚｵ飜ｸﾁｿ*/
		int FindSpan(const double x, const int degree, const int CtrlPtsNum, const varray<double>& knots)const;

		/*ｸ?ｾﾝｲﾎﾊ?ﾖｵ｣ｬｽﾚｵ耘ﾂｱ凜ｬｼﾆﾋ羹?ｺｯﾊ?
		u｣ｺｲﾎﾊ?ﾖｵ
		k｣ｺｽﾚｵ耘ﾂｱ?
		degree｣ｺｴﾎﾊ?
		knots｣ｺｽﾚｵ飜ｸﾁｿ
		N｣ｺｷｵｻﾘｵﾄｻ?ｺｯﾊ?*/
		void BasisFuns(const double u, const int k, const int degree, const varray<double>& knots,
			varray<double>& N)const;

		/*uｴｦﾋ?ﾓﾐｻ?ｺｯﾊ?
		u｣ｺｲﾎﾊ?ﾖｵ
		k｣ｺｽﾚｵ耘ﾂｱ?
		degree｣ｺｴﾎﾊ?
		knots｣ｺｽﾚｵ飜ｸﾁｿ
		ndu｣ｺｷｵｻﾘｵﾄﾋ?ﾓﾐｻ?ｺｯﾊ?*/
		void AllBasisFuns(const double u, const int k, const int degree,
			const varray<double>& knots, varray<varray<double>>& ndu)const;

		/*ｻ?ｺｯﾊ?nｽﾗｵｼ
		u｣ｺｲﾎﾊ?ﾖｵ
		k｣ｺｽﾚｵ耘ﾂｱ?
		degree｣ｺｴﾎﾊ?
		n｣ｺｵｼﾊｸｽﾗﾊ?
		knots｣ｺｽﾚｵ飜ｸﾁｿ
		basisDu｣ｺｻ?ｺｯﾊ?nｽﾗｵｼ*/
		void DerBasisFuns(const double u, const int k, const int degree, const int n, const varray<double>& knots,
			varray<varray<double>>& basisDu)const;
	};
class  Spline: public SplineBase
{
public:
	varray<Vec4> m_CtrlPts;
	DWORD m_mode;
	int m_Degree;
	varray<double> m_Knots;
public:
	
	enum SPLINE_MODE {
		//the segment between two control points in a spline with the following mode is decided by the two control points and its tangent.
		SPLMODE_LINER = 0,
        //t = [0,1] for the whole spline. 
        SPLMODE_BEZIER,          
        SPLMODE_CLOSED_BEZIER,	//not realized.
        //Spline. t = [0,1] for the every segment between two control points of the whole spline.
		SPLMODE_SPLINE,         
        SPLMODE_CLOSED_SPLINE,		
		//uniform Bspline.
		SPLMODE_2BSPLINE,
		SPLMODE_CLOSED_2BSPLINE,
        SPLMODE_3BSPLINE,
        SPLMODE_CLOSED_3BSPLINE,
		//rational Bspline
        SPLMODE_NBSPLINE,  //n = 2,3,4,......  Nｴﾎｾ?ﾔﾈﾓﾐﾀ?Bﾑ?ﾌ?   
		SPLMODE_CLOSED_NBSPLINE,  //ｷ箜ﾕｾ?ﾔﾈﾓﾐﾀ?Bﾑ?ﾌ?
		//non-uniform Bspline
        SPLMODE_NONUNI_NBSPLINE,
	};
public:
	bool				m_bSplineIsEdit;
	int					m_nSelected;
	int					m_nMaped;
	float				m_yCrossSection;
protected:
	int				m_bSplineGenerateMode;	//generate mode: 0:normal mode; 1:giving tangent of the control point directly for spline;2 for boundary spline					
	varray<Vec4>	m_vTangent;				//arry for tangent
	varray<float>   m_vTangentMagnitude0;   //the tangent of each control point i have two magnitude. magnitude0 for (i-1, i) and magnitude1 for(i, i+1)
	varray<float>   m_vTangentMagnitude1;   //		
	varray<Vec4>    m_RenderPt;    //ｵ耨ﾃﾓﾚﾏﾔﾊｾ｡｣
	
		/*ｼﾆﾋ縉ｽﾚｵ羝ﾔﾓｦｵﾄﾇ?ﾏﾟﾗ?ｱ?ｵ?
		u｣ｺｽﾚｵ羇ﾎﾊ?*/
public:
	Vec3 GetLinePoint(const double u)const;

	/*ｻ?ﾈ｡ﾏﾟﾉﾏｵﾄｵ?
	Unum｣ｺﾇ?ﾏﾟｵ飜?ﾁｿ
	linePts｣ｺﾇ?ﾏﾟｵ?	*/
	void CalLinePoint(const int Unum, varray<Vec3>& linePts)const;

	/*ﾇ?ﾏﾟﾊｸﾖｵｺｯﾊ?A(u)ﾋ?ﾓﾐpｽﾗｵｼ｣ｬﾊｸﾖｵｺｯﾊ?A(u)ｼｴNURBSﾇ?ﾏﾟｷﾖﾗﾓﾊｽ
	u｣ｺｲﾎﾊ?
	n｣ｺｵｼﾊｸｽﾗﾊ?
	Der｣ｺA(u)ﾋ?ﾓﾐpｽﾗｵｼ*/
	void PtsDerivsAu(const double u, const int n, varray<Vec4>& Der)const;

	/*ﾇ?ﾏﾟﾉﾏuｴｦﾋ?ﾓﾐnｽﾗｵｼ
	u｣ｺｲﾎﾊ?
	n｣ｺｵｼﾊｸｽﾗﾊ?
	Der｣ｺﾋ?ﾓﾐnｽﾗｵｼ */
	void PtsDerivs(const double u, const int n, varray<Vec4>& Der)const;

	/*ﾇ?ﾏﾟﾋ?ﾓﾐnｽﾗｵｼ
	step｣ｺｲﾎﾊ?ﾓ?ｵﾈｷﾖｲｽｽ?ﾁｿ
	Der｣ｺﾇ?ﾏﾟｵ羝ﾔﾓｦｵﾄnｽﾗｵｼ*/
	void CurveDerivs(const int n, const double step, varray<varray<Vec3>>& Der);

	/*ｽﾚｵ羇衒?
	u｣ｺﾐ靨ｪｲ衒?ｵﾄｽﾚｵ?
	r｣ｺｲ衒?ｴﾎﾊ?*/
	void KnotInsert(const double u, const int r);

	/*ｽﾚｵ翳･ｳ?
	u｣ｺﾐ靨ｪﾉｾｳ?ｵﾄｽﾚｵ?
	r｣ｺﾉｾｳ?ｴﾎﾊ?*/
	int KnotRemove(const double u, const int r);

	/*ｽﾚｵ耘ｸｻｯ
	u:ﾐ靨ｪｲ衒?ｵﾄｽﾚｵ?*/
	void KnotsRefine(const varray<double>& u);

	/* ﾇ?ﾏﾟﾉ?ｽﾗ
	degree｣ｺﾉ?ｽﾗｺ?ｴﾎﾊ?*/
	void DegreeElevate(const int degree);

	//ﾔﾚｽﾚｵ縉ｴｦｶﾔﾇ?ﾏﾟｽ?ﾐﾐｷﾖｸ?
	void Segmentation(const double u, varray<Spline>& lines);

	/*NURBSﾇ?ﾏﾟｷﾖｽ簧ｪBEZIER
	Qw｣ｺﾊ莎?ｵﾄBEZIERﾇ?ﾏﾟｶﾎｿﾘﾖﾆｵ?*/
	void Decompose(varray<varray<Vec4>>& Qw);



public:
	//ﾇ?ﾏﾟﾉ?ｳﾉ
	Spline();
	Spline(const Spline& spl);
	virtual ~Spline();
	enum{ CLSID = 0xd4fe5009 };

	virtual DWORD ClassID() const  {return (DWORD)CLSID;}
	Spline* Spline::CreateMe() const;//ｰﾑﾕ篋?ﾇ?ﾏﾟｸｴﾖﾆﾒｻｸ?ｳ?ﾀｴ
	void ClearStatus();//ﾇ蠡?ｺｯﾊ?
    
	//ﾇ?ﾏﾟｵﾄｴ豢｢ｺﾍｶﾁﾈ｡
	bool Save(BaseOStream& bs) const;
	bool Load(BaseIStream& bs);	
	bool Save0010(BaseOStream& bs) const;
	bool Load0010(BaseIStream& bs,float fCurVersion);  
	//bool Load0010(FILE* fp,float fCurVersion);  //just  version 1.005   ﾃｻﾓﾐﾊｵﾏﾖ
	
    //Bﾑ?ﾌ?ﾇ?ﾏﾟﾄ篌ﾏ
	void    FitSpline(varray<Vec4>& inputVecs,varray<float>& knots, int degreeNum,DWORD splineMode = SPLMODE_NBSPLINE);
	//ﾉ靹ﾃｺﾍｻ?ﾈ｡ｲﾙﾗ?
	void	SetTangentOfCtrlPnt(Vec4& vTangent, int i);
	void	SetSplineGenerateMd(int iMode);	
	void	SetMagnitudeOfTangent(int i, int j, float fMagnitude); //get the magnitude of control point, 20050125	
    void    SetCtrlPoint(int index,const Vec4& v);
	void    SetCtrlPointArr(const varray<Vec4>& ptsArr)    { m_CtrlPts = ptsArr;}
	void    SetSplineDegree(int n = 3);
	void    SetSplineKnot(const varray<double>& knots){m_Knots = knots;}

    DWORD   GetMode() const {return m_mode;}   
	int     GetSplineDegree() const {return m_Degree;}
	Vec4&   GetCtrlPoint(int i) {return m_CtrlPts[i];}
	const Vec4& GetCtrlPoint(int i) const {return m_CtrlPts[i];}
	int     GetCtrlPointCount() const {return static_cast<int>(m_CtrlPts.size());}
	Vec4    GetPoint(int i, double t, float tension=0.0);
	float   GetLength(int n);
	Vec4			GetTangentOfCtrlPnt(int i);
	varray<Vec4>&   GetCtrlPointArr()     {return m_CtrlPts;}
	int				GetSplineGenerateMd(void)const;
	int				GetTangentCount() const { return static_cast<int>(m_vTangent.size()); }
	float			GetMagnitudeOfTangent(int i, int j);
	void CruveReverse();//ﾇ?ﾏﾟｷｴﾗｪ
	void CreatLineWithTwoPoints(const Vec3 &p1, const Vec3 &p2, int degree = 2);

    //ｿﾘﾖﾆｵ羞ﾄﾒｻﾐｩﾆ萢?ｲﾙﾗ?
	void    AddTangent(Vec4 vTangent, float fMagnitude0, float fMagnitude1);
	void    ChangeMode(DWORD mode,int degree = 0);		
	void    AddCtrlPoint(const Vec4& v);
	void    AddKnots(const double knot);
	void    DelCtrlPoint(int index);    
	void    InsertCtrlPoint(int idx,Vec4 vt);	
	void    ClearCtrlPoint();  	

	void    resize(int sz) {m_CtrlPts.resize(sz);}
	int     size() const {return static_cast<int>(m_CtrlPts.size());}
	Vec4&   p(int i) {return m_CtrlPts[i];}//ｷｵｻﾘｿﾘﾖﾆｵ羞ﾄｵﾚiｸ?ｵ?
	const Vec4& p(int i) const {return m_CtrlPts[i];}
	Vec4&   pBack() {return m_CtrlPts.back();}
	const Vec4&  pBack()const {return m_CtrlPts.back();} 
};
}
#endif
