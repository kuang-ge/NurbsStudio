// Spline.h: Spline �N���X�̃C���^�[�t�F�C�X
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

		/*����ڵ��±�
		x���ڵ�ֵ
		degree������
		CtrlPtsNum�����Ƶ�����
		knots���ڵ�ʸ��*/
		int FindSpan(const double x, const int degree, const int CtrlPtsNum, const varray<double>& knots)const;

		/*���ݲ���ֵ���ڵ��±꣬���������
		u������ֵ
		k���ڵ��±�
		degree������
		knots���ڵ�ʸ��
		N�����صĻ�����*/
		void BasisFuns(const double u, const int k, const int degree, const varray<double>& knots,
			varray<double>& N)const;

		/*u�����л�����
		u������ֵ
		k���ڵ��±�
		degree������
		knots���ڵ�ʸ��
		ndu�����ص����л�����*/
		void AllBasisFuns(const double u, const int k, const int degree,
			const varray<double>& knots, varray<varray<double>>& ndu)const;

		/*������n�׵�
		u������ֵ
		k���ڵ��±�
		degree������
		n����ʸ����
		knots���ڵ�ʸ��
		basisDu��������n�׵�*/
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
        SPLMODE_NBSPLINE,  //n = 2,3,4,......  N�ξ�������B����   
		SPLMODE_CLOSED_NBSPLINE,  //��վ�������B����
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
	varray<Vec4>    m_RenderPt;    //��������ʾ��
	
		/*����u�ڵ��Ӧ�����������
		u���ڵ����*/
public:
	Vec3 GetLinePoint(const double u)const;

	/*��ȡ���ϵĵ�
	Unum�����ߵ�����
	linePts�����ߵ�	*/
	void CalLinePoint(const int Unum, varray<Vec3>& linePts)const;

	/*����ʸֵ����A(u)����p�׵���ʸֵ����A(u)��NURBS���߷���ʽ
	u������
	n����ʸ����
	Der��A(u)����p�׵�*/
	void PtsDerivsAu(const double u, const int n, varray<Vec4>& Der)const;

	/*������u������n�׵�
	u������
	n����ʸ����
	Der������n�׵� */
	void PtsDerivs(const double u, const int n, varray<Vec4>& Der)const;

	/*��������n�׵�
	step��������ȷֲ�����
	Der�����ߵ��Ӧ��n�׵�*/
	void CurveDerivs(const int n, const double step, varray<varray<Vec3>>& Der);

	/*�ڵ����
	u����Ҫ����Ľڵ�
	r���������*/
	void KnotInsert(const double u, const int r);

	/*�ڵ�ȥ��
	u����Ҫɾ���Ľڵ�
	r��ɾ������*/
	int KnotRemove(const double u, const int r);

	/*�ڵ�ϸ��
	u:��Ҫ����Ľڵ�*/
	void KnotsRefine(const varray<double>& u);

	/* ��������
	degree�����׺����*/
	void DegreeElevate(const int degree);

	//�ڽڵ�u�������߽��зָ�
	void Segmentation(const double u, varray<Spline>& lines);

	/*NURBS���߷ֽ�ΪBEZIER
	Qw�������BEZIER���߶ο��Ƶ�*/
	void Decompose(varray<varray<Vec4>>& Qw);



public:
	//��������
	Spline();
	Spline(const Spline& spl);
	virtual ~Spline();
	enum{ CLSID = 0xd4fe5009 };

	virtual DWORD ClassID() const  {return (DWORD)CLSID;}
	Spline* Spline::CreateMe() const;//��������߸���һ������
	void ClearStatus();//������
    
	//���ߵĴ洢�Ͷ�ȡ
	bool Save(BaseOStream& bs) const;
	bool Load(BaseIStream& bs);	
	bool Save0010(BaseOStream& bs) const;
	bool Load0010(BaseIStream& bs,float fCurVersion);  
	//bool Load0010(FILE* fp,float fCurVersion);  //just  version 1.005   û��ʵ��
	
    //B�����������
	void    FitSpline(varray<Vec4>& inputVecs,varray<float>& knots, int degreeNum,DWORD splineMode = SPLMODE_NBSPLINE);
	//���úͻ�ȡ����
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
	void CruveReverse();//���߷�ת
	void CreatLineWithTwoPoints(const Vec3 &p1, const Vec3 &p2, int degree = 2);

    //���Ƶ��һЩ��������
	void    AddTangent(Vec4 vTangent, float fMagnitude0, float fMagnitude1);
	void    ChangeMode(DWORD mode,int degree = 0);		
	void    AddCtrlPoint(const Vec4& v);
	void    AddKnots(const double knot);
	void    DelCtrlPoint(int index);    
	void    InsertCtrlPoint(int idx,Vec4 vt);	
	void    ClearCtrlPoint();  	

	void    resize(int sz) {m_CtrlPts.resize(sz);}
	int     size() const {return static_cast<int>(m_CtrlPts.size());}
	Vec4&   p(int i) {return m_CtrlPts[i];}//���ؿ��Ƶ�ĵ�i����
	const Vec4& p(int i) const {return m_CtrlPts[i];}
	Vec4&   pBack() {return m_CtrlPts.back();}
	const Vec4&  pBack()const {return m_CtrlPts.back();} 
};
}
#endif
