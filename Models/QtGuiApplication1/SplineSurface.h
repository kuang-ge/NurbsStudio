#pragma once

//#include "XFunction.h"
#include "spline.h"
#include "XBaseMesh.h"
#include "stdafx.h"
using namespace base;

namespace base {

class SplineSurface:public SplineBase
{
public:	
	enum SPLINESURFACE_MODE {
		SURFACEMODE_BSPLINE = 0,
		SURFACEMODE_BZIER,
		SURFACEMODE_NONUNI_BSPLINE,
	};
public:
	int					  m_uDegree;
	int					  m_vDegree;
	varray<Spline>        m_4sideSpline;  //four sided spline indexed u0,v0,u1,v1; �����߽�����β��ӡ�

public:
    varray<XBaseMesh>     m_vShowMeshes;
	varray<varray<Vec4>>  m_UVCtrlPts; //set all the control points.
	varray<Vec4>		  m_CtrlPts;
    varray<double>        m_uKnots;
	varray<double>        m_vKnots;
	int					  m_uNum;
	int					  m_vNum;
private:
    DWORD                 m_iSurfaceMode;

public:
	SplineSurface(void);
	SplineSurface(const SplineSurface& surface);
	SplineSurface& operator = (const SplineSurface& src);
	enum{ CLSID = 0xd4fe5339 };
	virtual DWORD ClassID() const  {return (DWORD)CLSID;}
	SplineSurface* SplineSurface::CreateMe() const;
public:
	~SplineSurface(void);
public:
	void AddCtrlPts(varray<varray<Vec4> >& CtrlPts);  //���ʰ���ij�������i��ʾV���кţ�j��ʾU��j�С�
	void AddKnotsVector(varray<double>& knots,bool uvflag);
	void AddKnotsVector(varray<float>& knots,bool uvflag);
	void ChangeSurfaceMode(DWORD iMode);
	void SetSurfaceDegree(int uDegree,int vDegree);
	void Generate4SidedBoudaryCurve();
	Vec4 GetUVPoint(float u,float v);
	void CreateShowMesh(bool bClear= false);
	void CopySurfaceWithoutMesh(SplineSurface& sp);
	void SetSurface(const int uDegree, const int vDegree, const int uNum, const int vNum,
		const varray<double>& uKnots, const varray<double>& vKnots);

	//���㣨u��v����Ӧ�������ϵĵ�
	Vec3 GetSurFacePoint(const double u, const double v)const;

	//�����ı�����Ƭ��ʾ����
	//num:�÷������������
	//quads:��Ƭ����
	//lines:�Ȳ�������
	threadParamSpline CalQuads(const int Unum, const int Vnum, varray<varray<Vec3>>& quads, varray<varray<Vec3>>& lines)const;

	/*Coons��ֵ
	EndgCtrlPts:�߽���Ƶ㣨e0,e1,e2,e3˳��*/
	void CoonsInterpolate(const varray<varray<Vec4>>& EdgeCtrlPts);

	/*Coons��ֵ
	EdgeLines:�߽����ߣ�e0,e1,e2,e3˳��*/
	void CoonsInterpolate(const varray<Spline>& EdgeLines);

	//��������
	//Udegree,Vdegree:���׺����
	void DegreeElevate(const int Udegree, const int Vdegree);

	//����ڵ����
	//Uknot,Vknot:��Ҫ����Ľڵ�
	void KnotsRefine(const varray<double>& Uknot, const varray<double>& Vknot);

	//����ֶ�ΪBezier����
	//dir:0=U����1=V����
	//QW��Bezier������Ƶ�
	void Decompose(const bool dir, varray<varray<varray<Vec4>>>& QW);

	//���ݿ��Ƶ��ά��ż���һά���
	int CtrlPtsIdx(const int uIdx, const int vIdx);

	//���Ƶ�����ת��ΪU-V
	void OrderCtrlPts();

	//���Ƶ�����ת��ΪU-V
	void OrderCtrlPts(SplineSurface& sf);

	//��ȡ�߽���,�����Ѿ���U-V�����
	void GetEdgeLines(varray<Spline>& EdgeLines);

//display
    /*void DisplaySplineSurface(int renderMode,bool bShowBoundary = false);
	void DrawSurfaceCtrlPts(COLORREF clr);*/
	void Clear();


public:
	Vec4 GetCtrlPt(int vi,int uj) const {return m_UVCtrlPts.at(vi).at(uj);}
	void SetCtrlPt(int vi,int uj,Vec4 vt)  {m_UVCtrlPts.at(vi).at(uj) = vt;}
	int  GetUCtrlPtNum(){return m_uNum;}
	int  GetVCtrlPtNum(){return m_vNum;}
    int  GetUDegree(){return m_uDegree;}
	int  GetVDegree(){return m_vDegree;}
//////////////////////////////////////////////////////////////////////
public:   //���º�����Ҫ����polyIGA�е���������
	bool IstwoSurfaceSame(SplineSurface& sf, bool notConsiderUVDirection = true);
	int  RotateTheWholeSurfaceAccordingtoBasesurface(int istate); //ͬʱ����һ��ֵ����ֵ��ӳ������ͻ�׼�����λ�ã���0,1,2,3λ�����С�
	void RotateBaseSurface();  
	int  RotateLastSuface(int istate);  //��ת�ⶥ���档Ҫ��istate�Ǹ����u0������λ�ã���ʱu0���Ѿ���ת��λ�������Ҫ���⴦��
	bool IsTwoVec4arraySameinOrder(const varray<Vec4>& arr1, const varray<Vec4>& arr2);
	int  GetTwoSurfaceEdgeshareState(const SplineSurface& sf);
	int  GetTwoSurfaceEdgeshareStateNew(const SplineSurface& sf);
	void RotateControlPoints(int imode);
	//���Ϻ�����Ҫ����polyIGA�е���������
/////////////////////////////////////////////////////////////////////
//B�������������
public:
	void CoonsInterpolate(const varray<varray<Vec4>>& EdgeCtrlPts, varray<varray<Vec4>>& coonsPatchCtrlpts);
	float GetDistanceBetweenTwoOrderSurface(SplineSurface& sf);
	float GetDistanceBetweenTwoSurface(SplineSurface& sf);
};
}