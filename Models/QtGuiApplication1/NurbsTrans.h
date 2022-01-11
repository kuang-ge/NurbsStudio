#pragma once
#include<iostream>
#include <vector>
#include"sisl.h"
#include "CNurbs.h"
#include "SplineVolume.h"

//��ά������ת��ΪSISL��ʽ��
double* Point3dToSislPoint(const Vec3 & inputpoint, int length = 3);

//NURBS����ת��ΪSISL��ʽ����
SISLCurve* NurbsLineToSislLine(const Spline &inputLine);
SISLCurve* NurbsLineToSislLine(const Spline &inputLine, int dim);

//sisl����תcnurbs��ʽ
Spline SislLineToNurbsLine(SISLCurve* cruve);

//���߲���u�������꼰������
//��������ά��(�������ʧ������Ϊ�շ���-1)
//u:���߲���ֵ
//derive:��������������
//der:�󵼴���
int CalLinePoint(SISLCurve *&curve, double u, std::vector<std::vector<double>> &mderive, int der = 1);

//һ���һNurbs������,���������ϵ�uֵ
double PointIntersectNurbsLine(double* &pt, SISLCurve* &pc, int length = 3);


//���ж���NURBS�����Ƿ���ڽ���
bool ISTwoNurbsLineIntersect(SISLCurve *&curve1, SISLCurve *&curve2);

bool ISTwoNurbsLineIntersect(const Spline &line1, const Spline &line2, int dim = 2);

//��Nurbs�����󽻣����Է��������ϵ�uֵ
int TwoNurbsLineIntersect(SISLCurve *&curve1, SISLCurve *&curve2, double *&intpar1, double *&intpar2);

//��Nurbs������,���������
//���Է��������ϵ�uֵ
//�����ཻ��������غ�������Ŀ
int TwoNurbsLineIntersectVer2(SISLCurve *&curve1, SISLCurve *&curve2, double *&intpar1, double *&intpar2, int &numintpt, int &numintcu);

//ֱ����Spline��
int StraLineIntersectNurbsLine(SISLCurve *&curve, const varray<Vec3> &strLine, double *&intpar, int &numintpt, int &numintcu);

//�ж�һ�����Ƿ����CNurbs�����ߵ�Ҫ��(SISL����Nurbs����в���),�˺�����ʱ����ʹ��
bool isInterPoint(double &u,SISLCurve * curve, const Spline* nurbsCurve);

//B�������
SISLCurve* FitBspline(const varray<Vec3>&p, int degree, int dim = 3);
Spline FitBsplineCnurbs(const varray<Vec3>&p, int degree, int dim = 3);

class NurbsTrans
{
public:
	//Spline���ͺ�Cnurbs���͸�ʽת��
	static Spline CnurbslineToSpline(const NurbsLine& nl);
	static NurbsLine SplineToCnurbsline(const Spline& sl);
	static varray<Spline> ClinesToSplines(const varray<NurbsLine>& nls);
	static varray<NurbsLine> SplinesToClines(const varray<Spline>& sls);
	static SplineSurface CnurbssurfToSplinesurf(const NurbsSurface& nsf);
	static NurbsSurface SplinesurfToCnurbssurf(const SplineSurface& sf);
	static varray<SplineSurface> CsurfsToSplinesurfs(const varray<NurbsSurface>& nsfs);
	static varray<NurbsSurface> SplinesurfsToCsurfs(const varray<SplineSurface>& sfs);
	static SplineVolume CnurbsvolToSplinevol(const NurbsVol& nvol);
	static NurbsVol SplinevolToCnurbsvol(const SplineVolume& vol);
	static varray<SplineVolume> CvolsToSplinevols(const varray<NurbsVol>& nvols);
	static varray<NurbsVol> SplinevolsToCvols(const varray<SplineVolume>& vols);

	//�������߲���
	//���߽�ά
	//mode:0��ʾ���ԭ������1��ʾ�����
	static void DimReduceNurbsLines(varray<Spline>& nl, varray<varray<Spline>>& allL, int mode = 0);

	//���潵ά
	//mode:0��ʾ���ԭ������1��ʾ�����
	static void DimReduceNurbsSurfs(varray<SplineSurface>& nsf, varray<varray<SplineSurface>>& allsf, int mode = 0);

};