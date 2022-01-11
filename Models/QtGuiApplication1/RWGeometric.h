#pragma once
#include "CNurbs.h"
#include <string>
#include"spline.h"
#include"SplineSurface.h"
#include"SplineVolume.h"
#include"XVec.h"

#include "Nurbs.h"

using namespace std;
class RWGeometric
{
public:
	//��ȡ3d��ɢ��
	int ReadPoint(const string& path, varray<varray<point3d>>& allPts);

	//��ȡ4d��ɢ��
	int ReadPoint(const string& path, varray<varray<point4d>>& allPts);

	//��ȡBezier����
	//������������
	int ReadBezierLine(const string& path, varray<BezierLine>& lines);

	//��ȡNURBS����
	//������������
	int ReadNurbsLine(const string& path, varray<NurbsLine>& lines);

	//��ȡNURBS����
	//������������
	int ReadNurbsSurface(const string& path, varray<NurbsSurface>& surfaces);

	//��ȡNURBS��
	//������ģ������
	int ReadNurbsVol(const string& path, varray<NurbsVol>& vols);

	int ReadNurbsVol(const string& path, varray<YN::NurbsVol>& vols);

	//��ȡ3d��ɢ��
	int WritePoint(const string& path, const varray<varray<point3d>>& allPts);

	//��ȡ4d��ɢ��
	int WritePoint(const string& path, const varray<varray<point4d>>& allPts);

	//д��NURBS����
	//����д������������
	int WriteBezierLine(const string& path, const varray<BezierLine>& lines);

	//д��NURBS����
	//����д������������
	int WriteNurbsLine(const string& path, const varray<NurbsLine>& lines);

	//д��NURBS����
	//����д������������
	int WriteNurbsSurface(const string& path, const varray<NurbsSurface>& surfaces);

	//д��NURBS��
	//����д������ģ������
	int WriteNurbsVol(const string& path, const varray<NurbsVol>& vols);

	//-------------------------------------------------------------//
		//��ȡNURBS����
	//������������
	int ReadSpline(const string& path, varray<Spline>& lines);

	//��ȡNURBS����
	//������������
	int ReadSplineSurface(const string& path, varray<SplineSurface>& surfaces);

	//��ȡNURBS��
	//������ģ������
	int ReadSplineVolume(const string& path, varray<SplineVolume>& vols);

	//д��NURBS����
	//����д������������
	int WriteSpline(const string& path, const varray<Spline>& lines);

	//д��NURBS����
	//����д������������
	int WriteSplineSurface(const string& path, const varray<SplineSurface>& surfaces);

	//д��NURBS��
	//����д������ģ������
	int WriteSplineVolume(const string& path, const varray<SplineVolume>& vols);


};