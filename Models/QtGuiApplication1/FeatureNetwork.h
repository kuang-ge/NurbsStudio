#pragma once
#include"CNurbs.h"
#include"globalFunc.h"
#include"pointXd.h"
#include"varray.h"
#include"RWGeometric.h"
#include "quadPart.h"
#include "NurbsTrans.h"
#include<vector>
#include<math.h>
#include<iostream>
#include<map>
#include<thread>
#include<fstream>
#include <algorithm>
//��������࣬�����˼��ֻ�е����ĳ�����״���Լ��������Ҷ�ֵ�NURBS��ģ�塣
//���в����������Ӳ��޸İ汾��
//�汾20201126

#define PI 3.1415926535
#define x_forward 1
#define x_backward -1
#define y_forward 2
#define y_backward -2
#define z_forward 3
#define z_backward -3
#define ERRF 0.0001
using namespace std;

//����Ԫ�ؼ���
//������
class Feature_Line :public Spline
{
public:
	Feature_Line() {}

	Feature_Line(const Spline& nl, bool seg = true);

	//ͨ�����㹹��ֱ����ʽ������
	Feature_Line(const Vec3& p1, const Vec3& p2, bool seg = true);

	//�������
	Spline GetLine();

public:
	//varray<int> m_subLinel;		//���߷ָ���Ӧ�Ķ���������
	bool m_seg;					//���߿ɷָ���
};

//������
class Feature_Surface
{
public:
	Feature_Surface() {}

	//���캯�������ɹ�����ɺ���������ݣ�����������
	Feature_Surface(const varray< Feature_Line>& fls, bool empty = false, int sfnum = 0);
	Feature_Surface(const varray< Feature_Line>& fls, const varray<Vec3>& ps, bool empty = false, int sfnum = 0);
	Feature_Surface(const SplineSurface& s);
	Feature_Surface(const varray<SplineSurface>& s);

	//�������
	varray<SplineSurface> GetSurfaces()const;


public:
	varray<Feature_Line> m_edges;		//�߽�������
	varray<Vec3> m_edgePoints;		//�߽綥��
	//varray<int> m_edgeNum;				//��Ӧ���߼����е����
	int m_sfNum;						//��������
	varray<SplineSurface> m_ns;			//������ɵ�NURBS����
	bool m_iscomplete;					//�Ƿ񴴽����
	bool m_isempty;						//�������Ƿ�Ϊ��(����)��Ĭ��Ϊfalse�ǿ�
};

class Feature_Vol
{
public:
	Feature_Vol() {}

	Feature_Vol(const Feature_Surface& fs, const Feature_Line& path, int mode, int dir = 0, double dis = 0);
	Feature_Vol(const varray<Feature_Surface>& fs, const Feature_Line& path, int mode, int dir = 0, double dis = 0);
	Feature_Vol(const Feature_Surface& fs, const varray<Feature_Line>& path, int mode, int dir = 0, double dis = 0);
	Feature_Vol(const varray<Feature_Surface>& fs, const varray<Feature_Line>& path, int mode, int dir = 0, double dis = 0);
	//���췽ʽ���캯��
	Feature_Vol(const Feature_Surface& fs, int mode = 1, int dir = 0, double dis = 0);
	Feature_Vol(const varray<Feature_Surface>& fs, int mode = 1, int dir = 0, double dis = 0);

	//����NURBS��
	//del:�����������Ƿ�ɾ����������·������
	void ConstructVold(bool del = false);

	//�Ƿ����NURBS�幹��
	bool Iscomplete()const;

	//���NURBS��
	varray<SplineVolume> GetVOls()const;

public:
	varray<Feature_Surface> m_surf;		//������
	varray<Feature_Line> m_path;		//·��(����õ�)
	//varray<int> m_pathNum;				//·�������߼����е����
	int m_mode;							//���ͷ���(����-1;ɨ��-2(���洹ֱ·��);ɨ��-3(���治��ֱ·��)����-4;��ֵ-5)
	int m_dir;							//ֵ��Ϊ0��ʾ���õ������죬(123->+xyz)
	double m_dis;						//���쳤��

private:
	varray<SplineVolume> m_fv;				//�������NURBS��
	bool m_iscomplete;					//������ɱ��
	friend class Model_Solution;
};





//����������
class Creat_Vol
{
public:
	SplineVolume NV;
	Spline NL;

	double m_PathLength;
	int mode;
	int flag = 0;
	void CreatPath()
	{
		Vec3 Ori;
		switch (mode) {
		case 1:
			Ori = { 1,0,0 };
			break;
		case -1:
			Ori = { -1,0,0 };
			break;
		case 2:
			Ori = { 0,1,0 };
			break;
		case -2:
			Ori = { 0,-1,0 };
			break;
		case 3:
			Ori = { 0,0,1 };
			break;
		case -3:
			Ori = { 0,0,-1 };
			break;
		}
		Vec3 Path = Ori * m_PathLength;
		Vec3 Path_Position0 = { 0,0,0 };
		Vec3 Path_Position1 = Path / 2;
		Vec3 Path_Position2 = Path;
		NL.m_CtrlPts.push_back(Path_Position0);
		NL.m_CtrlPts.push_back(Path_Position1);
		NL.m_CtrlPts.push_back(Path_Position2);
		NL.m_Degree = 2;
		NL.m_Knots.push_back(0);
		NL.m_Knots.push_back(0);
		NL.m_Knots.push_back(0);
		NL.m_Knots.push_back(1);
		NL.m_Knots.push_back(1);
		NL.m_Knots.push_back(1);
	}

	void CreatSweepVol(SplineSurface NS)
	{
		flag = 1;
		NV.CreateTransSweepSplineVolume(NL, NS);
	}

	void CreatLoftingVol(SplineSurface NS1, SplineSurface NS2)
	{
		flag = 2;
		varray<SplineSurface> NSS;
		NSS.push_back(NS1);
		NSS.push_back(NS2);
		NV.LoftingSplineVolume(NL, NSS);
	}

	void InitVol(SplineSurface NS, double length, int mode)
	{
		this->mode = mode;
		this->m_PathLength = length;
		CreatPath();
		CreatSweepVol(NS);
	}

	void InitVol(SplineSurface S1, SplineSurface S2, double length, int mode)
	{
		this->mode = mode;
		this->m_PathLength = length;
		CreatPath();
		CreatLoftingVol(S1, S2);
	}

};

//���ú���
class Model_Solution{

//���浥����ƽ��,mode1,2,3�ֱ�Ϊxyz����
public:
	//����Բ�����е�
	Vec4 GetCirclrPts(Vec3 p1, Vec3 p2, Vec3 p3)
	{
		double x1 = p1.x, y1 = p1.y;
		double x2 = p2.x, y2 = p2.y;
		double x3 = p3.x, y3 = p3.y;
		double a = 2 * (x2 - x1);
		double b = 2 * (y2 - y1);
		double c = x2 * x2 + y2 * y2 - x1 * x1 - y1 * y1;
		double d = 2 * (x3 - x2);
		double e = 2 * (y3 - y2);
		double f = x3 * x3 + y3 * y3 - x2 * x2 - y2 * y2;
		double x = (b*f - e * c) / (b*d - e * a);
		double y = (d*c - a * f) / (b*d - e * a);
		double r = sqrt((x - x1)*(x - x1) + (y - y1)*(y - y1));

		double k1 = -1 * (x1 - x) / (y1 - y);
		double b1 = y1 - k1 * x1;
		double k2 = -1 * (x3 - x) / (y3 - y);
		double b2 = y3 - k2 * x3;
		double Xres = (b2 - b1) / (k1 - k2);
		double Yres = k1 * Xres + b1;
		double l = sqrt((Xres - x1)*(Xres - x1) + (Yres - y1)*(Yres - y1));
		double beta0 = atan(l / r);
		return Vec4(Xres, Yres, 0, cos(beta0));
	}
	//����Բ���е㣬p1,p2Ϊ�˵㣬CptsΪԲ��
	Vec4 GetTangentPts(Vec4 p1, Vec4 p2, Vec4 Cpts)
	{
		double R = sqrt((Cpts.x - p1.x)*(Cpts.x - p1.x) + (Cpts.y - p1.y)*(Cpts.y - p1.y));
		double x1 = p1.x, y1 = p1.y;
		double x2 = p2.x, y2 = p2.y;
		double a = Cpts.x, b = Cpts.y;
		double k1 = -(x1 - a) / (y1 - b), k2 = -(x2 - a) / (y2 - b);
		double b1 = y1 - k1 * x1, b2 = y2 - k2 * x2;
		double x = -(b1 - b2) / (k1 - k2);
		double y = k1 * x + b1;
		double r0 = R;
		double l = sqrt((x - x1)*(x - x1) + (y - y1)*(y - y1));
		double delt0 = asin(l / r0);
		return Vec4(x, y, 0, cos(delt0));
	}


	static Vec4 NurbsCirle_mid_pt(Vec3& p1,Vec3& p2,double xt, double yt, double r) {
		double x0 = p1.x, y0 = p1.y;
		double x1 = p2.x, y1 = p2.y;
		double x, y;
		if (abs(y1 - yt) > 0.001 && abs(y0 - yt) > 0.001)
		{
			double k1 = -(x1 - xt) / (y1 - yt), k2 = -(x0 - xt) / (y0 - yt);
			x = (y0 - y1 + k1 * x1 - k2 * x0) / (k1 - k2);
			y = k1 * (x - x1) + y1;
		}
		else if (abs(y1 - yt) <= 0.001) {
			double k0 = -(x0 - xt) / (y0 - yt);
			x = x1;
			y = y0 + k0 * (x - x0);
		}
		else if (abs(y0 - yt) <= 0.001) {
			double k1 = -(x1 - xt) / (y1 - yt);
			x = x0;
			y = y1 + k1 * (x - x1);
		}
		double e = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1)) / 2;
		double f = sqrt((x0 - x)*(x0 - x) + (y0 - y)*(y0 - y));
		double w = e / f;
		return { x,y,0,w };
	}
	//�ж������Ƿ��غ�
	bool JudgeTwoPointsCoincide(const Vec3 & p1, const Vec3 & p2);

	//����ɱպ϶���ε��������߰���ʱ������,������ʱ�����������
	//polLines:���������α����߼���
	varray<Vec3> OrderLinesAntioclock(varray<Spline> &polLines)
	{
		int line1 = -1, line2 = -1;									//��ѡȡ�ĵ�һ���������
		Vec3 tmpPoint;
		Vec3 secondPoint1, secondPoint2;								//��ѡȡ�ĵڶ����㣬����ѡȡ��Ҫ��
		varray<Vec3> orderedPoints;
		if (polLines.size() < 4)										//���������߼�����������4������ӡerror������
		{
			cout << "OrderLinesAntioclock��LinesSize<4" << endl;
			system("pause");
			return orderedPoints;
		}

		tmpPoint = polLines[0].m_CtrlPts[0];
		for (auto i : polLines)											//�ҵ�y����ֵ��С�ĵ���Ϊ��ʼ��
		{
			if (i.m_CtrlPts[0].y < tmpPoint.y)
			{
				tmpPoint = i.m_CtrlPts[0];
			}
			if ((i.m_CtrlPts.end() - 1)->y < tmpPoint.y)
			{
				tmpPoint = *(i.m_CtrlPts.end() - 1);
			}
		}
		orderedPoints.push_back(tmpPoint);								//����ʼ���������õļ���

		for (int i = 0; i < polLines.size(); i++)						//�ҵ���ʼ������������������
		{
			if (JudgeTwoPointsCoincide(tmpPoint, polLines[i].m_CtrlPts[0]))
			{
				if (line1 == -1 && line2 == -1)
				{
					secondPoint1 = *(polLines[i].m_CtrlPts.end() - 1);
					line1 = i;
				}
				else
				{
					secondPoint2 = *(polLines[i].m_CtrlPts.end() - 1);
					line2 = i;
				}
			}
			else if (JudgeTwoPointsCoincide(tmpPoint, *(polLines[i].m_CtrlPts.end() - 1)))
			{
				if (line1 == -1 && line2 == -1)
				{
					secondPoint1 = polLines[i].m_CtrlPts[0];
					line1 = i;
				}
				else
				{
					secondPoint2 = polLines[i].m_CtrlPts[0];
					line2 = i;
				}
			}
			if (line1 != -1 && line2 != -1)
			{
				break;
			}
		}

		Vec3 vec1 = secondPoint1 - tmpPoint;							//��ʼ����������1
		Vec3 vec2 = secondPoint2 - tmpPoint;							//��ʼ����������2
		vec1 = vec1.Normalize();
		vec2 = vec2.Normalize();
		Vec3 cross = vec1.CrossVecX(vec2);
		if (cross.z > 0)												//secondPoint1Ϊ�����
		{
			tmpPoint = secondPoint1;
			orderedPoints.push_back(tmpPoint);
			//��������˳��
			Spline tmpLine = polLines[line1];
			polLines[line1] = polLines[0];
			polLines[0] = tmpLine;
		}
		else if (cross.z < 0)											//secondPoint2Ϊ�����
		{
			tmpPoint = secondPoint2;
			orderedPoints.push_back(tmpPoint);
			Spline tmpLine = polLines[line2];
			polLines[line2] = polLines[0];
			polLines[0] = tmpLine;
		}
		else if (fabs(cross.z - 0) < 1e-5)								//��������180��
		{
			if (vec1.x > vec2.x)
			{
				tmpPoint = secondPoint1;
				orderedPoints.push_back(tmpPoint);
				Spline tmpLine = polLines[line1];
				polLines[line1] = polLines[0];
				polLines[0] = tmpLine;
			}
			else
			{
				tmpPoint = secondPoint2;
				orderedPoints.push_back(tmpPoint);
				Spline tmpLine = polLines[line2];
				polLines[line2] = polLines[0];
				polLines[0] = tmpLine;
			}
		}

		//ѭ�������ʱ��˳���
		while (orderedPoints.size() < polLines.size())
		{
			for (int i = orderedPoints.size() - 1; i < polLines.size(); i++)
			{
				if (JudgeTwoPointsCoincide(tmpPoint, polLines[i].m_CtrlPts[0]))
				{
					tmpPoint = *(polLines[i].m_CtrlPts.end() - 1);
					//��������
					Spline tmpLine = polLines[i];
					polLines[i] = polLines[orderedPoints.size() - 1];
					polLines[orderedPoints.size() - 1] = tmpLine;
					//�����µĵ�
					orderedPoints.push_back(tmpPoint);
					break;
				}
				else if (JudgeTwoPointsCoincide(tmpPoint, *(polLines[i].m_CtrlPts.end() - 1)))
				{
					tmpPoint = polLines[i].m_CtrlPts[0];
					//��������
					Spline tmpLine = polLines[i];
					polLines[i] = polLines[orderedPoints.size() - 1];
					polLines[orderedPoints.size() - 1] = tmpLine;
					//�����µĵ�
					orderedPoints.push_back(tmpPoint);
					break;
				}
				if (i == polLines.size() - 1)
				{
					cout << "OrderLinesAntioclock��δ�ҵ���������" << endl;	//������������δ�ҵ���������һ������
					system("pause");
				}
			}
		}

		//�������һ�����Ƿ���ȷ
		if (JudgeTwoPointsCoincide(orderedPoints[0], (polLines.end() - 1)->m_CtrlPts[0])
			|| JudgeTwoPointsCoincide(orderedPoints[0], *((polLines.end() - 1)->m_CtrlPts.end() - 1)))
		{
			cout << "���һ������Ϊ��Ӧ����" << endl;
		}

		return orderedPoints;
	}

	//���ߵ�����ƽ��
	void Trans(Spline &SL, double distance, int mode)
	{

		if (mode == 1) {
			for (int i = 0; i < SL.m_CtrlPts.size(); i++)
			{
				SL.m_CtrlPts[i].x += distance;
			}
		}
		else if (mode == 2) {
			for (int i = 0; i < SL.m_CtrlPts.size(); i++)
			{
				SL.m_CtrlPts[i].y += distance;
			}
		}
		else if (mode == 3) {
			for (int i = 0; i < SL.m_CtrlPts.size(); i++)
			{
				SL.m_CtrlPts[i].z += distance;
			}
		}
		else if (mode == -1) {
			for (int i = 0; i < SL.m_CtrlPts.size(); i++)
			{
				SL.m_CtrlPts[i].x -= distance;
			}
		}
		else if (mode == -2) {
			for (int i = 0; i < SL.m_CtrlPts.size(); i++)
			{
				SL.m_CtrlPts[i].y -= distance;
			}
		}
		else if (mode == -3) {
			for (int i = 0; i < SL.m_CtrlPts.size(); i++)
			{
				SL.m_CtrlPts[i].z -= distance;
			}
		}
	}

	//���ߵ�����ƽ��
	void Trans(varray<Spline> &SL, double distance, int mode)
	{
		for (int j = 0; j < SL.size(); j++) {
			Trans(SL[j], distance, mode);
		}
	}

	//���浥����ƽ��
	void Trans(SplineSurface &NS, double distance, int mode)
	{

		if (mode == 1) {
			for (int i = 0; i < NS.m_CtrlPts.size(); i++)
			{
				NS.m_CtrlPts[i].x += distance;
			}
		}
		else if (mode == 2) {
			for (int i = 0; i < NS.m_CtrlPts.size(); i++)
			{
				NS.m_CtrlPts[i].y += distance;
			}
		}
		else if (mode == 3) {
			for (int i = 0; i < NS.m_CtrlPts.size(); i++)
			{
				NS.m_CtrlPts[i].z += distance;
			}
		}
		else if (mode == -1) {
			for (int i = 0; i < NS.m_CtrlPts.size(); i++)
			{
				NS.m_CtrlPts[i].x -= distance;
			}
		}
		else if (mode == -2) {
			for (int i = 0; i < NS.m_CtrlPts.size(); i++)
			{
				NS.m_CtrlPts[i].y -= distance;
			}
		}
		else if (mode == -3) {
			for (int i = 0; i < NS.m_CtrlPts.size(); i++)
			{
				NS.m_CtrlPts[i].z -= distance;
			}
		}
	}
	//���浥����ƽ��
	void Trans(varray<SplineSurface> &NS, double distance, int mode)
	{
		for (int j = 0; j < NS.size(); j++) {
			if (mode == 1) {
				for (int i = 0; i < NS[j].m_CtrlPts.size(); i++)
				{
					NS[j].m_CtrlPts[i].x += distance;
				}
			}
			else if (mode == 2) {
				for (int i = 0; i < NS[j].m_CtrlPts.size(); i++)
				{
					NS[j].m_CtrlPts[i].y += distance;
				}
			}
			else if (mode == 3) {
				for (int i = 0; i < NS[j].m_CtrlPts.size(); i++)
				{
					NS[j].m_CtrlPts[i].z += distance;
				}
			}
			else if (mode == -1) {
				for (int i = 0; i < NS[j].m_CtrlPts.size(); i++)
				{
					NS[j].m_CtrlPts[i].x -= distance;
				}
			}
			else if (mode == -2) {
				for (int i = 0; i < NS[j].m_CtrlPts.size(); i++)
				{
					NS[j].m_CtrlPts[i].y -= distance;
				}
			}
			else if (mode == -3) {
				for (int i = 0; i < NS[j].m_CtrlPts.size(); i++)
				{
					NS[j].m_CtrlPts[i].z -= distance;
				}
			}
		}
	}

	//���嵥����ƽ�ƣ���ȫ��һ�£���������������ƽ��
	void Trans(SplineVolume& NV, double distance, int mode)
	{

		if (mode == 1) {
			for (int i = 0; i < NV.m_CtrlPts.size(); i++)
			{
				NV.m_CtrlPts[i].x += distance;
			}
		}
		else if (mode == 2) {
			for (int i = 0; i < NV.m_CtrlPts.size(); i++)
			{
				NV.m_CtrlPts[i].y += distance;
			}
		}
		else if (mode == 3) {
			for (int i = 0; i < NV.m_CtrlPts.size(); i++)
			{
				NV.m_CtrlPts[i].z += distance;
			}
		}
		else if (mode == -1) {
			for (int i = 0; i < NV.m_CtrlPts.size(); i++)
			{
				NV.m_CtrlPts[i].x -= distance;
			}
		}
		else if (mode == -2) {
			for (int i = 0; i < NV.m_CtrlPts.size(); i++)
			{
				NV.m_CtrlPts[i].y -= distance;
			}
		}
		else if (mode == -3) {
			for (int i = 0; i < NV.m_CtrlPts.size(); i++)
			{
				NV.m_CtrlPts[i].z -= distance;
			}
		}
	}

	//���嵥����ƽ�� ��ȫ��һ�£���������������ƽ��
	void Trans(varray<SplineVolume>& NVS, double distance, int mode)
	{
		for (int j = 0; j < NVS.size(); j++) {
			if (mode == 1) {
				for (int i = 0; i < NVS[j].m_CtrlPts.size(); i++)
				{
					NVS[j].m_CtrlPts[i].x += distance;
				}
			}
			else if (mode == 2) {
				for (int i = 0; i < NVS[j].m_CtrlPts.size(); i++)
				{
					NVS[j].m_CtrlPts[i].y += distance;
				}
			}
			else if (mode == 3) {
				for (int i = 0; i < NVS[j].m_CtrlPts.size(); i++)
				{
					NVS[j].m_CtrlPts[i].z += distance;
				}
			}
			else if (mode == -1) {
				for (int i = 0; i < NVS[j].m_CtrlPts.size(); i++)
				{
					NVS[j].m_CtrlPts[i].x -= distance;
				}
			}
			else if (mode == -2) {
				for (int i = 0; i < NVS[j].m_CtrlPts.size(); i++)
				{
					NVS[j].m_CtrlPts[i].y -= distance;
				}
			}
			else if (mode == -3) {
				for (int i = 0; i < NVS[j].m_CtrlPts.size(); i++)
				{
					NVS[j].m_CtrlPts[i].z -= distance;
				}
			}
		}
	}

	//������ת��anglesΪ�Ƕȣ�mode123�ֱ�Ϊ��xyz����
	void Rolate(Spline &NL, double angles, int mode)
	{
		if (mode == 1)
		{
			for (int i = 0; i < NL.m_CtrlPts.size(); i++)
			{
				double w = NL.m_CtrlPts[i].w;
				NL.m_CtrlPts[i] = NL.m_CtrlPts[i].RotateX(angles);
				NL.m_CtrlPts[i].w = w;
			}
		}
		if (mode == 2)
		{
			for (int i = 0; i < NL.m_CtrlPts.size(); i++)
			{
				double w = NL.m_CtrlPts[i].w;
				NL.m_CtrlPts[i] = NL.m_CtrlPts[i].RotateY(angles);
				NL.m_CtrlPts[i].w = w;
			}
		}
		if (mode == 3)
		{
			for (int i = 0; i < NL.m_CtrlPts.size(); i++)
			{
				double w = NL.m_CtrlPts[i].w;
				NL.m_CtrlPts[i] = NL.m_CtrlPts[i].RotateZ(angles);
				NL.m_CtrlPts[i].w = w;
			}
		}
	}



	//���浥������ת��������ת����Ӧ����ƽ��֮ǰ��anglesΪ�Ƕȣ�mode123�ֱ�Ϊ��xyz����
	void Rolate(SplineSurface &NS, double angles, int mode)
	{
		if (mode == 1)
		{
			for (int i = 0; i < NS.m_CtrlPts.size(); i++)
			{
				double w = NS.m_CtrlPts[i].w;
				NS.m_CtrlPts[i] = NS.m_CtrlPts[i].RotateX(angles);
				NS.m_CtrlPts[i].w = w;
			}
		}
		if (mode == 2)
		{
			for (int i = 0; i < NS.m_CtrlPts.size(); i++)
			{
				double w = NS.m_CtrlPts[i].w;
				NS.m_CtrlPts[i] = NS.m_CtrlPts[i].RotateY(angles);
				NS.m_CtrlPts[i].w = w;
			}
		}
		if (mode == 3)
		{
			for (int i = 0; i < NS.m_CtrlPts.size(); i++)
			{
				double w = NS.m_CtrlPts[i].w;
				NS.m_CtrlPts[i] = NS.m_CtrlPts[i].RotateZ(angles);
				NS.m_CtrlPts[i].w = w;
			}
		}
	}

	//���浥������ת��������ת����Ӧ����ƽ��֮ǰ��anglesΪ�Ƕȣ�mode123�ֱ�Ϊ��xyz����
	void Rolate(varray<SplineSurface> &NS, double angles, int mode)
	{
		for (int j = 0; j < NS.size(); j++) {
			if (mode == 1)
			{
				for (int i = 0; i < NS[j].m_CtrlPts.size(); i++)
				{
					double w = NS[j].m_CtrlPts[i].w;
					NS[j].m_CtrlPts[i] = NS[j].m_CtrlPts[i].RotateX(angles);
					NS[j].m_CtrlPts[i].w = w;
				}
			}
			if (mode == 2)
			{
				for (int i = 0; i < NS[j].m_CtrlPts.size(); i++)
				{
					double w = NS[j].m_CtrlPts[i].w;
					NS[j].m_CtrlPts[i] = NS[j].m_CtrlPts[i].RotateY(angles);
					NS[j].m_CtrlPts[i].w = w;
				}
			}
			if (mode == 3)
			{
				for (int i = 0; i < NS[j].m_CtrlPts.size(); i++)
				{
					double w = NS[j].m_CtrlPts[i].w;
					NS[j].m_CtrlPts[i] = NS[j].m_CtrlPts[i].RotateZ(angles);
					NS[j].m_CtrlPts[i].w = w;
				}
			}
		}
	}

	void Rolate(varray<SplineVolume>& NV, double angles, int mode) {
		for (int j = 0; j < NV.size(); j++) {
			if (mode == 1)
			{
				for (int i = 0; i < NV[j].m_CtrlPts.size(); i++)
				{
					double w = NV[j].m_CtrlPts[i].w;
					NV[j].m_CtrlPts[i] = NV[j].m_CtrlPts[i].RotateX(angles);
					NV[j].m_CtrlPts[i].w = w;
				}
			}
			if (mode == 2)
			{
				for (int i = 0; i < NV[j].m_CtrlPts.size(); i++)
				{
					double w = NV[j].m_CtrlPts[i].w;
					NV[j].m_CtrlPts[i] = NV[j].m_CtrlPts[i].RotateY(angles);
					NV[j].m_CtrlPts[i].w = w;
				}
			}
			if (mode == 3)
			{
				for (int i = 0; i < NV[j].m_CtrlPts.size(); i++)
				{
					double w = NV[j].m_CtrlPts[i].w;
					NV[j].m_CtrlPts[i] = NV[j].m_CtrlPts[i].RotateZ(angles);
					NV[j].m_CtrlPts[i].w = w;
				}
			}
		}
	}

	
	//Add-Start
	//��λ��
	void MovePoint(Vec3& p, double dis, int mode);
	void MovePoint(Vec4& p, double dis, int mode);
	void MovePoints(varray<Vec3>& ps, double dis, int mode);
	void MovePoints(varray<Vec4>& ps, double dis, int mode);


	//��λ��
	void MoveLine(Spline& nl, double dis, int mode);
	void MoveLine(Feature_Line& fl, double dis, int mode);
	void MoveLine(Spline& nl, const Vec3& p1, const Vec3& p2);
	void MoveLine(Feature_Line& fl, const Vec3& p1, const Vec3& p2);

	void MoveLines(varray<Spline>& nls, double dis, int mode);
	void MoveLines(varray<Feature_Line>& fls, double dis, int mode);
	void MoveLines(varray<Spline>& nls, const Vec3& p1, const Vec3& p2);
	void MoveLines(varray<Feature_Line>& fls, const Vec3& p1, const Vec3& p2);

	//��λ��
	void MoveSurface(SplineSurface& sf, double dis, int mode);
	void MoveSurfaces(varray<SplineSurface>& sfs, double dis, int mode);
	void MoveSurface(SplineSurface& sf, const Vec3& p1, const Vec3& p2);
	void MoveSurfaces(varray<SplineSurface>& sfs, const Vec3& p1, const Vec3& p2);

	void MoveSurface(Feature_Surface& fs, double dis, int mode);
	void MoveSurfaces(varray<Feature_Surface>& fss, double dis, int mode);
	void MoveSurface(Feature_Surface& fs, const Vec3& p1, const Vec3& p2);
	void MoveSurfaces(varray<Feature_Surface>& fss, const Vec3& p1, const Vec3& p2);

	//��λ��
	void MoveVol(SplineVolume& v, double dis, int mode);
	void MoveVols(varray<SplineVolume>& vs, double dis, int mode);

	void MoveVol(Feature_Vol& fv, double dis, int mode);
	void MoveVols(varray<Feature_Vol>& fvs, double dis, int mode);

	//�㾵��
	//p1:ԭʼ��
	//p2��������
	void MirrorPoint(const Vec3& p1, Vec3& p2, int mode);
	void MirrorPoint(const Vec4& p1, Vec4& p2, int mode);

	//һ��㾵��
	void MirrorPoints(const varray<Vec3>& p1, varray<Vec3>& p2, int mode);
	void MirrorPoints(const varray<Vec4>& p1, varray<Vec4>& p2, int mode);

	//���߾���
	void MirrorLine(const Spline& l1, Spline& l2, int mode);
	void MirrorLines(const varray<Spline>& l1, varray<Spline>& l2, int mode);
	void MirrorLine(const Feature_Line& fl1, Feature_Line& fl2, int mode);
	void MirrorLines(const varray<Feature_Line>& fl1, varray<Feature_Line>& fl2, int mode);

	//���澵��
	void MirrorSuface(const SplineSurface& s1, SplineSurface& s2, int mode);
	void MirrorSufaces(const varray<SplineSurface>& s1, varray<SplineSurface>& s2, int mode);
	void MirrorSuface(const Feature_Surface& fs1, Feature_Surface& fs2, int mode);
	void MirrorSufaces(const varray<Feature_Surface>& fs1, varray<Feature_Surface>& fs2, int mode);

	//�徵��
	void MirrorVol(const SplineVolume& v1, SplineVolume& v2, int mode);
	void MirrorVols(const varray<SplineVolume>& v1, varray<SplineVolume>& v2, int mode);
	void MirrorVol(const Feature_Vol& fv1, Feature_Vol& fv2, int mode);
	void MirrorVols(const varray<Feature_Vol>& fv1, varray<Feature_Vol>& fv2, int mode);

	//һ�����ת
	void RotatePoints(varray<Vec3>& ps, double ang, int mode = 3);
	void RotatePoints(varray<Vec4>& ps, double ang, int mode = 3);

	//������ת
	void RotateLine(Spline& nl, double ang, int mode = 3);
	void RotateLine(Feature_Line& fl, double ang, int mode = 3);

	//һ��������ת
	void RotateLines(varray<Spline>& nl, double ang, int mode = 3);
	void RotateLines(varray<Feature_Line>& fl, double ang, int mode = 3);

	//������ת
	void RotateSurface(SplineSurface& ns, double ang, int mode = 3);
	void RotateSurface(Feature_Surface& fs, double ang, int mode = 3);

	//һ��������ת
	void RotateSurfaces(varray<SplineSurface>& ns, double ang, int mode = 3);
	void RotateSurfaces(varray<Feature_Surface>& fs, double ang, int mode = 3);

	//����ת
	void RotateVol(Feature_Vol& fv, double ang, int mode = 3);
	void RotateVols(varray<Feature_Vol>& fvs, double ang, int mode = 3);

	//������
	//num:������Ŀ
	void RoArrayPoint(const Vec3& p, varray<Vec3>& points, int num, int mode = 3);
	void RoArrayPoint(const Vec4& p, varray<Vec4>& points, int num, int mode = 3);

	//������
	void RoArrayLine(const Spline& nl, varray<Spline>& nls, int num, int mode = 3);
	void RoArrayLine(const Feature_Line& fl, varray<Feature_Line>& fls, int num, int mode = 3);
	void RoArrayLines(const varray<Spline>& nl, varray<varray<Spline>>& nls, int num, int mode = 3);
	void RoArrayLines(const varray<Feature_Line>& fl, varray<varray<Feature_Line>>& fls, int num, int mode = 3);

	//������
	void RoArraySurface(const SplineSurface& sf, varray<SplineSurface>& sfs, int num, int mode = 3);
	void RoArraySurface(const Feature_Surface& fs, varray<Feature_Surface>& fss, int num, int mode = 3);
	void RoArraySurfaces(const varray<SplineSurface>& sf, varray<varray<SplineSurface>>& sfs, int num, int mode = 3);
	void RoArraySurfaces(const varray<Feature_Surface>& fs, varray<varray<Feature_Surface>>& fss, int num, int mode = 3);

	//������
	//num:������Ŀ
	void RoArray(varray< SplineVolume>& svol, int num, int mode) {
		varray<SplineVolume> sv = svol;
		double angle = PI * 2.0 / num;
		if (mode == 3) {
			for (int i = 1; i < num; ++i) {
				Rolate(sv, angle, mode);
				for (auto&v : sv) {
					svol.push_back(v);
				}
			}
		}
	}

	void RoArrayVol(const Feature_Vol& fv, varray<Feature_Vol>& fvs, int num, int mode = 3);
	void RoArrayVols(const varray<Feature_Vol>& fv, varray<varray<Feature_Vol>>& fvs, int num, int mode = 3);
	void RoArrayVols(const varray<Feature_Vol>& fv, varray<Feature_Vol>& fvs, int num, int mode = 3);

	//��ά
	//del:�Ƿ�ɾ����ά����
	void DimReduce(varray<varray<Feature_Vol>>& fvs1, varray<Feature_Vol>& fvs2, bool del = false);

	//��ȡ�����߼���������x/y/zֵ����ĳ����������
	varray<Feature_Line> GetFlWithValue(varray<Feature_Line>& fls, double val, int mode);

	//�����߼��ϰ�x/y/z��������
	void SortFeatureLines(varray<Feature_Line>& fls, int mode);


	//COONS��ֵ��������
	void OrderCoonsLines(varray<Spline>& coonsline);

	//��ȡһ��NURBS���ĳһ�������(unfinished)
	varray<SplineSurface> GetSurfsWithVols(const varray<SplineVolume>& vols, int dir);

	//����ȥ��
	void DecoincideNurbsLine(varray<Spline>& nl);

	//�趨�ʷֲ�������
	void SetQuadParameter(const varray<Feature_Surface>& fss, varray<varray<Spline>>& edgelines,
		varray<varray<int>>& seg, varray<bool>& genus, varray<int>& sfnum);

	//��ȡ�ʷ���ɵ�����������
	void GetQuadPartData(varray<Feature_Surface>& fss, SfCtainTreeNode* root);

	//���������漯�ϣ�ִ���ʷֲ���
	void QuadPartFeatureSurfaces(varray<Feature_Surface>& fss, varray<Spline>& conline);
	void QuadPartFeatureSurfaces(Feature_Surface& fs);

	//Add-End





	varray<SplineVolume> CreatSweepVol(const varray<SplineSurface>& NSf, double lengths, int modes)
	{
		varray<SplineVolume> NVS;
		for (int i = 0; i < NSf.size(); i++)
		{
			Creat_Vol CV;
			CV.InitVol(NSf[i], lengths, modes);
			NVS.push_back(CV.NV);
		}
		return NVS;
	}

	varray<SplineVolume> CreatSweepVol_6pts(const varray<SplineSurface>& NSf, double lengths, int modes)
	{
		Vec3 Ori;
		Spline NL;
		varray<SplineVolume> NVS;
		switch (modes) {
		case 1:
			Ori = { 1,0,0 };
			break;
		case -1:
			Ori = { -1,0,0 };
			break;
		case 2:
			Ori = { 0,1,0 };
			break;
		case -2:
			Ori = { 0,-1,0 };
			break;
		case 3:
			Ori = { 0,0,1 };
			break;
		case -3:
			Ori = { 0,0,-1 };
			break;
		}
		Vec3 Path = Ori * lengths;
		Vec3 Path_Position0 = { 0,0,0 };
		Vec3 Path_Position1 = Path / 5;
		Vec3 Path_Position2 = 2 * Path / 5;
		Vec3 Path_Position3 = 3 * Path / 5;
		Vec3 Path_Position4 = 4 * Path / 5;
		Vec3 Path_Position5 = 5 * Path / 5;
		Vec3 Path_Position6 = Path;
		NL.m_CtrlPts.push_back(Path_Position0);
		NL.m_CtrlPts.push_back(Path_Position1);
		NL.m_CtrlPts.push_back(Path_Position2);
		NL.m_CtrlPts.push_back(Path_Position3);
		NL.m_CtrlPts.push_back(Path_Position4);
		NL.m_CtrlPts.push_back(Path_Position5);
		NL.m_Degree = 2;
		NL.m_Knots.push_back(0);
		NL.m_Knots.push_back(0);
		NL.m_Knots.push_back(0);
		NL.m_Knots.push_back(0.25);
		NL.m_Knots.push_back(0.5);
		NL.m_Knots.push_back(0.75);
		NL.m_Knots.push_back(1);
		NL.m_Knots.push_back(1);
		NL.m_Knots.push_back(1);
		for (int i = 0; i < NSf.size(); i++)
		{
			SplineVolume NV;
			NV.CreateTransSweepSplineVolume(NL, NSf[i]);
			NVS.push_back(NV);
		}
		return NVS;
	}

	void KnotRefine(varray<SplineVolume>& new_All,
		varray<double>& uk, varray<double>& vk, varray<double>& wk)
	{
		for (int i = 0; i < new_All.size(); i++)
		{
			new_All[i].KnotsRefine(uk, vk, wk);
		}
	}

	varray<varray<SplineSurface>> GetSurfaces(const varray<SplineVolume>& NVS)			//�����е�����ȡ����
	{
		SplineSurface U0;		//0
		SplineSurface U1;		//1
		SplineSurface V0;		//2
		SplineSurface V1;		//3
		SplineSurface W0;		//4
		SplineSurface W1;		//5
		varray<varray<SplineSurface>> Res;
		varray<SplineSurface> Sum;
		for (int i = 0; i < NVS.size(); i++) {
			int m = NVS[i].m_vNum;
			int n = NVS[i].m_uNum;
			int p = NVS[i].m_wNum;
			SplineVolume Vtemp = NVS[i];
			int idx = 0;
			for (int x = 0; x < m*n; x++)
			{
				W0.m_CtrlPts.push_back(NVS[i].m_CtrlPts[x]);
			}
			W0.SetSurface(Vtemp.m_uDegree, Vtemp.m_vDegree, Vtemp.m_uNum, Vtemp.m_vNum, Vtemp.m_uKnots, Vtemp.m_vKnots);
			for (int x = m * n*(p - 1); x < m*n*p; x++)
			{
				W1.m_CtrlPts.push_back(NVS[i].m_CtrlPts[x]);
			}
			W1.SetSurface(Vtemp.m_uDegree, Vtemp.m_vDegree, Vtemp.m_uNum, Vtemp.m_vNum, Vtemp.m_uKnots, Vtemp.m_vKnots);

			for (int y = 0; y <= (p - 1)*m*n; y += m * n) {
				for (int x = y; x <= y + m * n - m; x += m) {
					V0.m_CtrlPts.push_back(NVS[i].m_CtrlPts[x]);
				}
			}
			V0.SetSurface(Vtemp.m_wDegree, Vtemp.m_uDegree, Vtemp.m_wNum, Vtemp.m_uNum, Vtemp.m_wKnots, Vtemp.m_uKnots);
			for (int y = m - 1; y <= m * n*(p - 1) + m - 1; y += m * n)
				for (int x = y; x <= y + m * n - m; x += m) {
					V1.m_CtrlPts.push_back(NVS[i].m_CtrlPts[x]);
				}
			V1.SetSurface(Vtemp.m_wDegree, Vtemp.m_uDegree, Vtemp.m_wNum, Vtemp.m_uNum, Vtemp.m_wKnots, Vtemp.m_uKnots);

			for (int y = 0; y < p; y++) {
				for (int x = y * m*n; x < y*m*n + m; x++)
				{
					U0.m_CtrlPts.push_back(NVS[i].m_CtrlPts[x]);
				}
			}
			U0.SetSurface(Vtemp.m_wDegree, Vtemp.m_vDegree, Vtemp.m_wNum, Vtemp.m_vNum, Vtemp.m_wKnots, Vtemp.m_vKnots);

			for (int y = 0; y < p; y++) {
				for (int x = y * m*n + m * n - m; x < y*m*n + m * n; x++)
				{
					U1.m_CtrlPts.push_back(NVS[i].m_CtrlPts[x]);
				}
			}
			U1.SetSurface(Vtemp.m_wDegree, Vtemp.m_vDegree, Vtemp.m_wNum, Vtemp.m_vNum, Vtemp.m_wKnots, Vtemp.m_vKnots);


			Sum.push_back(V1);
			V1.m_CtrlPts.clear();
			Sum.push_back(V0);
			V0.m_CtrlPts.clear();

			Sum.push_back(U0);
			U0.m_CtrlPts.clear();

			Sum.push_back(U1);
			U1.m_CtrlPts.clear();

			Sum.push_back(W0);
			W0.m_CtrlPts.clear();

			Sum.push_back(W1);
			W1.m_CtrlPts.clear();

			Res.push_back(Sum);
			Sum.clear();
		}
		return Res;
	}

	varray<SplineSurface> GetSurfaces(SplineVolume NVS)			//�����е�����ȡ����
	{
		SplineSurface U0;		//0
		SplineSurface U1;		//1
		SplineSurface V0;		//2
		SplineSurface V1;		//3
		SplineSurface W0;		//4
		SplineSurface W1;		//5
		varray<varray<SplineSurface>> Res;
		varray<SplineSurface> Sum;
		int m = NVS.m_vNum;
		int n = NVS.m_uNum;
		int p = NVS.m_wNum;
		SplineVolume Vtemp = NVS;
		int idx = 0;
		for (int x = 0; x < m*n; x++)
		{
			W0.m_CtrlPts.push_back(NVS.m_CtrlPts[x]);
		}
		W0.SetSurface(Vtemp.m_uDegree, Vtemp.m_vDegree, Vtemp.m_uNum, Vtemp.m_vNum, Vtemp.m_uKnots, Vtemp.m_vKnots);
		for (int x = m * n*(p - 1); x < m*n*p; x++)
		{
			W1.m_CtrlPts.push_back(NVS.m_CtrlPts[x]);
		}
		W1.SetSurface(Vtemp.m_uDegree, Vtemp.m_vDegree, Vtemp.m_uNum, Vtemp.m_vNum, Vtemp.m_uKnots, Vtemp.m_vKnots);

		for (int y = 0; y <= (p - 1)*m*n; y += m * n) {
			for (int x = y; x <= y + m * n - m; x += m) {
				V0.m_CtrlPts.push_back(NVS.m_CtrlPts[x]);
			}
		}
		V0.SetSurface(Vtemp.m_wDegree, Vtemp.m_uDegree, Vtemp.m_wNum, Vtemp.m_uNum, Vtemp.m_wKnots, Vtemp.m_uKnots);
		for (int y = m - 1; y <= m * n*(p - 1) + m - 1; y += m * n)
			for (int x = y; x <= y + m * n - m; x += m) {
				V1.m_CtrlPts.push_back(NVS.m_CtrlPts[x]);
			}
		V1.SetSurface(Vtemp.m_wDegree, Vtemp.m_uDegree, Vtemp.m_wNum, Vtemp.m_uNum, Vtemp.m_wKnots, Vtemp.m_uKnots);

		for (int y = 0; y < p; y++) {
			for (int x = y * m*n; x < y*m*n + m; x++)
			{
				U0.m_CtrlPts.push_back(NVS.m_CtrlPts[x]);
			}
		}
		U0.SetSurface(Vtemp.m_wDegree, Vtemp.m_vDegree, Vtemp.m_wNum, Vtemp.m_vNum, Vtemp.m_wKnots, Vtemp.m_vKnots);

		for (int y = 0; y < p; y++) {
			for (int x = y * m*n + m * n - m; x < y*m*n + m * n; x++)
			{
				U1.m_CtrlPts.push_back(NVS.m_CtrlPts[x]);
			}
		}
		U1.SetSurface(Vtemp.m_wDegree, Vtemp.m_vDegree, Vtemp.m_wNum, Vtemp.m_vNum, Vtemp.m_wKnots, Vtemp.m_vKnots);


		Sum.push_back(V1);
		V1.m_CtrlPts.clear();
		Sum.push_back(V0);
		V0.m_CtrlPts.clear();

		Sum.push_back(U0);
		U0.m_CtrlPts.clear();

		Sum.push_back(U1);
		U1.m_CtrlPts.clear();

		Sum.push_back(W0);
		W0.m_CtrlPts.clear();

		Sum.push_back(W1);
		W1.m_CtrlPts.clear();

		return Sum;
	}
	varray<varray<Spline>> GetBoundryLines(const varray<SplineSurface>& NSf)
	{
		varray<varray<Spline>> ans;
		varray<Spline> cur_ans;
		Spline L;
		varray<double> k;
		k.push_back(0);
		k.push_back(0);
		k.push_back(0);
		k.push_back(1);
		k.push_back(1);
		k.push_back(1);
		L.m_Degree = 2;
		L.m_Knots = k;
		for (auto NS : NSf) {
			L.m_CtrlPts.clear();
			L.m_CtrlPts.push_back(NS.m_CtrlPts[0]);
			L.m_CtrlPts.push_back(NS.m_CtrlPts[1]);
			L.m_CtrlPts.push_back(NS.m_CtrlPts[2]);
			cur_ans.push_back(L);
			L.m_CtrlPts.clear();
			L.m_CtrlPts.push_back(NS.m_CtrlPts[0]);
			L.m_CtrlPts.push_back(NS.m_CtrlPts[3]);
			L.m_CtrlPts.push_back(NS.m_CtrlPts[6]);
			cur_ans.push_back(L);
			L.m_CtrlPts.clear();
			L.m_CtrlPts.push_back(NS.m_CtrlPts[6]);
			L.m_CtrlPts.push_back(NS.m_CtrlPts[7]);
			L.m_CtrlPts.push_back(NS.m_CtrlPts[8]);
			cur_ans.push_back(L);
			L.m_CtrlPts.clear();
			L.m_CtrlPts.push_back(NS.m_CtrlPts[2]);
			L.m_CtrlPts.push_back(NS.m_CtrlPts[5]);
			L.m_CtrlPts.push_back(NS.m_CtrlPts[8]);
			cur_ans.push_back(L);
			ans.push_back(cur_ans);
			cur_ans.clear();
		}
		return ans;
	}

};

//Բ��
class Cir_arc
{
public:
	Cir_arc() {};
	Cir_arc(double r, double ang, bool seg = true);

	Feature_Line Getarc();
	Feature_Line GetStdarc();				//���Բ����p0��X����
	Feature_Line GetRoarc(double roang);	//���Բ������StdArcΪ������תroang�Ƕ�
private:
	double m_ang;							//��ת�Ƕ�
	Feature_Line m_arc;						//Բ��������
};


//Բ��������
class Circle_2
{
public:
	Circle_2() {}

	//���캯��
	//seg:�߽��߿ɷָ���(Ĭ�ϲ��ɷָ�)
	//empty:�������Ƿ�Ϊ��
	//complete:�������Ƿ���Ҫ�������(ֱ�����NURBS���洴��)
	Circle_2(double r, bool seg = false, bool empty = false, bool complete = true);

	Feature_Surface GetCircle();		//���Բ��������

private:
	Model_Solution ms;
	Feature_Surface m_circle;			//�������������
};

//Բ����������
class Ring
{
public:
	Ring() {}

	//���캯��
	//Ĭ��90��
	Ring(double r, double R, double angle = PI*0.5);

	Feature_Surface GetRing();						//���������

private:
	Model_Solution ms;
	Feature_Surface m_ring;							//�������������	
};

//������Բ��
class Cir_Ring
{
public:
	Cir_Ring() {}

	Cir_Ring(double r, double R);

	Feature_Surface GetCirRing();

private:
	Model_Solution ms;
	Feature_Surface m_cir_ring;

};

//����
class Rectangle_2
{
public:
	Rectangle_2() {}

	//���캯��
	//length:��
	//width:��
	//empty:Ĭ�Ϸǿ�
	//complete:Ĭ�����NURBS���洴��
	Rectangle_2(double length, double width, bool empty = false, bool complete = true);

	//�������������
	Feature_Surface GetRectangle();

	//������������棬���½���ԭ���غ�
	Feature_Surface GetStdRectangle();

private:
	Model_Solution ms;
	Feature_Surface m_rectangle;	//�������������
	double m_length;				//��
	double m_width;					//��
};

//���Ρ�Բ������
class Rec_Cir
{
public:
	Rec_Cir() {}

	//���캯��
	//length:��
	//width:��
	//r:Բ�뾶
	Rec_Cir(double length, double width, double r);

	//������Ρ�Բ������
	Feature_Surface GetRecCir();

	//������Ρ�Բ�����棬���½���ԭ���غ�
	Feature_Surface GetStdRecCir();

private:
	double m_length;				//���γ�
	double m_width;					//���ο�
	Model_Solution ms;
	Feature_Surface m_rec_cir;		//�������������
};



//ƽ��Բ����
//��Բ��Ĭ��Բ�ĵ�Ϊ����ԭ�㣬������x��н�90�ȣ���ʱ�룩
//�Ƕ��û����ƣ���Χ��0��PI��
class Cirle_arc {
public:
	//�ײ�����
	class L_Feature {
	public:
		varray<Vec4> m_Cpts;			//���Ƶ�
		Vec3 p11 = { -7.071,7.071,0 };
		Vec3 p12 = { 7.071,7.071,0 };
		Vec3 p21 = { 0,5,0 };
		Vec3 p22 = { 0,10,0 };
		Vec3 p23 = { 0,15,0 };

		double m_udegree = 2;				//u�������
		double m_vdegree = 2;				//v�������
		int m_uNum = 3;
		int m_vNum = 3;

		varray<double> m_uknocts;			//u�ڵ�ʸ��
		varray<double> m_vknocts;			//v�ڵ�ʸ��

		void InitDegree2(double udegree = 2,
			double vdegree = 2,
			vector<double> uknocts = { 0,0,0,1,1,1 },
			vector<double> vknocts = { 0,0,0,1,1,1 },
			int uNum = 3,
			int vNum = 3)
		{
			this->m_udegree = udegree;
			this->m_vdegree = vdegree;
			this->m_uNum = uNum;
			this->m_vNum = vNum;
			for (int i = 0; i < uknocts.size(); i++)
				this->m_uknocts.push_back(uknocts[i]);
			for (int i = 0; i < vknocts.size(); i++)
				this->m_vknocts.push_back(vknocts[i]);
		}

		void Get2DegreeCtrlPts(Vec3 p1, Vec3 p2, double alph, double r)
		{
			double L = abs(p1.x);

			double yd = r / cos(alph);

			double x = r * tan(alph);

			double wd = abs(L) / x;

			Vec4 ctrpt0 = { 0,yd,0,wd };

			m_Cpts.push_back(p1);
			m_Cpts.push_back(ctrpt0);
			m_Cpts.push_back(p2);
		}


	};

	//�в�����
	class M_Feature {
	public:
		//��������C1
		double m_r1 = 10;							//�뾶
		double m_cangle1 = PI / 2;					//Բ�Ľ�

	//��������C2
		double m_L2 = 10;
		double m_L21 = 5;
		double m_L22 = 5;
	};

	//�߲�����
	class H_Feature {
	public:
		//Vec3 ctpt;								//���ĵ�
		double m_r = 5;								//�ڱ߽�뾶
		double m_R = 15;							//��߽�뾶
		double m_cangle = PI / 2;						//���Ľ�
	};

	L_Feature lfeature;
	M_Feature mfeature;
	H_Feature hfeature;

	SplineSurface Circular;
	Vec3 m_Position;

	double angle;

	//����ӳ��
	//FF1 �߲㵽��������C1��C2

	void f11(double angle)
	{
		mfeature.m_cangle1 = angle;
	}

	void f12(double r)
	{
		double R = hfeature.m_R;

		mfeature.m_r1 = (R + r) / 2;

		mfeature.m_L2 = R - r;

		mfeature.m_L21 = mfeature.m_L2 / 2;

		mfeature.m_L22 = mfeature.m_L21;

	}

	void f13(double R)
	{
		double r = hfeature.m_r;

		mfeature.m_r1 = (R + r) / 2;

		mfeature.m_L2 = R - r;

		mfeature.m_L21 = mfeature.m_L2 / 2;

		mfeature.m_L22 = mfeature.m_L21;
	}


	//����ӳ��
	//FF2 �в㵽�ײ���ӳ��

	void f21(double r1)
	{
		Vec3 p11;
		Vec3 p12;
		Vec3 p22;
		Vec3 p21;
		Vec3 p23;
		double alph = PI / 2 - mfeature.m_cangle1 / 2;

		p22.y = r1;

		p12.x = r1 * cos(alph);
		p12.y = r1 * sin(alph);

		p11.x = -p12.x;
		p11.y = p12.y;

		p21.y = p22.y - mfeature.m_L21;

		p23.y = p22.y + mfeature.m_L22;


		lfeature.p11 = p11;
		lfeature.p12 = p12;
		lfeature.p21 = p21;
		lfeature.p22 = p22;
		lfeature.p23 = p23;
	}

	void f22(double cangle)
	{
		Vec3 p12;
		Vec3 p11;
		double alph = cangle / 2;

		p12.x = mfeature.m_r1*sin(alph);
		p12.y = mfeature.m_r1*cos(alph);

		p11.x = -p12.x;
		p11.y = p12.y;

		lfeature.p11 = p11;
		lfeature.p12 = p12;
	}

	void f23(double L2)
	{
		Vec3 p23;
		Vec3 p21;

		p23.y = lfeature.p22.y + (lfeature.p23.y - lfeature.p22.y) / mfeature.m_L2*L2;

		p21.y = 2 * lfeature.p22.y - p23.y;

		lfeature.p23 = p23;
		lfeature.p21 = p21;
	}

	void f24(double L21)
	{
		Vec3 p23;

		p23.y = lfeature.p22.y + (lfeature.p23.y - lfeature.p22.y) / mfeature.m_L21*L21;

		lfeature.p23 = p23;

	}

	void f25(double L22)
	{
		Vec3 p21;

		p21.y = lfeature.p22.y - (lfeature.p22.y - lfeature.p21.y) / mfeature.m_L22 * L22;

		lfeature.p21 = p21;

	}
	void fH2M()
	{
		f11(hfeature.m_cangle);
		f12(hfeature.m_r);
		f13(hfeature.m_R);
	}

	void fM2L()
	{
		f22(mfeature.m_cangle1);
		f21(mfeature.m_r1);
		f23(mfeature.m_L2);
		f24(mfeature.m_L21);
		f25(mfeature.m_L22);

	}

	void SetPosition(double x, double y, double z,double alph)
	{
		m_Position = { x,y,z };
		double delt = alph - PI / 2;
		for (int i = 0; i < Circular.m_CtrlPts.size(); i++)
		{
			double x0 = Circular.m_CtrlPts[i].x;
			double y0 = Circular.m_CtrlPts[i].y;
			double z0 = Circular.m_CtrlPts[i].z;

			Circular.m_CtrlPts[i].x = x0 * cos(delt) - y0 * sin(delt);
			Circular.m_CtrlPts[i].y = x0 * sin(delt) + y0 * cos(delt);

			Circular.m_CtrlPts[i].x += x;
			Circular.m_CtrlPts[i].y += y;
			Circular.m_CtrlPts[i].z += z;
		}
		angle = alph;
	}

	void Get2DegreeCtrPts() {
		varray<Vec3> CenterPts;
		double alph = PI / 2 - hfeature.m_cangle / 2;
		Vec3 p21 = lfeature.p21;
		Vec3 p22 = lfeature.p22;
		Vec3 p23 = lfeature.p23;
		Vec3 p11 = lfeature.p11;
		Vec3 p12 = lfeature.p12;
		double r = hfeature.m_r;
		double L1 = mfeature.m_L21;
		double L2 = mfeature.m_L22;
		CenterPts.push_back(p21);
		CenterPts.push_back(p22);
		CenterPts.push_back(p23);
		Vec3 ctpt0;
		ctpt0.x = p12.x - L1 * cos(alph);
		ctpt0.y = p12.y - L1 * sin(alph);
		Vec3 ctpt1 = { -1 * ctpt0.x,ctpt0.y,0 };
		lfeature.Get2DegreeCtrlPts(ctpt0, ctpt1, hfeature.m_cangle / 2, hfeature.m_r);
		lfeature.Get2DegreeCtrlPts(lfeature.p12, lfeature.p11, hfeature.m_cangle / 2, mfeature.m_r1);
		Vec3 ctpt2;
		ctpt2.x = p12.x + L2 * cos(alph);
		ctpt2.y = p12.y + L2 * sin(alph);
		Vec3 ctpt3 = { -1 * ctpt2.x,ctpt2.y ,0 };
		lfeature.Get2DegreeCtrlPts(ctpt2, ctpt3, hfeature.m_cangle / 2, hfeature.m_R);

		lfeature.InitDegree2();
	}

	Cirle_arc(double r, double R, double cangle, double x, double y, double z,double alph)
	{
		//���ݲ���
		hfeature.m_r = r;
		hfeature.m_R = R;
		hfeature.m_cangle = cangle;

		//ӳ��
		fH2M();
		fM2L();

		//������Ƭ
		Get2DegreeCtrPts();

		Circular.SetSurface(lfeature.m_udegree, lfeature.m_vdegree,
			lfeature.m_uNum, lfeature.m_vNum,
			lfeature.m_uknocts, lfeature.m_vknocts
		);

		Circular.m_CtrlPts = lfeature.m_Cpts;
		//��λ
		SetPosition(x, y, z,alph);

	}
	Cirle_arc(){}
	SplineSurface getSuf() { return Circular; }
};



//�����α�Բ������ʣ�ಿ��
//C1�е�Ϊ����ԭ��
class Cirle_Triangle
{
public:
	
	class h_feature {
	public:
		double m_L1 = 5;
		double m_L2 = 5;
		double m_H = 5;
		double m_RX = 0;
		double m_RY = 5;
		double m_R = 2.5*sqrt(2);
	};

	h_feature hfeature;
	
	class m_feature {
	public:
		//C1
		double m_L11 = 5;
		double m_L12 = 5;
		double m_H11 = 5;
		//C2
		double m_RX2 = 0;
		double m_RY2 = 5;
		double m_R = 2.5 * sqrt(2);

	};
	m_feature mfeature;

	class l_feature {
	public:
		varray<Spline> NLS;

		Vec4 P11 = {-5,0,0};
		Vec4 P12 = { 0,0,0 };
		Vec4 P13 = { 5,0,0 };
		Vec4 P21 = { -2.5,2.5,0 };
		Vec4 P22 = { 2.5,2.5,0 };
		double m_udegree = 2;				//u�������
		double m_vdegree = 2;				//v�������
		int m_uNum = 3;
		int m_vNum = 3;

		varray<double> m_uknocts;			//u�ڵ�ʸ��
		varray<double> m_vknocts;			//v�ڵ�ʸ��

		void InitDegree2(double udegree = 2,
			double vdegree = 2,
			vector<double> uknocts = { 0,0,0,1,1,1 },
			vector<double> vknocts = { 0,0,0,1,1,1 },
			int uNum = 3,
			int vNum = 3)
		{
			NLS.resize(4);
			this->m_udegree = udegree;
			this->m_vdegree = vdegree;
			this->m_uNum = uNum;
			this->m_vNum = vNum;
			for (int i = 0; i < uknocts.size(); i++)
				this->m_uknocts.push_back(uknocts[i]);
			for (int i = 0; i < vknocts.size(); i++)
				this->m_vknocts.push_back(vknocts[i]);
			for (int i = 0; i < 4; i++)
			{
				NLS[i].m_Degree = udegree;
				NLS[i].m_Knots = this->m_uknocts;
			}

			
		}

		

		void Get2DegreeCtrlPts(Vec3 p1, Vec3 p2, double rx,double ry,double r)
		{
			Vec4 ctrpt0;
			double x1 = p1.x;
			double y1 = p1.y;
			double x2 = p2.x;
			double y2 = p2.y;
			double a1 = 2 * (x1 - rx);
			double b1 = 2 * (y1 - ry);
			double a2 = 2 * (x2 - rx);
			double b2 = 2 * (y2 - ry);
			double k1 = -1 * a1 / b1;
			double k2 = -1 * a2 / b2;
			ctrpt0.x = (k1*x1 - k2 * x2 + y2 - y1) / (k1 - k2);
			ctrpt0.y = k1 * (ctrpt0.x - x1) + y1;
			Vec3 p0 = p2 - p1;
			double CL0 = p0.Magnitude();
			double alph = asin(CL0 / 2 / r);
			ctrpt0.w = cos(alph);
			NLS[2].m_CtrlPts.push_back(p1);
			NLS[2].m_CtrlPts.push_back(ctrpt0);
			NLS[2].m_CtrlPts.push_back(p2);
		}


	};
	l_feature lfeature;
	SplineSurface CircuSquare;
	Vec3 m_Position = {0,0,0};

	double angle = 0;

	void f11(double L1)
	{
		mfeature.m_L11 = L1;

	}

	void f12(double L2)
	{
		mfeature.m_L12 = L2;
	}

	void f13(double H)
	{
		mfeature.m_H11 = H;
	}

	void f14(double RX)
	{
		mfeature.m_RX2 = RX;
	}

	void f15(double RY)
	{
		mfeature.m_RY2 = RY;
	}

	void f16(double R)
	{
		mfeature.m_R = R;
	}

	void f21(double L11)
	{
		Vec4 P11 = lfeature.P11;
		Vec4 P21 = lfeature.P21;
		P11.x = -1 * L11;
		mfeature.m_L11 = L11;
		double x0 = P11.x;
		double y0 = P11.y;
		double xr = mfeature.m_RX2;
		double yr = mfeature.m_RY2;
		double R = mfeature.m_R;
		if (L11 > ERRF) {
			double  k = mfeature.m_H11 / L11;
			double a = 1 + k * k;
			double b = 2 * k*(y0 - yr) - (2 * xr + 2 * x0*k*k);
			double c = xr * xr + k * k*x0*x0 - 2 * k*(y0 - yr)*x0 + (y0 - yr)*(y0 - yr) - R * R;
			double x1 = ((-1 * b) - sqrt(b*b - 4 * a*c)) / (2 * a);
			double x2 = ((-1 * b) + sqrt(b*b - 4 * a*c)) / (2 * a);
			lfeature.P21.x = x1 < x2 ? x1 : x2;
			lfeature.P21.y = k * (lfeature.P21.x - x0) + y0;
			
		}
		else {
			lfeature.P21.x = 0;
			lfeature.P21.y = yr - sqrt(R*R - xr * xr);
		}
		lfeature.P11 = P11;
	}

	void f22(double L12)
	{
		Vec4 P13 = lfeature.P13;
		Vec4 P22 = lfeature.P22;
		P13.x = L12;
		mfeature.m_L12 = L12;
		double x0 = P13.x;
		double y0 = P13.y;
		double xr = mfeature.m_RX2;
		double yr = mfeature.m_RY2;
		double R = mfeature.m_R;
		if (L12 > ERRF) {
			double  k = -1 * mfeature.m_H11 / L12;

			double a = 1 + k * k;
			double b = 2 * k*(y0 - yr) - (2 * xr + 2 * x0*k*k);
			double c = xr * xr + k * k*x0*x0 - 2 * k*(y0 - yr)*x0 + (y0 - yr)*(y0 - yr) - R * R;
			double x1 = ((-1 * b) - sqrt(b*b - 4 * a*c)) / (2 * a);
			double x2 = ((-1 * b) + sqrt(b*b - 4 * a*c)) / (2 * a);
			lfeature.P22.x = x1 > x2 ? x1 : x2;
			lfeature.P22.y = k * (lfeature.P22.x - x0) + y0;
			
		}
		else {
			lfeature.P22.x = 0;
			lfeature.P22.y = yr - sqrt(R*R - xr * xr);
		}
		lfeature.P13 = P13;
	}

	void f23(double H)
	{
		mfeature.m_H11 = H;
		Vec4 P11 = lfeature.P11;
		Vec4 P21 = lfeature.P21;
		double L11 = abs(P11.x);
		double x0 = P11.x;
		double y0 = P11.y;
		double xr = mfeature.m_RX2;
		double yr = mfeature.m_RY2;
		double R = mfeature.m_R;
		if (L11 > ERRF) {
			double  k = mfeature.m_H11 / L11;

			double a = 1 + k * k;
			double b = 2 * k*(y0 - yr) - (2 * xr + 2 * x0*k*k);
			double c = xr * xr + k * k*x0*x0 - 2 * k*(y0 - yr)*x0 + (y0 - yr)*(y0 - yr) - R * R;
			double x1 = ((-1 * b) - sqrt(b*b - 4 * a*c)) / (2 * a);
			double x2 = ((-1 * b) + sqrt(b*b - 4 * a*c)) / (2 * a);
			lfeature.P21.x = x1 < x2 ? x1 : x2;
			lfeature.P21.y = k * (lfeature.P21.x - x0) + y0;
		}
		else
		{
			lfeature.P21.x = 0;
			lfeature.P21.y = yr - sqrt(R*R - xr * xr);
		}
		Vec4 P13 = lfeature.P13;
		Vec4 P22 = lfeature.P22;
		double L12 = P13.x;
		double nx0 = P13.x;
		double ny0 = P13.y;
		double nxr = mfeature.m_RX2;
		double nyr = mfeature.m_RY2;
		double nR = mfeature.m_R;
		if (L12 > ERRF) {
			double  nk = -1 * mfeature.m_H11 / L12;

			double na = 1 + nk * nk;
			double nb = 2 * nk*(ny0 - nyr) - 1 * (2 * nxr + 2 * nx0*nk*nk);
			double nc = nxr * nxr + nk * nk*nx0*nx0 - 2 * nk*(ny0 - nyr)*nx0 + (ny0 - nyr)*(ny0 - nyr) - nR * nR;
			double nx1 = ((-1 * nb) - sqrt(nb*nb - 4 * na*nc)) / (2 * na);
			double nx2 = ((-1 * nb) + sqrt(nb*nb - 4 * na*nc)) / (2 * na);
			lfeature.P22.x = nx1 > nx2 ? nx1 : nx2;
			lfeature.P22.y = nk * (lfeature.P22.x - nx0) + ny0;

		}
		else {
			lfeature.P22.x = 0;
			lfeature.P22.y = nyr - sqrt(nR*nR - nxr * nxr);
		}
	}

	void f24(double RX)
	{
		mfeature.m_RX2 = RX;
		Vec4 P11 = lfeature.P11;
		Vec4 P21 = lfeature.P21;
		double L11 = abs(P11.x);
		double x0 = P11.x;
		double y0 = P11.y;
		double xr = mfeature.m_RX2;
		double yr = mfeature.m_RY2;
		double R = mfeature.m_R;
		if (L11 > ERRF) {
			double  k = mfeature.m_H11 / L11;

			double a = 1 + k * k;
			double b = 2 * k*(y0 - yr) - (2 * xr + 2 * x0*k*k);
			double c = xr * xr + k * k*x0*x0 - 2 * k*(y0 - yr)*x0 + (y0 - yr)*(y0 - yr) - R * R;
			double x1 = ((-1 * b) - sqrt(b*b - 4 * a*c)) / (2 * a);
			double x2 = ((-1 * b) + sqrt(b*b - 4 * a*c)) / (2 * a);
			lfeature.P21.x = x1 < x2 ? x1 : x2;
			lfeature.P21.y = k * (lfeature.P21.x - x0) + y0;
		}
		else
		{
			lfeature.P21.x = 0;
			lfeature.P21.y = yr - sqrt(R*R - xr * xr);
		}
		Vec4 P13 = lfeature.P13;
		Vec4 P22 = lfeature.P22;
		double L12 = P13.x;
		double nx0 = P13.x;
		double ny0 = P13.y;
		double nxr = mfeature.m_RX2;
		double nyr = mfeature.m_RY2;
		double nR = mfeature.m_R;
		if (L12 > ERRF) {
			double  nk = -1 * mfeature.m_H11 / L12;

			double na = 1 + nk * nk;
			double nb = 2 * nk*(ny0 - nyr) - 1 * (2 * nxr + 2 * nx0*nk*nk);
			double nc = nxr * nxr + nk * nk*nx0*nx0 - 2 * nk*(ny0 - nyr)*nx0 + (ny0 - nyr)*(ny0 - nyr) - nR * nR;
			double nx1 = ((-1 * nb) - sqrt(nb*nb - 4 * na*nc)) / (2 * na);
			double nx2 = ((-1 * nb) + sqrt(nb*nb - 4 * na*nc)) / (2 * na);
			lfeature.P22.x = nx1 > nx2 ? nx1 : nx2;
			lfeature.P22.y = nk * (lfeature.P22.x - nx0) + ny0;

		}
		else {
			lfeature.P22.x = 0;
			lfeature.P22.y = nyr - sqrt(nR*nR - nxr * nxr);
		}
	}

	void f25(double RY)
	{
		mfeature.m_RY2 = RY;
		Vec4 P11 = lfeature.P11;
		Vec4 P21 = lfeature.P21;
		double L11 = abs(P11.x);
		double x0 = P11.x;
		double y0 = P11.y;
		double xr = mfeature.m_RX2;
		double yr = mfeature.m_RY2;
		double R = mfeature.m_R;
		if (L11 > ERRF) {
			double  k = mfeature.m_H11 / L11;

			double a = 1 + k * k;
			double b = 2 * k*(y0 - yr) - (2 * xr + 2 * x0*k*k);
			double c = xr * xr + k * k*x0*x0 - 2 * k*(y0 - yr)*x0 + (y0 - yr)*(y0 - yr) - R * R;
			double x1 = ((-1 * b) - sqrt(b*b - 4 * a*c)) / (2 * a);
			double x2 = ((-1 * b) + sqrt(b*b - 4 * a*c)) / (2 * a);
			lfeature.P21.x = x1 < x2 ? x1 : x2;
			lfeature.P21.y = k * (lfeature.P21.x - x0) + y0;
		}
		else
		{
			lfeature.P21.x = 0;
			lfeature.P21.y = yr - sqrt(R*R - xr * xr);
		}
		Vec4 P13 = lfeature.P13;
		Vec4 P22 = lfeature.P22;
		double L12 = P13.x;
		double nx0 = P13.x;
		double ny0 = P13.y;
		double nxr = mfeature.m_RX2;
		double nyr = mfeature.m_RY2;
		double nR = mfeature.m_R;
		if (L12 > ERRF) {
			double  nk = -1 * mfeature.m_H11 / L12;

			double na = 1 + nk * nk;
			double nb = 2 * nk*(ny0 - nyr) - 1 * (2 * nxr + 2 * nx0*nk*nk);
			double nc = nxr * nxr + nk * nk*nx0*nx0 - 2 * nk*(ny0 - nyr)*nx0 + (ny0 - nyr)*(ny0 - nyr) - nR * nR;
			double nx1 = ((-1 * nb) - sqrt(nb*nb - 4 * na*nc)) / (2 * na);
			double nx2 = ((-1 * nb) + sqrt(nb*nb - 4 * na*nc)) / (2 * na);
			lfeature.P22.x = nx1 > nx2 ? nx1 : nx2;
			lfeature.P22.y = nk * (lfeature.P22.x - nx0) + ny0;

		}
		else {
			lfeature.P22.x = 0;
			lfeature.P22.y = nyr - sqrt(nR*nR - nxr * nxr);
		}
	}

	void f26(double R1)
	{
		mfeature.m_R = R1;
		Vec4 P11 = lfeature.P11;
		Vec4 P21 = lfeature.P21;
		double L11 = abs(P11.x);
		double x0 = P11.x;
		double y0 = P11.y;
		double xr = mfeature.m_RX2;
		double yr = mfeature.m_RY2;
		double R = mfeature.m_R;
		if (L11 > ERRF) {
			double  k = mfeature.m_H11 / L11;

			double a = 1 + k * k;
			double b = 2 * k*(y0 - yr) - (2 * xr + 2 * x0*k*k);
			double c = xr * xr + k * k*x0*x0 - 2 * k*(y0 - yr)*x0 + (y0 - yr)*(y0 - yr) - R * R;
			double x1 = ((-1 * b) - sqrt(b*b - 4 * a*c)) / (2 * a);
			double x2 = ((-1 * b) + sqrt(b*b - 4 * a*c)) / (2 * a);
			lfeature.P21.x = x1 < x2 ? x1 : x2;
			lfeature.P21.y = k * (lfeature.P21.x - x0) + y0;
		}
		else
		{
			lfeature.P21.x = 0;
			lfeature.P21.y = yr - sqrt(R*R - xr * xr);
		}
		Vec4 P13 = lfeature.P13;
		Vec4 P22 = lfeature.P22;
		double L12 = P13.x;
		double nx0 = P13.x;
		double ny0 = P13.y;
		double nxr = mfeature.m_RX2;
		double nyr = mfeature.m_RY2;
		double nR = mfeature.m_R;
		if (L12 > ERRF) {
			double  nk = -1 * mfeature.m_H11 / L12;

			double na = 1 + nk * nk;
			double nb = 2 * nk*(ny0 - nyr) - 1 * (2 * nxr + 2 * nx0*nk*nk);
			double nc = nxr * nxr + nk * nk*nx0*nx0 - 2 * nk*(ny0 - nyr)*nx0 + (ny0 - nyr)*(ny0 - nyr) - nR * nR;
			double nx1 = ((-1 * nb) - sqrt(nb*nb - 4 * na*nc)) / (2 * na);
			double nx2 = ((-1 * nb) + sqrt(nb*nb - 4 * na*nc)) / (2 * na);
			lfeature.P22.x = nx1 > nx2 ? nx1 : nx2;
			lfeature.P22.y = nk * (lfeature.P22.x - nx0) + ny0;

		}
		else {
			lfeature.P22.x = 0;
			lfeature.P22.y = nyr - sqrt(nR*nR - nxr * nxr);
		}
	}

	void fH2M()
	{
		f11(hfeature.m_L1);
		f12(hfeature.m_L2);
		f13(hfeature.m_H);
		f14(hfeature.m_RX);
		f15(hfeature.m_RY);
		f16(hfeature.m_R);
	}

	void fM2L()
	{
		f21(mfeature.m_L11);
		f22(mfeature.m_L12);
		f23(mfeature.m_H11);
		f24(mfeature.m_RX2);
		f25(mfeature.m_RY2);
		f26(mfeature.m_R);
	}

	void SetPosition(double x, double y, double z, double alph)
	{
		m_Position = { x,y,z };
		double delt = alph - PI / 2;
		for (int i = 0; i < CircuSquare.m_CtrlPts.size(); i++)
		{
			double x0 = CircuSquare.m_CtrlPts[i].x;
			double y0 = CircuSquare.m_CtrlPts[i].y;
			double z0 = CircuSquare.m_CtrlPts[i].z;

			CircuSquare.m_CtrlPts[i].x = x0 * cos(delt) - y0 * sin(delt);
			CircuSquare.m_CtrlPts[i].y = x0 * sin(delt) + y0 * cos(delt);

			CircuSquare.m_CtrlPts[i].x += x;
			CircuSquare.m_CtrlPts[i].y += y;
			CircuSquare.m_CtrlPts[i].z += z;
		}
		angle = alph;
	}

	void Get2DegreeCtrPts() {

		lfeature.NLS.clear();
		lfeature.m_uknocts.clear();
		lfeature.m_vknocts.clear();
		lfeature.InitDegree2();
		

		lfeature.NLS[0].m_CtrlPts.push_back(lfeature.P11);
		Vec4 pt = (lfeature.P11 + lfeature.P13) / 2;
		lfeature.NLS[0].m_CtrlPts.push_back(pt);
		lfeature.NLS[0].m_CtrlPts.push_back(lfeature.P13);

		Vec4 P2t = (lfeature.P11 + lfeature.P21) / 2;
		lfeature.NLS[1].m_CtrlPts.push_back(lfeature.P11);
		lfeature.NLS[1].m_CtrlPts.push_back(P2t);
		lfeature.NLS[1].m_CtrlPts.push_back(lfeature.P21);

		lfeature.Get2DegreeCtrlPts(lfeature.P21,lfeature.P22,mfeature.m_RX2,mfeature.m_RY2,mfeature.m_R);

		Vec4 P3t = (lfeature.P13 + lfeature.P22) / 2;
		lfeature.NLS[3].m_CtrlPts.push_back(lfeature.P13);
		lfeature.NLS[3].m_CtrlPts.push_back(P3t);
		lfeature.NLS[3].m_CtrlPts.push_back(lfeature.P22);

		
	}

	Cirle_Triangle(double H , double L1 , 
		double L2 , double Rx , double Ry , 
		double R ,
		double x ,double y ,double z ,double alph )
	{
		hfeature.m_H = H;
		hfeature.m_L1 = L1;
		hfeature.m_L2 = L2;
		hfeature.m_R = R;
		hfeature.m_RX = Rx;
		hfeature.m_RY = Ry;

		fH2M();
		fM2L();

		Get2DegreeCtrPts();

		CircuSquare.CoonsInterpolate(lfeature.NLS);

		SetPosition(x, y, z, alph);

	}
	
	Cirle_Triangle(){}

	void Init()
	{
		fH2M();
		fM2L();

		Get2DegreeCtrPts();

		CircuSquare.CoonsInterpolate(lfeature.NLS);

		SetPosition(m_Position.x, m_Position.y, m_Position.z, angle);
	}
	SplineSurface getSuf() { return CircuSquare; }

};


//���α�Բ������ʣ�ಿ��
class Cirle_Rectangle {
public:

	SplineSurface CirRect;
	Vec3 m_Position = { 0,0,0 };
	double angle = 0;

	class h_feature {
	public:
		double m_L1 = 2;
		double m_L2 = 2;
		double m_xr = 0;
		double m_yr = 10;
		double m_R = 2;
	};

	h_feature hfeature;

	class m_feature {
	public:
		double m_L11 = 2;
		double m_L12 = 2;
		double m_xr2 = 0;
		double m_yr2 = 10;
		double m_R = 2;
	};

	m_feature mfeature;

	class l_feature {
	public:
		Vec4 P11 = { -2,0,0 };
		Vec4 P12 = { 2,0,0 };
		Vec4 P21 = { -2,10,0 };
		Vec4 P22 = { 2,10,0 };

		varray<Spline> NLS;
		double m_udegree = 2;				//u�������
		double m_vdegree = 2;				//v�������
		int m_uNum = 3;
		int m_vNum = 3;

		varray<double> m_uknocts;			//u�ڵ�ʸ��
		varray<double> m_vknocts;			//v�ڵ�ʸ��


		//��ʼ��
		void InitDegree2(double udegree = 2,
			double vdegree = 2,
			vector<double> uknocts = { 0,0,0,1,1,1 },
			vector<double> vknocts = { 0,0,0,1,1,1 },
			int uNum = 3,
			int vNum = 3)
		{
			NLS.resize(4);
			this->m_udegree = udegree;
			this->m_vdegree = vdegree;
			this->m_uNum = uNum;
			this->m_vNum = vNum;
			for (int i = 0; i < uknocts.size(); i++)
				this->m_uknocts.push_back(uknocts[i]);
			for (int i = 0; i < vknocts.size(); i++)
				this->m_vknocts.push_back(vknocts[i]);
			for (int i = 0; i < 4; i++)
			{
				NLS[i].m_Degree = udegree;
				NLS[i].m_Knots = this->m_uknocts;
			}
		}

		//��ȡԲ�����Ŀ��Ƶ�
		void Get2DegreeCtrlPts(Vec3 p1, Vec3 p2, double rx, double ry, double r)
		{
			Vec4 ctrpt0;
			Vec3 p0;
			double x1 = p1.x;
			double y1 = p1.y;
			double x2 = p2.x;
			double y2 = p2.y;
			double a1 = 2 * (x1 - rx);
			double b1 = 2 * (y1 - ry);
			double a2 = 2 * (x2 - rx);
			double b2 = 2 * (y2 - ry);
			double k1 = -1 * a1 / b1;
			double k2 = -1 * a2 / b2;
			ctrpt0.x = (k1*x1 - k2 * x2 + y2 - y1) / (k1 - k2);
			ctrpt0.y = k1 * (ctrpt0.x - x1) + y1;
			p0 = p2 - p1;
			double CL0 = p0.Magnitude();
			double alph = asin(CL0 / 2 / r);
			ctrpt0.w = cos(alph);
			NLS[2].m_CtrlPts.push_back(p1);
			NLS[2].m_CtrlPts.push_back(ctrpt0);
			NLS[2].m_CtrlPts.push_back(p2);
		}
	};

	l_feature lfeature;

	

	void f1()
	{
		mfeature.m_L11 = hfeature.m_L1;
		mfeature.m_L12 = hfeature.m_L2;
		mfeature.m_R = hfeature.m_R;
		mfeature.m_xr2 = hfeature.m_xr;
		mfeature.m_yr2 = hfeature.m_yr;
	}


	double CalCirPts(double x) {
		double x0 = mfeature.m_xr2;
		double y0 = mfeature.m_yr2;
		double R = mfeature.m_R;
		double y = y0 - sqrt(R*R - (x - x0)*(x - x0));
		return y;
	}

	void f2()
	{
		lfeature.P11.x = -1 * mfeature.m_L11;
		lfeature.P21.x = -1 * mfeature.m_L11;
		lfeature.P12.x = mfeature.m_L12;
		lfeature.P22.x = mfeature.m_L12;
		double x1 = lfeature.P11.x;
		lfeature.P21.y = CalCirPts(x1);
		double x2 = lfeature.P12.x;
		lfeature.P22.y = CalCirPts(x2);
	}

	void SetPosition(double x, double y, double z, double alph)
	{
		m_Position = { x,y,z };
		double delt = alph - PI / 2;
		for (int i = 0; i < CirRect.m_CtrlPts.size(); i++)
		{
			double x0 = CirRect.m_CtrlPts[i].x;
			double y0 = CirRect.m_CtrlPts[i].y;
			double z0 = CirRect.m_CtrlPts[i].z;

			CirRect.m_CtrlPts[i].x = x0 * cos(delt) - y0 * sin(delt);
			CirRect.m_CtrlPts[i].y = x0 * sin(delt) + y0 * cos(delt);

			CirRect.m_CtrlPts[i].x += x;
			CirRect.m_CtrlPts[i].y += y;
			CirRect.m_CtrlPts[i].z += z;
		}
		angle = alph;
	}

	void Get2DegreeCtrPts() {

		lfeature.NLS.clear();
		lfeature.m_uknocts.clear();
		lfeature.m_vknocts.clear();
		lfeature.InitDegree2();


		Vec4 tp1 = (lfeature.P11 + lfeature.P12) / 2;
		lfeature.NLS[0].m_CtrlPts.push_back(lfeature.P11);
		lfeature.NLS[0].m_CtrlPts.push_back(tp1);
		lfeature.NLS[0].m_CtrlPts.push_back(lfeature.P12);

		Vec4 P2t = (lfeature.P11 + lfeature.P21) / 2;
		lfeature.NLS[1].m_CtrlPts.push_back(lfeature.P11);
		lfeature.NLS[1].m_CtrlPts.push_back(P2t);
		lfeature.NLS[1].m_CtrlPts.push_back(lfeature.P21);

		lfeature.Get2DegreeCtrlPts(lfeature.P21, lfeature.P22, mfeature.m_xr2, mfeature.m_yr2, mfeature.m_R);

		Vec4 P3t = (lfeature.P12 + lfeature.P22) / 2;
		lfeature.NLS[3].m_CtrlPts.push_back(lfeature.P12);
		lfeature.NLS[3].m_CtrlPts.push_back(P3t);
		lfeature.NLS[3].m_CtrlPts.push_back(lfeature.P22);


	}

	Cirle_Rectangle( double L1,
		double L2, double Rx, double Ry,
		double R,
		double x, double y, double z, double alph)
	{
		hfeature.m_L1 = L1;
		hfeature.m_L2 = L2;
		hfeature.m_R = R;
		hfeature.m_xr = Rx;
		hfeature.m_yr = Ry;

		f1();
		f2();

		Get2DegreeCtrPts();

		CirRect.CoonsInterpolate(lfeature.NLS);

		SetPosition(x, y, z, alph);

	}
	Cirle_Rectangle(){}

	void Init()
	{
		f1();
		f2();

		Get2DegreeCtrPts();

		CirRect.CoonsInterpolate(lfeature.NLS);

		SetPosition(m_Position.x, m_Position.y, m_Position.z, angle);
	}
	
	SplineSurface getSuf() { return CirRect; }

};
class RecTangle {
//����ϵԭ��Ϊ�±߽��е�
public:
	double L = 10;
	double H = 20;

	Vec4 P11 = { -5,0,0 };
	Vec4 P12 = { 5,0,0 };
	Vec4 P21 = { -5,20,0 };
	Vec4 P22 = { 5,20,0 };
public:
	SplineSurface Rect;
	Vec3 m_Position;
	double angle;
public:
	void SetPosition(double x, double y, double z, double alph)
	{
		m_Position = { x,y,z };
		angle = alph;
		double delt = PI / 2 - alph;
		for (int i = 0; i < Rect.m_CtrlPts.size(); i++)
		{
			double x0 = Rect.m_CtrlPts[i].x;
			double y0 = Rect.m_CtrlPts[i].y;
			double z0 = Rect.m_CtrlPts[i].z;

			Rect.m_CtrlPts[i].x = x0 * cos(delt) - y0 * sin(delt);
			Rect.m_CtrlPts[i].y = x0 * sin(delt) + y0 * cos(delt);

			Rect.m_CtrlPts[i].x += x;
			Rect.m_CtrlPts[i].y += y;
			Rect.m_CtrlPts[i].z += z;
		}
	}
public:
	RecTangle (double L0,double H0,double x,double y,double z,double alph)
	{
		H = H0;
		L = L0;

		P11 = { -L / 2,0,0 };
		P12 = { L / 2,0,0 };
		P21 = { -L / 2,H,0 };
		P22 = { L / 2,H,0 };

		varray<Spline> NL(4);
		for (int i = 0; i < 4; i++) {
			NL[i].m_Degree = 2;
			NL[i].m_Knots.push_back(0);
			NL[i].m_Knots.push_back(0);
			NL[i].m_Knots.push_back(0);
			NL[i].m_Knots.push_back(1);
			NL[i].m_Knots.push_back(1);
			NL[i].m_Knots.push_back(1);
		}

		Vec4 P1t = (P11 + P12) / 2;
		NL[0].m_CtrlPts.push_back(P11);
		NL[0].m_CtrlPts.push_back(P1t);
		NL[0].m_CtrlPts.push_back(P12);

		Vec4 P2t = (P11 + P21) / 2;
		NL[1].m_CtrlPts.push_back(P11);
		NL[1].m_CtrlPts.push_back(P2t);
		NL[1].m_CtrlPts.push_back(P21);

		Vec4 P3t = (P21 + P22) / 2;
		NL[2].m_CtrlPts.push_back(P21);
		NL[2].m_CtrlPts.push_back(P3t);
		NL[2].m_CtrlPts.push_back(P22);

		Vec4 P4t = (P12 + P22) / 2;
		NL[3].m_CtrlPts.push_back(P12);
		NL[3].m_CtrlPts.push_back(P4t);
		NL[3].m_CtrlPts.push_back(P22);

		Rect.CoonsInterpolate(NL);

		SetPosition(x, y, z, alph);

	}

	RecTangle()
	{
		P11 = { -L / 2,0,0 };
		P12 = { L / 2,0,0 };
		P21 = { -L / 2,H,0 };
		P22 = { L / 2,H,0 };

		varray<Spline> NL(4);
		for (int i = 0; i < 4; i++) {
			NL[i].m_Degree = 2;
			NL[i].m_Knots.push_back(0);
			NL[i].m_Knots.push_back(0);
			NL[i].m_Knots.push_back(0);
			NL[i].m_Knots.push_back(1);
			NL[i].m_Knots.push_back(1);
			NL[i].m_Knots.push_back(1);
		}

		Vec4 P1t = (P11 + P12) / 2;
		NL[0].m_CtrlPts.push_back(P11);
		NL[0].m_CtrlPts.push_back(P1t);
		NL[0].m_CtrlPts.push_back(P12);

		Vec4 P2t = (P11 + P21) / 2;
		NL[1].m_CtrlPts.push_back(P11);
		NL[1].m_CtrlPts.push_back(P2t);
		NL[1].m_CtrlPts.push_back(P21);

		Vec4 P3t = (P21 + P22) / 2;
		NL[2].m_CtrlPts.push_back(P21);
		NL[2].m_CtrlPts.push_back(P3t);
		NL[2].m_CtrlPts.push_back(P22);

		Vec4 P4t = (P12 + P22) / 2;
		NL[3].m_CtrlPts.push_back(P12);
		NL[3].m_CtrlPts.push_back(P4t);
		NL[3].m_CtrlPts.push_back(P22);

		Rect.CoonsInterpolate(NL);

		SetPosition(m_Position.x, m_Position.y, m_Position.z, angle);
	}

	SplineSurface getSuf() { return Rect; }
};

//�������Σ��е�Ϊ�µ�������
class trapezoid {
private:
	
	double m_upL;
	double m_btL;
	double m_H;
	Vec3 m_Position;
	double angle;
	SplineSurface suf_trapeziod;
private:
	void SetPosition()
	{
		for (auto& pt : suf_trapeziod.m_CtrlPts)
		{
			double delt = PI / 2 - angle;
			double x0 = pt.x;
			double y0 = pt.y;
			double z0 = pt.z;

			pt.x = x0 * cos(delt) - y0 * sin(delt);
			pt.y = x0 * sin(delt) + y0 * cos(delt);

			pt.x += m_Position.x;
			pt.y += m_Position.y;
			pt.z += m_Position.z;
		}
	}
public:
	trapezoid(double upL,double btL,double H,double x, double y,double z, double alph):
		m_upL(upL), m_btL(btL), m_H(H), m_Position{x,y,z},angle(alph)
	{
		varray<double> k;
		k.push_back(0);
		k.push_back(0);
		k.push_back(0);
		k.push_back(1);
		k.push_back(1);
		k.push_back(1);

		suf_trapeziod.SetSurface(2, 2, 3, 3, k, k);
		
		Vec4 pt = { -1*m_btL / 2,0,0,1 };
		suf_trapeziod.m_CtrlPts.push_back(pt);
		pt.x = 0;
		suf_trapeziod.m_CtrlPts.push_back(pt);
		pt.x = m_btL/2;
		suf_trapeziod.m_CtrlPts.push_back(pt);
		pt.x = -(m_btL/2 + m_upL/2) / 2;
		pt.y = m_H / 2;
		suf_trapeziod.m_CtrlPts.push_back(pt);
		pt.x = 0;
		suf_trapeziod.m_CtrlPts.push_back(pt);
		pt.x = (m_btL/2 + m_upL/2) / 2;
		suf_trapeziod.m_CtrlPts.push_back(pt);
		pt.x = -m_upL / 2;
		pt.y = m_H;
		suf_trapeziod.m_CtrlPts.push_back(pt);
		pt.x = 0;
		suf_trapeziod.m_CtrlPts.push_back(pt);
		pt.x = m_upL / 2;
		suf_trapeziod.m_CtrlPts.push_back(pt);

		SetPosition();
	}
	SplineSurface getSuf() { return suf_trapeziod; }


};

class Shape {
	//varray<SplineVolume> Shape_vols;
	static const int mode = 0;
public:
	varray<SplineVolume> virtual getVols() = 0;
	static int getmode() { return  mode; }
};

//�������
class BearingBlock{
private:

	RWGeometric RW;
	double		L = 10;			//�ײ���
	double		H = 20;			//�ײ���
	double		t_base = 2;			//�ײ���
	double		R = 1;			//�ײ�Բ�װ뾶
	int			mode_base = -3;			//�ײ�����
	double		H1 = 10;			//Բ�ĵ��ײ��ľ���
	double		R1 = 4;			//��пװ뾶
	double		r1 = 2;			//��п��ھ�
	double		t_back = 2;			//��Բ����
	double		t_front = 6;			//ǰԲ����	

	string path;
	varray<SplineVolume>		All;
	varray<SplineVolume>		new_All;
	varray<SplineSurface>	NSf;
	varray<Cirle_Triangle>	Css0;
	varray<Cirle_Triangle>	Css1;
	varray<Cirle_Rectangle>	CRs;
	varray<Cirle_arc>		Carcs;
	varray<RecTangle>		Rects;
	varray<double>			L0s;
	varray<double>			H0s;
	Model_Solution m;
	
private:
	varray<double> newUk;
	varray<double> newVk;
	//ѡ�����εĳ��Ϳ�
	void SectRectLH()
	{
		RecTangle* ptr;
		for (int i = 0; i < Rects.size(); i++)
		{
			ptr = &Rects[i];
			L0s.push_back(ptr->L);
			H0s.push_back(ptr->H);
		}
	}

	//�趨���ο�λ��
	void SetPositionRect()
	{
		double h0 = Rects[2].H;
		Rects[0].m_Position.x = Rects[0].L / 2;
		Rects[0].m_Position.y = h0;
		Rects[1].m_Position.x = Rects[0].L + Rects[1].L / 2;
		Rects[1].m_Position.y = h0;
		Rects[2].m_Position.x = Rects[0].L + Rects[2].L / 2;
		Rects[2].m_Position.y = 0;
		Rects[3].m_Position.x = Rects[0].L + Rects[1].L + Rects[3].L / 2;
		Rects[3].m_Position.y = h0;
		Rects[4].m_Position.x = Rects[3].m_Position.x;
		Rects[4].m_Position.y = 0;
		Rects[5].m_Position.x = Rects[0].L + Rects[1].L + Rects[3].L + Rects[5].L / 2;
		Rects[5].m_Position.y = h0;
		Rects[6].m_Position.x = Rects[5].m_Position.x;
		Rects[6].m_Position.y = 0;
		Rects[7].m_Position.x = Rects[0].L + Rects[1].L + Rects[3].L + Rects[5].L + Rects[7].L / 2;
		Rects[7].m_Position.y = h0;
		for (int i = 0; i < Rects.size(); i++)
			Rects[i] = RecTangle();
	}

	//�趨���ǻ���λ��
	void SetPositionCss0()
	{
		double L0 = L0s[0] + L0s[1] + L0s[3] + L0s[5] + L0s[7];
		Css0[0].m_Position.x = Css0[0].hfeature.m_L1;
		Css0[0].m_Position.y = 0;
		Css0[1].m_Position.x = Css0[0].hfeature.m_L1 + Css0[0].hfeature.m_L2;
		Css0[1].m_Position.y = Css0[0].hfeature.m_H;
		Css0[2].m_Position.x = Css0[0].hfeature.m_L1;
		Css0[2].m_Position.y = 2 * Css0[0].hfeature.m_H;
		Css0[3].m_Position.x = 0;
		Css0[3].m_Position.y = Css0[0].hfeature.m_H;
		double delt = L - Css0[0].hfeature.m_L1 - Css0[0].hfeature.m_L2;
		for (int i = 4; i < 8; i++)
		{
			Css0[i] = Css0[i - 4];
			Css0[i].m_Position.x += delt;
		}

		for (int i = 0; i < Css0.size(); i++)
			Css0[i].Init();
	}

	//���γ��Ϳ����ӳ��
	void Map2LH()
	{
		RecTangle* ptr;
		varray<SplineVolume> tempV;
		varray<SplineSurface> tempS;
		double L0 = L0s[0] + L0s[1] + L0s[3] + L0s[5] + L0s[7];
		double H0 = H0s[1] + H0s[2];
		double kl = L / L0;
		double kh = H / H0;
		for (int i = 0; i < Rects.size(); i++)
		{
			ptr = &Rects[i];
			ptr->H = ptr->H * kh;
			ptr->L = ptr->L * kl;
		}
		SetPositionRect();
		for (int i = 0; i < Rects.size(); i++)
		{
			tempS.push_back(Rects[i].Rect);
		}

		tempV = m.CreatSweepVol(tempS, t_base, mode_base);

		for (int i = 0; i < tempV.size(); i++)
			new_All.push_back(tempV[i]);
	}


	void Map2Cs()
	{
		Cirle_Triangle* Csptr = nullptr;
		varray<SplineVolume> tempV;
		varray<SplineSurface> tempS;

		double L0 = L0s[0] + L0s[1] + L0s[3] + L0s[5] + L0s[7];
		double H0 = H0s[1] + H0s[2];
		double kl = L / L0;
		double kh = H / H0;

		for (int i = 0; i < Css0.size(); i++)
		{
			Csptr = &Css0[i];
			if (i % 2 == 0)
			{
				Csptr->hfeature.m_L1 = Csptr->hfeature.m_L1*kl;
				Csptr->hfeature.m_L2 = Csptr->hfeature.m_L2 *kl;
				Csptr->hfeature.m_H = Csptr->hfeature.m_H*kh;
			}
			else
			{
				Csptr->hfeature.m_L1 = Csptr->hfeature.m_L1 *kh;
				Csptr->hfeature.m_L2 = Csptr->hfeature.m_L2 *kh;
				Csptr->hfeature.m_H = Csptr->hfeature.m_H*kl;

			}
			Csptr->hfeature.m_RX = 0;
			Csptr->hfeature.m_RY = Csptr->hfeature.m_H;
			Csptr->hfeature.m_R = R;
		}
		SetPositionCss0();
		for (int i = 0; i < Css0.size(); i++)
		{
			tempS.push_back(Css0[i].CircuSquare);
		}
		tempV = m.CreatSweepVol(tempS, this->t_base, this->mode_base);
		for (int i = 0; i < tempV.size(); i++)
			new_All.push_back(tempV[i]);
	}

	void Map2CRs()
	{
		varray<SplineVolume> tempV;
		double hk = H1 * Css1[0].hfeature.m_H / Css1[0].hfeature.m_RY;
		double z = Rects[0].H;
		int modez = y_forward;
		Css1[0].hfeature.m_L1 = Rects[1].L;
		Css1[0].hfeature.m_L2 = 0;
		Css1[0].m_Position.x = Rects[0].L + Rects[1].L;
		Css1[0].m_Position.y = Rects[2].H;
		Css1[0].hfeature.m_R = R1;
		Css1[0].hfeature.m_RX = Rects[3].L / 2;
		Css1[0].hfeature.m_RY = H1;

		double lh = Css1[0].hfeature.m_L1 + Rects[3].L / 2;
		double k = R1 / lh;
		double a = 1 - k * k;
		double b = -2 * k*k*H1;
		double c = -1 * (k*k*H1*H1 + R1 * R1);
		double x1 = (-1 * b + sqrt(b*b - 4 * a*c)) / (2 * a);
		double x2 = (-1 * b - sqrt(b*b - 4 * a*c)) / (2 * a);
		x1 = x1 > x2 ? x1 : x2;

		double hx = Css1[0].hfeature.m_L1 / lh * (x1 + H1);

		Css1[0].hfeature.m_H = hx;

		Css1[0].Init();
		m.Trans(Css1[0].CircuSquare, -1 * Rects[2].H, y_forward);
		m.Rolate(Css1[0].CircuSquare, PI / 2, 1);
		m.Trans(Css1[0].CircuSquare, 1 * Rects[2].H, y_forward);

		Css1[1].hfeature.m_L1 = 0;
		Css1[1].hfeature.m_L2 = Rects[5].L;
		Css1[1].hfeature.m_H = hx;
		Css1[1].m_Position.x = Rects[0].L + Rects[1].L + Rects[3].L;
		Css1[1].m_Position.y = Rects[2].H;
		Css1[1].hfeature.m_R = R1;
		Css1[1].hfeature.m_RX = -1 * Rects[3].L / 2;
		Css1[1].hfeature.m_RY = H1;
		Css1[1].Init();
		m.Trans(Css1[1].CircuSquare, -1 * Rects[2].H, y_forward);
		m.Rolate(Css1[1].CircuSquare, PI / 2, 1);
		m.Trans(Css1[1].CircuSquare, 1 * Rects[2].H, y_forward);
		varray<SplineSurface> tempS;
		tempS.push_back(Css1[0].CircuSquare);
		tempS.push_back(Css1[1].CircuSquare);
		tempV = m.CreatSweepVol(tempS, z, modez);
		new_All.push_back(tempV[0]);
		new_All.push_back(tempV[1]);


		CRs[0].hfeature.m_L1 = Rects[3].L / 2;
		CRs[0].hfeature.m_L2 = Rects[3].L / 2;
		CRs[0].hfeature.m_R = R1;
		CRs[0].hfeature.m_xr = 0;
		CRs[0].hfeature.m_yr = H1;
		CRs[0].m_Position = Rects[3].m_Position;
		CRs[0].Init();
		m.Trans(CRs[0].CirRect, -1 * Rects[2].H, y_forward);
		m.Rolate(CRs[0].CirRect, PI / 2, 1);
		m.Trans(CRs[0].CirRect, 1 * Rects[2].H, y_forward);
		tempS.clear();
		tempV.clear();
		tempS.push_back(CRs[0].CirRect);
		tempV = m.CreatSweepVol(tempS, z, modez);
		new_All.push_back(tempV[0]);
	}

	void Map2Cas()
	{
		varray<SplineSurface> tempS;
		varray<SplineVolume> tempV;
		double Rx0 = Rects[0].L + Rects[1].L + Rects[3].L / 2;
		double Ry0 = Rects[2].H;
		double tz = Rects[3].H;
		double angle0 = 2 * asin(Rects[3].L / 2 / R1);

		Vec3 lss = Css1[0].lfeature.P21 - Css1[0].lfeature.P22;

		double Lss = lss.Magnitude();

		double angle1 = 2 * asin(Lss / 2 / R1);
		double delt = asin(Rects[3].L / 2 / R1) + angle1 / 2;

		double angle2 = angle1;

		double angle3 = 2 * PI - angle0 - angle1 - angle2;

		Cirle_arc cir = Cirle_arc(r1, R1, angle3, 0, 0, 0, PI / 2);
		Carcs.push_back(cir);
		cir = Cirle_arc(r1, R1, angle1, 0, 0, 0, delt - PI / 2);
		Carcs.push_back(cir);
		cir = Cirle_arc(r1, R1, angle0, 0, 0, 0, -PI / 2);
		Carcs.push_back(cir);
		cir = Cirle_arc(r1, R1, angle2, 0, 0, 0, -delt - PI / 2);
		Carcs.push_back(cir);
		for (int i = 0; i < Carcs.size(); i++)
		{
			m.Rolate(Carcs[i].Circular, PI / 2, 1);
			m.Trans(Carcs[i].Circular, Rx0, 1);
			m.Trans(Carcs[i].Circular, Ry0, 2);
			m.Trans(Carcs[i].Circular, H1, 3);
			tempS.push_back(Carcs[i].Circular);

		}
		tempV = m.CreatSweepVol(tempS, tz, y_forward);
		for (int i = 0; i < tempV.size(); i++)
			new_All.push_back(tempV[i]);

		tempS.clear();
		tempV.clear();
		cir = Cirle_arc(r1, R1, angle3, 0, 0, 0, PI / 2);
		Carcs.push_back(cir);
		cir = Cirle_arc(r1, R1, angle1, 0, 0, 0, delt - PI / 2);
		Carcs.push_back(cir);

		cir = Cirle_arc(r1, R1, angle0, 0, 0, 0, -PI / 2);
		Carcs.push_back(cir);

		cir = Cirle_arc(r1, R1, angle2, 0, 0, 0, -delt - PI / 2);
		Carcs.push_back(cir);
		for (int i = 4; i < Carcs.size(); i++)
		{
			m.Rolate(Carcs[i].Circular, PI / 2, 1);
			m.Trans(Carcs[i].Circular, Rx0, 1);
			m.Trans(Carcs[i].Circular, Ry0 + Rects[1].H, 2);
			m.Trans(Carcs[i].Circular, H1, 3);
			tempS.push_back(Carcs[i].Circular);

		}
		tempV = m.CreatSweepVol(tempS, t_back, y_forward);
		for (int i = 0; i < tempV.size(); i++)
			new_All.push_back(tempV[i]);

		tempS.clear();
		tempV.clear();
		cir = Cirle_arc(r1, R1, angle3, 0, 0, 0, PI / 2);
		Carcs.push_back(cir);

		cir = Cirle_arc(r1, R1, angle1, 0, 0, 0, delt - PI / 2);
		Carcs.push_back(cir);

		cir = Cirle_arc(r1, R1, angle0, 0, 0, 0, -PI / 2);
		Carcs.push_back(cir);

		cir = Cirle_arc(r1, R1, angle2, 0, 0, 0, -delt - PI / 2);
		Carcs.push_back(cir);

		for (int i = 8; i < Carcs.size(); i++)
		{
			m.Rolate(Carcs[i].Circular, PI / 2, 1);
			m.Trans(Carcs[i].Circular, Rx0, 1);
			m.Trans(Carcs[i].Circular, Ry0, 2);
			m.Trans(Carcs[i].Circular, H1, 3);
			tempS.push_back(Carcs[i].Circular);

		}
		tempV = m.CreatSweepVol(tempS, t_front, y_backward);
		for (int i = 0; i < tempV.size(); i++)
			new_All.push_back(tempV[i]);

	}

	//����ӳ��
	void Map2Rib()
	{
		SplineVolume NVl;
		SplineSurface NSf;
		varray<SplineSurface> n1;
		Cirle_Rectangle CRt0;
		CRt0 = CRs[0];
		varray<Spline> NLS;
		NLS = CRt0.lfeature.NLS;
		//��ת
		for (int i = 0; i < NLS.size(); i++)
		{
			for (int j = 0; j < NLS[i].m_CtrlPts.size(); j++)
			{
				double w = NLS[i].m_CtrlPts[j].w;
				Vec4 p = NLS[i].m_CtrlPts[j].RotateX(PI / 2);
				p.w = w;
				NLS[i].m_CtrlPts[j] = p;
			}
		}

		for (int i = 0; i < 3; i++)
		{
			NLS[0].m_CtrlPts[i].x += Rects[1].L + Rects[0].L + Rects[3].L / 2;
			NLS[2].m_CtrlPts[i].x += Rects[1].L + Rects[0].L + Rects[3].L / 2;
			NLS[2].m_CtrlPts[i].y += Rects[2].H - t_front;
		}

		NLS[1].m_CtrlPts.clear();
		Vec4 pt1 = (NLS[0].m_CtrlPts[0] + NLS[2].m_CtrlPts[0]) / 2;
		NLS[1].m_CtrlPts.push_back(NLS[0].m_CtrlPts[0]);
		NLS[1].m_CtrlPts.push_back(pt1);
		NLS[1].m_CtrlPts.push_back(NLS[2].m_CtrlPts[0]);

		NLS[3].m_CtrlPts.clear();
		Vec4 pt2 = (NLS[0].m_CtrlPts[2] + NLS[2].m_CtrlPts[2]) / 2;
		NLS[3].m_CtrlPts.push_back(NLS[0].m_CtrlPts[2]);
		NLS[3].m_CtrlPts.push_back(pt2);
		NLS[3].m_CtrlPts.push_back(NLS[2].m_CtrlPts[2]);

		NSf.CoonsInterpolate(NLS);
		/*varray<SplineSurface> NSfs;
		NSfs.push_back(NSf);
		RW.WriteSplineSurface("D:\\r\\������.txt", NSfs);*/
		n1.push_back(CRt0.CirRect);
		n1.push_back(NSf);

		Vec3 p1 = { 0,0,0 };
		Vec3 p2 = { 0,-0.5,0 };
		Vec3 p3 = { 0,-1,0 };

		Spline NL;
		NL = NLS[0];
		NL.m_CtrlPts.clear();
		NL.m_CtrlPts.push_back(p1);
		NL.m_CtrlPts.push_back(p2);
		NL.m_CtrlPts.push_back(p3);


		NVl.LoftingSplineVolume(NL, n1);
		NVl.m_CtrlPts.clear();
		NVl.m_wNum = 3;
		NVl.m_wKnots.clear();
		NVl.m_wKnots.push_back(0);
		NVl.m_wKnots.push_back(0);
		NVl.m_wKnots.push_back(0);
		NVl.m_wKnots.push_back(1);
		NVl.m_wKnots.push_back(1);
		NVl.m_wKnots.push_back(1);

		for (int i = 0; i < CRt0.CirRect.m_CtrlPts.size(); i++)
		{
			NVl.m_CtrlPts.push_back(CRt0.CirRect.m_CtrlPts[i]);
		}
		for (int i = 0; i < CRt0.CirRect.m_CtrlPts.size(); i++)
		{
			NVl.m_CtrlPts.push_back((CRt0.CirRect.m_CtrlPts[i]+ NSf.m_CtrlPts[i])/2);
		}
		for (int i = 0; i < CRt0.CirRect.m_CtrlPts.size(); i++)
		{
			NVl.m_CtrlPts.push_back(NSf.m_CtrlPts[i]);
		}


		new_All.push_back(NVl);
	}

	void InitCRs()
	{
		varray<SplineSurface> NS2;
		Cirle_Rectangle CRs0;
		Cirle_Triangle Cs0, Cs1;
		Cs0 = Cirle_Triangle(23.077, 5, 0, 1.5, 15, 3.18, -1.5, 0, 0, PI / 2);
		Cs1 = Cirle_Triangle(23.077, 0, 5, -1.5, 15, 3.18, 1.5, 0, 0, PI / 2);
		CRs0 = Cirle_Rectangle(1.5, 1.5, 0, 15, 3.18, 0, 0, 0, PI / 2);
		NS2.push_back(Cs0.CircuSquare);
		NS2.push_back(Cs1.CircuSquare);
		NS2.push_back(CRs0.CirRect);
		Css1.push_back(Cs0);
		Css1.push_back(Cs1);
		CRs.push_back(CRs0);
		m.Rolate(NS2, PI / 2, 1);
		varray<SplineVolume> NVS1 = m.CreatSweepVol(NS2, 5, 2);
		m.Trans(NVS1, 16.5, x_forward);
		m.Trans(NVS1, 20, y_forward);
		for (int i = 0; i < NVS1.size(); i++)
			All.push_back(NVS1[i]);

	}

	void InitRect()
	{
		varray<SplineVolume> NVS;
		NVS.clear();
		varray<RecTangle> Rect(8);

		Rect[0] = RecTangle(10, 5, 5, 20, 0, PI / 2);
		Rect[1] = RecTangle(5, 5, 12.5, 20, 0, PI / 2);
		Rect[2] = RecTangle(5, 20, 12.5, 0, 0, PI / 2);
		Rect[3] = RecTangle(3, 5, 16.5, 20, 0, PI / 2);
		Rect[4] = RecTangle(3, 20, 16.5, 0, 0, PI / 2);
		Rect[5] = RecTangle(5, 5, 20.5, 20, 0, PI / 2);
		Rect[6] = RecTangle(5, 20, 20.5, 0, 0, PI / 2);
		Rect[7] = RecTangle(10, 5, 28, 20, 0, PI / 2);
		Rects = Rect;
		varray<SplineSurface> NSS;
		for (int i = 0; i < Rect.size(); i++)
		{
			NSS.push_back(Rect[i].Rect);
		}
		NVS = m.CreatSweepVol(NSS, t_base, 3);
		RWGeometric RW;
		//RW.WriteSplineVolume("D:\\r\\77777.txt", NVS);
		for (int i = 0; i < NVS.size(); i++)
			All.push_back(NVS[i]);

		//ѡ�����Ϳ����ӳ��

		SectRectLH();


		Cirle_Triangle C0, C1, C2, C3;
		C0 = Cirle_Triangle(10, 5, 5, 0, 10, 3, 0, 0, 0, PI / 2);
		C1 = Cirle_Triangle(5, 10, 10, 0, 5, 3, 5, 10, 0, PI);
		C2 = Cirle_Triangle(10, 5, 5, 0, 10, 3, 0, 20, 0, -PI / 2);
		C3 = Cirle_Triangle(5, 10, 10, 0, 5, 3, -5, 10, 0, 0);
		Css0.push_back(C0);
		Css0.push_back(C1);
		Css0.push_back(C2);
		Css0.push_back(C3);
		for (int i = 0; i < Css0.size(); i++)
		{
			NSf.push_back(Css0[i].CircuSquare);
			Css0[i].m_Position.x += 5;
		}
		m.Trans(NSf, 5, x_forward);

		varray<SplineSurface> NS1 = NSf;

		m.Trans(NS1, 23, x_forward);
		varray<Cirle_Triangle> temp = Css0;
		for (int i = 0; i < NS1.size(); i++)
		{
			temp[i].m_Position.x += 23;
			Css0.push_back(temp[i]);
		}
		temp.clear();
		NVS = m.CreatSweepVol(NSf, t_base, 3);
		for (int i = 0; i < NVS.size(); i++)
		{
			All.push_back(NVS[i]);
		}
		NVS = m.CreatSweepVol(NS1, t_base, 3);
		for (int i = 0; i < NVS.size(); i++)
		{
			All.push_back(NVS[i]);
		}
	}

	void WriteSplineVolumes(string path)
	{
		RWGeometric RW;
		RW.WriteSplineVolume(path, All);
	}

	void WriteNewSplineVolumes(string path)
	{
		RWGeometric RW;
		RW.WriteSplineVolume(path, new_All);
	}
public:
	BearingBlock(double L = 18, double H = 12, double R = 1, double t = 1,
				 double H1 = 6, double R1 = 2, double r1 = 1, double t_b = 1, 
				 double t_f = 5):
		L(L),H(H),R(R),t_base(t),H1(H1),R1(R1),r1(r1),t_back(t_b),t_front(t_f)
	{
		InitRect();
		InitCRs();
		Map2Cs();
		Map2LH();
		Map2CRs();
		Map2Cas();
		Map2Rib();
		AdjSortOfVol();
	}
	varray<SplineVolume>getVols(){ return new_All; }

private:
	void AdjSortOfVol()
	{
		for (int i = 0; i < new_All.size(); i++)
		{
			if (i == 16 || i == 17 || i == 18 || i == 31)
			{
				double n = new_All[i].m_uNum;
				double m = new_All[i].m_vNum;
				double l = new_All[i].m_wNum;
				varray<Vec4> temp = new_All[i].m_CtrlPts;
				new_All[i].m_CtrlPts.clear();
				for (int i0 = n - 1; i0 >= 0; i0--)
				{
					for (int j = 0; j < m; j++) {
						for (int k = 0; k < l; k++)
						{
							new_All[i].m_CtrlPts.push_back(temp[i0 + j * m + k * l*m]);
						}
					}
				}
				if (i == 31)
				{
					varray<Vec4> temp = new_All[i].m_CtrlPts;
					new_All[i].m_CtrlPts.clear();
					for (int k = 0; k < n*m*l; k +=3)
					{
						new_All[i].m_CtrlPts.push_back(temp[k + 2]);
						new_All[i].m_CtrlPts.push_back(temp[k + 1]);
						new_All[i].m_CtrlPts.push_back(temp[k]);
					}
				}
			}
			if (i >= 19 && i <= 30)
			{
				double n = new_All[i].m_uNum;
				double m = new_All[i].m_vNum;
				double l = new_All[i].m_wNum;
				varray<Vec4> temp = new_All[i].m_CtrlPts;
				new_All[i].m_CtrlPts.clear();
				for (int i0 = 0; i0 < n*m; i0 += m)
				{
					for (int k = 0; k < n; k++) {
						for (int j = 0; j < l; j++)
						{
							new_All[i].m_CtrlPts.push_back(temp[i0 + j * n*m + k]);
						}
					}
				}
			}


		}
	}

	void KnotRefine()
	{
		varray<double> k;
		k.push_back(0.25);
		k.push_back(0.25);
		k.push_back(0.5);
		k.push_back(0.5);
		k.push_back(0.75);
		k.push_back(0.75);

		for (int i = 0; i < new_All.size(); i++)
		{
			new_All[i].KnotsRefine(k, k, k);
		}
	}

	//double		L = 10;			//�ײ���
	//double		H = 20;			//�ײ���
	//double		t_base = 2;			//�ײ���
	//double		R = 1;			//�ײ�Բ�װ뾶
	//int			mode_base = -3;			//�ײ�����
	//double		H1 = 10;			//Բ�ĵ��ײ��ľ���
	//double		R1 = 4;			//��пװ뾶
	//double		r1 = 2;			//��п��ھ�
	//double		t_back = 2;			//��Բ����
	//double		t_front = 6;			//ǰԲ����

};

//Ҷ��ҶƬ��
class VaneWhell {
private:
	double r = 10;//��챰뾶

public:

};

//������
class Box {
public:
	//�߲����
	double m_length;	//���峤
	double m_width;		//�����
	double m_height;	//�����
	double m_t;			//����ں�
	double m_ht;		//����Ǻ�
	double m_dt;		//�׺�
	double m_len;		//�����淽�׳�
	double m_R;			//����Բ�״�뾶
	double m_r;			//����Բ��С�뾶
	//double m_uplen;		//�����Ϸ��ױ߳�
	double m_dut;		//x��������
	double m_ddt;		//z��������
	double m_dr;		//�ݿװ뾶
	double m_gr;		//ԲͲ�뾶
	double m_grd;		//ԲͲ���
	double m_grt;		//ԲͲ��


private:
	varray<SplineVolume> Boxes;
	varray<SplineSurface> CirTs;
	varray<SplineSurface> Traps;
	//temp
	varray<SplineSurface> patch;
	
public:
	Box(double length, double width, double height, double t, double ht,double dt, double len, 
		double R,double r,double dut,double ddt, double dr,
		double gr, double grd, double grt) :
		m_length(length), m_width(width), m_height(height), m_t(t), m_ht(ht), m_dt(dt),m_len(len),
		m_R(R), m_r(r) ,m_dut(dut),m_ddt(ddt), m_dr(dr),
		m_gr(gr),m_grd(grd),m_grt(grt){
		//����
		trapezoid* tr = new trapezoid(m_len, m_length / 2, (m_height - m_len) / 2, m_length / 4, 0, 0, PI / 2);
		patch.push_back(tr->getSuf());
		tr = new trapezoid(m_len, m_height, (m_length / 2 - m_len) / 2, m_length / 2, m_height / 2, 0, 0);
		patch.push_back(tr->getSuf());
		tr = new trapezoid(m_len, m_length / 2, (m_height - m_len) / 2, m_length / 4, height, 0, 3 * PI / 2);
		patch.push_back(tr->getSuf());
		tr = new trapezoid(m_len, m_height, (m_length / 2 - m_len) / 2, 0, m_height / 2,0, -PI);
		patch.push_back(tr->getSuf());
		
		Cirle_Triangle* cirt = new Cirle_Triangle(m_height / 2, m_length / 4, m_length / 4, 0, m_height / 2, m_R, 3 * m_length / 4, 0 , 0, PI / 2);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_length / 4, m_height / 2, m_height / 2, 0, m_length / 4, m_R, m_length, m_height / 2, 0, PI);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_height / 2, m_length / 4, m_length / 4, 0, m_height / 2, m_R, 3 * m_length / 4, m_height, 0, -PI / 2);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_length / 4, m_height / 2, m_height / 2, 0, m_length / 4, m_R, m_length/2, m_height / 2, 0, 0);
		patch.push_back(cirt->getSuf());

		Cirle_arc* cira = new Cirle_arc(m_r,m_R,2*atan(m_length/2/m_height), 3 * m_length / 4, m_height/2, 0, PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2*atan(2*m_height/m_length), 3*m_length/4, m_height / 2, 0, PI);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2*atan(m_length/2/m_height), 3 * m_length / 4, m_height/2, 0, -PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2*atan(2*m_height/m_length), 3*m_length / 4, m_height / 2, 0, 0);
		patch.push_back(cira->getSuf());
		
		delete cirt;
		cirt = nullptr;
		delete tr;
		tr = nullptr;
		delete cira;
		cira = nullptr;
		Model_Solution M;
		varray<SplineVolume> NV = M.CreatSweepVol(patch, m_t, z_forward);
		for (int i = 0; i < NV.size();i++) {
			Boxes.push_back(NV[i]);
		}
		
		//����
		varray<SplineVolume> NVS = Boxes;
		M.Rolate(NVS, PI, 2);
		M.Trans(NVS, m_length, x_forward);
		M.Trans(NVS, -1*m_width, z_forward);
		for (auto& V : NVS) {
			Boxes.push_back(V);
		}
		//����
	/*	patch.clear();
		tr = new trapezoid(m_uplen, m_length / 2, (m_width - m_uplen) / 2, m_length / 4, 0, 0, PI / 2);
		patch.push_back(tr->getSuf());
		tr = new trapezoid(m_uplen, m_width, (m_length / 2 - m_uplen) / 2, m_length / 2, m_width / 2, 0, 0);
		patch.push_back(tr->getSuf());
		tr = new trapezoid(m_uplen, m_length / 2, (m_width - m_uplen) / 2, m_length / 4, m_width, 0, 3 * PI / 2);
		patch.push_back(tr->getSuf());
		tr = new trapezoid(m_uplen, m_width, (m_length / 2 - m_uplen) / 2, 0, m_width / 2, 0, -PI);
		patch.push_back(tr->getSuf());
		
		RecTangle* rt = new RecTangle(m_length / 2, m_width, 3 * m_length / 4, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		
		NV = M.CreatSweepVol(patch, m_ht, z_forward);

		M.Rolate(NV, -PI / 2, 1);
		M.Trans(NV, m_height, y_forward);

		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		delete tr;
		tr = nullptr;*/

		//����
		patch.clear();
		NV.clear();
		RecTangle* rt = new RecTangle(m_width, m_height, 0, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 2);
		NV = M.CreatSweepVol(patch, -m_t, 1);
		M.Trans(NV, -m_width / 2, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		patch.clear();
		NV.clear();
		rt = new RecTangle(m_width, m_height, 0, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, PI / 2, 2);
		M.Trans(patch, m_length, 1);
		NV = M.CreatSweepVol(patch, m_t, 1);
		M.Trans(NV, -m_width / 2, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		//����
		patch.clear();
		NV.clear();
		rt = new RecTangle(m_length / 2, m_width, m_length / 4, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		rt = new RecTangle(m_length / 2, m_width, 3 * m_length / 4, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		NV = M.CreatSweepVol(patch, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		//ǰ�����Ӵ�
		patch.clear();
		NV.clear();
		rt = new RecTangle(m_length / 2, m_dt, m_length / 4, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		rt = new RecTangle(m_length / 2, m_dt, 3 * m_length / 4, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		NV = M.CreatSweepVol(patch, m_t, 3);
		M.Trans(NV, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_dt, 2);
		M.Trans(NV, m_height, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		M.Trans(NV, -m_t - m_width, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		M.Trans(NV, -m_dt - m_height, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		//�������Ӵ�
		patch.clear();
		NV.clear();
		rt = new RecTangle(m_t, m_dt, 0, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		NV = M.CreatSweepVol(patch, -m_width, 3);
		M.Trans(NV, -m_t / 2, 1);
		M.Trans(NV, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_height + m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_length + m_t, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -m_height- m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		rt = new RecTangle(m_t, m_t, 0, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		M.Trans(patch, -m_t / 2, 1);
		M.Trans(patch, m_t, 3);
		NV = M.CreatSweepVol(patch, m_height, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		NV = M.CreatSweepVol(patch, m_dt, 2);
		M.Trans(NV, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_dt + m_height, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		int blen = Boxes.size();
		NV.clear();
		NV.push_back(Boxes[blen - 3]);
		NV.push_back(Boxes[blen - 2]);
		NV.push_back(Boxes[blen - 1]);
		M.Trans(NV, m_length + m_t, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_width + m_t), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_length + m_t), 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		rt = new RecTangle(m_dut, m_t, 0, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		M.Trans(patch, -m_t-m_dut/2, 1);
		NV = M.CreatSweepVol(patch, m_dt, 2);
		M.Trans(NV, -m_dt, 2);
		M.Trans(NV, m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_width + m_t), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_dut + 2 * m_t + m_length, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_width + m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		rt = new RecTangle(m_dut, m_width, 0, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		NV = M.CreatSweepVol(patch, m_dt, 2);
		M.Trans(NV, -m_dut/2-m_t, 1);
		M.Trans(NV, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_height + m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_dut + 2 * m_t + m_length, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_height + m_dt), 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		cirt = new Cirle_Triangle(m_t / 2, m_dut / 2, m_dut / 2, 0, m_t / 2, m_dr, 0, 0, 0, PI / 2);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_dut / 2, m_t / 2, m_t / 2, 0, m_dut / 2, m_dr, m_dut/2, m_t/2, 0, PI);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_t / 2, m_dut / 2, m_dut / 2, 0, m_t / 2, m_dr, 0, m_t, 0, -PI / 2);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_dut / 2, m_t / 2, m_t / 2, 0, m_dut / 2, m_dr, -m_dut/2, m_t/2, 0, 0);
		patch.push_back(cirt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		NV = M.CreatSweepVol(patch, m_ht, 2);
		M.Trans(NV, m_height, 2);
		M.Trans(NV, -m_dut/2 - m_t, 1);
		M.Trans(NV, m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -m_width - m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_dut + 2 * m_t + m_length, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_t + m_width, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		rt = new RecTangle(m_length / 2, m_ddt, m_length/4, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		NV = M.CreatSweepVol(patch, m_dt, 2);
		M.Trans(NV, -m_dt, 2);
		M.Trans(NV, m_ddt+m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_length / 2, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(2 * m_t + m_width + m_ddt), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -m_length / 2, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		rt = new RecTangle(m_t, m_ddt, -m_t / 2, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		NV = M.CreatSweepVol(patch, m_dt, 2);
		M.Trans(NV, m_ddt + m_t, 3);
		M.Trans(NV, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_ddt + m_width + 2 * m_t), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_length + m_t, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_ddt + m_width + 2 * m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		cirt = new Cirle_Triangle(m_ddt / 2, m_dut / 2, m_dut / 2, 0, m_ddt / 2, m_dr, 0, 0, 0, PI / 2);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_dut / 2, m_ddt / 2, m_ddt / 2, 0, m_dut / 2, m_dr, m_dut / 2, m_ddt / 2, 0, PI);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_ddt / 2, m_dut / 2, m_dut / 2, 0, m_ddt / 2, m_dr, 0, m_ddt, 0, -PI / 2);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_dut / 2, m_ddt / 2, m_ddt / 2, 0, m_dut / 2, m_dr, -m_dut / 2, m_ddt / 2, 0, 0);
		patch.push_back(cirt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		NV = M.CreatSweepVol(patch, m_dt, 2);
		M.Trans(NV, -(m_t + m_dut/2), 1);
		M.Trans(NV, m_t + m_ddt, 3);
		M.Trans(NV, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_ddt + 2 * m_t + m_width), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_length + 2 * m_t + m_dut,1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_ddt + 2 * m_t + m_width, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		patch.clear();
		NV.clear();
		cira = new Cirle_arc(m_r, m_R, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, PI);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, -PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, 0);
		patch.push_back(cira->getSuf());

		NV = M.CreatSweepVol(patch, m_grd, 3);
		M.Trans(NV, m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		cira = new Cirle_arc(m_r, m_R, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, PI);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, -PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, 0);
		patch.push_back(cira->getSuf());

		NV = M.CreatSweepVol(patch, m_grt, 3);
		M.Trans(NV, m_t + m_grd, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		
		patch.clear();
		NV.clear();
		cira = new Cirle_arc(m_R, m_gr, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_R, m_gr, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, PI);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_R, m_gr, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, -PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_R, m_gr, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, 0);
		patch.push_back(cira->getSuf());
		NV = M.CreatSweepVol(patch, m_grt, 3);
		M.Trans(NV, m_t + m_grd, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		blen = Boxes.size();
		for (int i = blen - 1; i >= blen - 12; i--)	NV.push_back(Boxes[i]);
		M.Trans(NV, m_grd + m_grt, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_grd + m_grt, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_grd + m_grt, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		//���ӵĲ���
		//M.Trans(NV, -2 * atan(2 * m_height / m_length), 1);
		M.Trans(NV, 0.5*m_width, 3);
		M.Trans(NV, -0.5*m_length, 1);
		
		
		M.Rolate(NV, PI, 2);
		M.Trans(NV, 0.5*m_length, 1);
		M.Trans(NV, -0.5*m_width, 3);

		//���ӵ�end

		//M.Rolate(NV, PI, 2);

		//M.Trans(NV, -m_grd, 3);
		//M.Trans(NV, m_length, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, (m_grd + m_grt), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, (m_grd + m_grt), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, (m_grd + m_grt), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}



		delete rt;
		rt = nullptr;
		delete cirt;
		cirt = nullptr;

	}
	Box(){}

	varray<SplineVolume> getVols() { return Boxes; };	 
};

class Box2 {
public:
	//�߲����
	double m_length;	//���峤
	double m_width;		//�����
	double m_height;	//�����
	double m_t;			//����ں�
	double m_ht;		//����Ǻ�
	double m_dt;		//�׺�
	double m_len;		//�����淽�׳�
	double m_R;			//����Բ�״�뾶
	double m_r;			//����Բ��С�뾶
	//double m_uplen;		//�����Ϸ��ױ߳�
	double m_dut;		//x��������
	double m_ddt;		//z��������
	double m_dr;		//�ݿװ뾶
	double m_gr;		//ԲͲ�뾶
	double m_grd;		//ԲͲ���
	double m_grt;		//ԲͲ��


private:
	varray<SplineVolume> Boxes;
	varray<SplineSurface> CirTs;
	varray<SplineSurface> Traps;
	//temp
	varray<SplineSurface> patch;

public:
	Box2(double length=8, double width=3.5, double height=4, double t=0.5, double ht=0.5, double dt=0.5, double len=1.5,
		double R=1.5, double r=1, double dut=0.5, double ddt=0.5, double dr=0.15,
		double gr=2.5, double grd=0.6, double grt=0.5) :
		m_length(length), m_width(width), m_height(height), m_t(t), m_ht(ht), m_dt(dt), m_len(len),
		m_R(R), m_r(r), m_dut(dut), m_ddt(ddt), m_dr(dr),
		m_gr(gr), m_grd(grd), m_grt(grt) {
		//����
		trapezoid* tr = new trapezoid(m_len, m_length / 2, (m_height - m_len) / 2, m_length / 4, 0, 0, PI / 2);
		patch.push_back(tr->getSuf());
		tr = new trapezoid(m_len, m_height, (m_length / 2 - m_len) / 2, m_length / 2, m_height / 2, 0, 0);
		patch.push_back(tr->getSuf());
		tr = new trapezoid(m_len, m_length / 2, (m_height - m_len) / 2, m_length / 4, height, 0, 3 * PI / 2);
		patch.push_back(tr->getSuf());
		tr = new trapezoid(m_len, m_height, (m_length / 2 - m_len) / 2, 0, m_height / 2, 0, -PI);
		patch.push_back(tr->getSuf());

		Cirle_Triangle* cirt = new Cirle_Triangle(m_height / 2, m_length / 4, m_length / 4, 0, m_height / 2, m_R, 3 * m_length / 4, 0, 0, PI / 2);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_length / 4, m_height / 2, m_height / 2, 0, m_length / 4, m_R, m_length, m_height / 2, 0, PI);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_height / 2, m_length / 4, m_length / 4, 0, m_height / 2, m_R, 3 * m_length / 4, m_height, 0, -PI / 2);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_length / 4, m_height / 2, m_height / 2, 0, m_length / 4, m_R, m_length / 2, m_height / 2, 0, 0);
		patch.push_back(cirt->getSuf());

		Cirle_arc* cira = new Cirle_arc(m_r, m_R, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, PI);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, -PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, 0);
		patch.push_back(cira->getSuf());

		delete cirt;
		cirt = nullptr;
		delete tr;
		tr = nullptr;
		delete cira;
		cira = nullptr;
		Model_Solution M;
		varray<SplineVolume> NV = M.CreatSweepVol(patch, m_t, z_forward);
		for (int i = 0; i < NV.size(); i++) {
			Boxes.push_back(NV[i]);
		}

		//����
		varray<SplineVolume> NVS = Boxes;
		M.Rolate(NVS, PI, 2);
		M.Trans(NVS, m_length, x_forward);
		M.Trans(NVS, -1 * m_width, z_forward);
		for (auto& V : NVS) {
			Boxes.push_back(V);
		}
		//����
	/*	patch.clear();
		tr = new trapezoid(m_uplen, m_length / 2, (m_width - m_uplen) / 2, m_length / 4, 0, 0, PI / 2);
		patch.push_back(tr->getSuf());
		tr = new trapezoid(m_uplen, m_width, (m_length / 2 - m_uplen) / 2, m_length / 2, m_width / 2, 0, 0);
		patch.push_back(tr->getSuf());
		tr = new trapezoid(m_uplen, m_length / 2, (m_width - m_uplen) / 2, m_length / 4, m_width, 0, 3 * PI / 2);
		patch.push_back(tr->getSuf());
		tr = new trapezoid(m_uplen, m_width, (m_length / 2 - m_uplen) / 2, 0, m_width / 2, 0, -PI);
		patch.push_back(tr->getSuf());

		RecTangle* rt = new RecTangle(m_length / 2, m_width, 3 * m_length / 4, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());

		NV = M.CreatSweepVol(patch, m_ht, z_forward);

		M.Rolate(NV, -PI / 2, 1);
		M.Trans(NV, m_height, y_forward);

		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		delete tr;
		tr = nullptr;*/

		//����
		patch.clear();
		NV.clear();
		RecTangle* rt = new RecTangle(m_width, m_height, 0, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 2);
		NV = M.CreatSweepVol(patch, -m_t, 1);
		M.Trans(NV, -m_width / 2, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		patch.clear();
		NV.clear();
		rt = new RecTangle(m_width, m_height, 0, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, PI / 2, 2);
		M.Trans(patch, m_length, 1);
		NV = M.CreatSweepVol(patch, m_t, 1);
		M.Trans(NV, -m_width / 2, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		//����
		patch.clear();
		NV.clear();
		rt = new RecTangle(m_length / 2, m_width, m_length / 4, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		rt = new RecTangle(m_length / 2, m_width, 3 * m_length / 4, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		NV = M.CreatSweepVol(patch, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		//ǰ�����Ӵ�
		patch.clear();
		NV.clear();
		rt = new RecTangle(m_length / 2, m_dt, m_length / 4, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		rt = new RecTangle(m_length / 2, m_dt, 3 * m_length / 4, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		NV = M.CreatSweepVol(patch, m_t, 3);
		M.Trans(NV, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_dt, 2);
		M.Trans(NV, m_height, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		M.Trans(NV, -m_t - m_width, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		M.Trans(NV, -m_dt - m_height, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		//�������Ӵ�
		patch.clear();
		NV.clear();
		rt = new RecTangle(m_t, m_dt, 0, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		NV = M.CreatSweepVol(patch, -m_width, 3);
		M.Trans(NV, -m_t / 2, 1);
		M.Trans(NV, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_height + m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_length + m_t, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -m_height - m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		rt = new RecTangle(m_t, m_t, 0, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		M.Trans(patch, -m_t / 2, 1);
		M.Trans(patch, m_t, 3);
		NV = M.CreatSweepVol(patch, m_height, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		NV = M.CreatSweepVol(patch, m_dt, 2);
		M.Trans(NV, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_dt + m_height, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		int blen = Boxes.size();
		NV.clear();
		NV.push_back(Boxes[blen - 3]);
		NV.push_back(Boxes[blen - 2]);
		NV.push_back(Boxes[blen - 1]);
		M.Trans(NV, m_length + m_t, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_width + m_t), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_length + m_t), 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		rt = new RecTangle(m_dut, m_t, 0, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		M.Trans(patch, -m_t - m_dut / 2, 1);
		NV = M.CreatSweepVol(patch, m_dt, 2);
		M.Trans(NV, -m_dt, 2);
		M.Trans(NV, m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_width + m_t), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_dut + 2 * m_t + m_length, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_width + m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		rt = new RecTangle(m_dut, m_width, 0, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		NV = M.CreatSweepVol(patch, m_dt, 2);
		M.Trans(NV, -m_dut / 2 - m_t, 1);
		M.Trans(NV, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_height + m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_dut + 2 * m_t + m_length, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_height + m_dt), 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		cirt = new Cirle_Triangle(m_t / 2, m_dut / 2, m_dut / 2, 0, m_t / 2, m_dr, 0, 0, 0, PI / 2);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_dut / 2, m_t / 2, m_t / 2, 0, m_dut / 2, m_dr, m_dut / 2, m_t / 2, 0, PI);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_t / 2, m_dut / 2, m_dut / 2, 0, m_t / 2, m_dr, 0, m_t, 0, -PI / 2);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_dut / 2, m_t / 2, m_t / 2, 0, m_dut / 2, m_dr, -m_dut / 2, m_t / 2, 0, 0);
		patch.push_back(cirt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		NV = M.CreatSweepVol(patch, m_ht, 2);
		M.Trans(NV, m_height, 2);
		M.Trans(NV, -m_dut / 2 - m_t, 1);
		M.Trans(NV, m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -m_width - m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_dut + 2 * m_t + m_length, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_t + m_width, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		rt = new RecTangle(m_length / 2, m_ddt, m_length / 4, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		NV = M.CreatSweepVol(patch, m_dt, 2);
		M.Trans(NV, -m_dt, 2);
		M.Trans(NV, m_ddt + m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_length / 2, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(2 * m_t + m_width + m_ddt), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -m_length / 2, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		rt = new RecTangle(m_t, m_ddt, -m_t / 2, 0, 0, PI / 2);
		patch.push_back(rt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		NV = M.CreatSweepVol(patch, m_dt, 2);
		M.Trans(NV, m_ddt + m_t, 3);
		M.Trans(NV, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_ddt + m_width + 2 * m_t), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_length + m_t, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_ddt + m_width + 2 * m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		cirt = new Cirle_Triangle(m_ddt / 2, m_dut / 2, m_dut / 2, 0, m_ddt / 2, m_dr, 0, 0, 0, PI / 2);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_dut / 2, m_ddt / 2, m_ddt / 2, 0, m_dut / 2, m_dr, m_dut / 2, m_ddt / 2, 0, PI);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_ddt / 2, m_dut / 2, m_dut / 2, 0, m_ddt / 2, m_dr, 0, m_ddt, 0, -PI / 2);
		patch.push_back(cirt->getSuf());
		cirt = new Cirle_Triangle(m_dut / 2, m_ddt / 2, m_ddt / 2, 0, m_dut / 2, m_dr, -m_dut / 2, m_ddt / 2, 0, 0);
		patch.push_back(cirt->getSuf());
		M.Rolate(patch, -PI / 2, 1);
		NV = M.CreatSweepVol(patch, m_dt, 2);
		M.Trans(NV, -(m_t + m_dut / 2), 1);
		M.Trans(NV, m_t + m_ddt, 3);
		M.Trans(NV, -m_dt, 2);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_ddt + 2 * m_t + m_width), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_length + 2 * m_t + m_dut, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_ddt + 2 * m_t + m_width, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		patch.clear();
		NV.clear();
		cira = new Cirle_arc(m_r, m_R, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, PI);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, -PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, 0);
		patch.push_back(cira->getSuf());

		NV = M.CreatSweepVol(patch, m_grd, 3);
		M.Trans(NV, m_t, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		cira = new Cirle_arc(m_r, m_R, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, PI);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, -PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_r, m_R, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, 0);
		patch.push_back(cira->getSuf());

		NV = M.CreatSweepVol(patch, m_grt, 3);
		M.Trans(NV, m_t + m_grd, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		cira = new Cirle_arc(m_R, m_gr, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_R, m_gr, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, PI);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_R, m_gr, 2 * atan(m_length / 2 / m_height), 3 * m_length / 4, m_height / 2, 0, -PI / 2);
		patch.push_back(cira->getSuf());
		cira = new Cirle_arc(m_R, m_gr, 2 * atan(2 * m_height / m_length), 3 * m_length / 4, m_height / 2, 0, 0);
		patch.push_back(cira->getSuf());
		NV = M.CreatSweepVol(patch, m_grt, 3);
		M.Trans(NV, m_t + m_grd, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		patch.clear();
		NV.clear();
		blen = Boxes.size();
		for (int i = blen - 1; i >= blen - 12; i--)	NV.push_back(Boxes[i]);
		M.Trans(NV, m_grd + m_grt, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_grd + m_grt, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, m_grd + m_grt, 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}

		M.Rolate(NV, PI, 2);

		M.Trans(NV, -m_grd, 3);
		M.Trans(NV, m_length, 1);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_grd + m_grt), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_grd + m_grt), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}
		M.Trans(NV, -(m_grd + m_grt), 3);
		for (auto &vol : NV) {
			Boxes.push_back(vol);
		}



		delete rt;
		rt = nullptr;
		delete cirt;
		cirt = nullptr;

	}
	Box2() {}

	varray<SplineVolume> getVols() { return Boxes; };
};


class Reducer {
private:
	//�߲����
	double	 m_gear_r1;					//����ְ뾶
	double	 m_gear_r2;					//С���ְ뾶
	double	 m_shaft_r1;				//�������뾶(��Բ)
	double	 m_shaft_r2;				//С������뾶
	double	 m_shaft_r3;				//С������뾶
	double	 m_shaft_l1;				//��������С������ĳ���L1
	double	m_shaft_l2;					//С������L2
	double	m_shaft_h;					//������ĸ߶�h
	double	m_shaft_t;					//����Ͳ��t
	double	m_box_L;					//�����L
	double	m_box_t;					//�����t
private:
	varray<SplineVolume>			front;
	varray<SplineVolume>			reducer;
	varray<SplineVolume>			tempVols;
	varray<SplineSurface>		outer;
	varray<SplineSurface>		mid;
	varray<SplineSurface>		ribs;
	varray<SplineSurface>		outer_boundary;
	varray<double>				m_knots;
	Spline					path;
	SplineSurface				gear_cir;
	varray<SplineSurface>		NS;
	double						m_degree = 2;
	double						m_ctrlNum = 3;
	Model_Solution				M;
private:
	//�м��
	
	double m_box_w;					//���峤w
	double box_outer_t = 0.4*m_shaft_r2;	//������Ե��
	double box_mid_t = m_shaft_r2;		//�����в�����̨��
	double box_rib_t = 0.28;					//����
	double d1 = m_shaft_l1;		//��Բ��СԲ�����ľ�
	double d2 = m_shaft_l2;		//СԲ��СԲ�����ľ�			
	//�ײ�

	Vec4 p1 = { 0,box_outer_t / 2 + m_shaft_h,0 };
	Vec4 p0, p2, p3,p5;
	Vec4 p4;
	Vec4 b1;
	Vec4 b0,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12;
	Vec4 c0, c1, c2, c3, c4, c5, c6, c7;
public:
	//Reducer(){}
	Reducer(double gr1 = 5, double gr2 = 3, double sr1 = 2, double sr2 = 1.5, double sr3 = 1.5,
		double sl1 = 5, double sl2 = 5, double sh = 10, double bl = 6, double bt = 0.5, double st = 0.5):
		m_gear_r1(gr1), m_gear_r2(gr2), m_shaft_l1(sl1), m_shaft_l2(sl2),
		m_shaft_r1(sr1), m_shaft_r2(sr2), m_shaft_r3(sr3),
		m_shaft_h(sh), m_box_L(bl), m_box_t(bt), m_shaft_t(st)
	{
		m_knots.push_back(0);
		m_knots.push_back(0);
		m_knots.push_back(0);
		m_knots.push_back(1);
		m_knots.push_back(1);
		m_knots.push_back(1);
		double x1 = sqrt((m_gear_r1)*(m_gear_r1)-(box_outer_t / 2)*(box_outer_t / 2));
		double y1 = m_shaft_h;
		double x2 = x1 + m_shaft_l1;
		double y2 = m_shaft_h;
		double x3 = x1 + m_shaft_l1 + m_shaft_l2;
		double y3 = y1;
		double ceta = asin((m_gear_r1 - m_gear_r2) / (d1 + d2));
		double alph = PI / 2 - ceta;
		m_box_w = x3 + sqrt(m_gear_r2 * m_gear_r2 - box_outer_t * box_outer_t / 4);
		b1 = { x1 - sqrt(m_shaft_r1*m_shaft_r1 - (box_outer_t*box_outer_t / 4)),p1.y,0 };
		p2 = { m_gear_r1*cos(alph) + x1, m_gear_r1*sin(alph) + y1,0 };
		p4 = { m_box_w,box_outer_t / 2 + m_shaft_h,0 };
		p0 = Model_Solution::NurbsCirle_mid_pt(p1, p2, x1, y1, m_gear_r1);
		b2 = { m_shaft_r1*cos(alph) + x1, m_shaft_r1*sin(alph) + y1,0 };
		b0 = Model_Solution::NurbsCirle_mid_pt(b1, b2, x1, y1, m_shaft_r1);
		varray<Vec4> pt1;
		pt1.push_back(p1);
		pt1.push_back(p0);
		pt1.push_back(p2);
		varray<Vec4> pt2;
		pt2.push_back(b1);
		pt2.push_back(b0);
		pt2.push_back(b2);
		gear_cir = gear_create(pt1, pt2);
		NS.push_back(gear_cir);
		tempVols.clear();
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		p3 = { m_gear_r2*cos(alph) + x3, m_gear_r2*sin(alph) + y3 ,0 };
		p5 = Model_Solution::NurbsCirle_mid_pt(p3, p4, x3, y3, m_gear_r2);
		b4 = { m_shaft_r3*cos(alph) + x3 ,m_shaft_r3*sin(alph) + y3,0 };
		b5 = { x3 + sqrt(m_shaft_r3*m_shaft_r3 - (box_outer_t / 2)*(box_outer_t / 2)),y3 + box_outer_t / 2,0 };
		b3 = Model_Solution::NurbsCirle_mid_pt(b4, b5, x3, y3, m_shaft_r3);
		pt1.clear(), pt2.clear();
		pt1.push_back(b5);
		pt1.push_back(b3);
		pt1.push_back(b4);
		pt2.push_back(p4);
		pt2.push_back(p5);
		pt2.push_back(p3);
		NS.clear();
		gear_cir = gear_create(pt1, pt2);
		NS.push_back(gear_cir);
		tempVols.clear();
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		b6 = { x2 - sqrt(m_shaft_r2*m_shaft_r2 - (box_mid_t / 2)*(box_mid_t / 2)),y3 + box_mid_t / 2,0 };
		b7 = { 2 * x2 - b6.x, b6.y,0 };
		b8 = Model_Solution::NurbsCirle_mid_pt(b6, b7, x2, y2, m_shaft_r2);
		pt1.clear(), pt2.clear();
		pt1.push_back(b6);
		pt1.push_back(b8);
		pt1.push_back(b7);
		pt2.push_back(p2);
		pt2.push_back(p3);
		NS.clear();
		tempVols.clear();
		SplineSurface Tri_cir = Tri_circle_create(pt1, pt2);
		NS.push_back(Tri_cir);
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		b9 = { x1 + sqrt(m_shaft_r1*m_shaft_r1 - (box_mid_t / 2)*(box_mid_t / 2)),y1 + box_mid_t / 2,0 };
		b10 = Model_Solution::NurbsCirle_mid_pt(b2, b9, x1, y1, m_shaft_r1);
		pt1.clear(), pt2.clear();
		pt1.push_back(b2);
		pt1.push_back(b10);
		pt1.push_back(b9);
		pt2.push_back(p2);
		pt2.push_back(b6);
		Tri_cir = Tri_circle_create(pt1, pt2);
		NS.clear();
		tempVols.clear();
		NS.push_back(Tri_cir);
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		b11 = { x3 - sqrt(m_shaft_r3*m_shaft_r3 - box_mid_t * box_mid_t / 4),y3 + box_mid_t / 2,0 };
		b12 = Model_Solution::NurbsCirle_mid_pt(b11, b4, x3, y3, m_shaft_r3);
		pt1.clear(), pt2.clear();
		NS.clear();
		tempVols.clear();
		pt1.push_back(b11);
		pt1.push_back(b12);
		pt1.push_back(b4);
		pt2.push_back(b7);
		pt2.push_back(p3);
		Tri_cir = Tri_circle_create(pt1, pt2);
		NS.push_back(Tri_cir);
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		c0 = { 0,m_shaft_h,0 };
		c1 = { x1 - m_shaft_r1, y1, 0 };
		Vec4 d1 = Model_Solution::NurbsCirle_mid_pt(c1, b1, x1, y1, m_shaft_r1);
		pt1.clear(), pt2.clear();
		pt1.push_back(c1);
		pt1.push_back(d1);
		pt1.push_back(b1);
		pt2.push_back(c0);
		pt2.push_back(p1);
		NS.clear();
		tempVols.clear();
		Tri_cir = Tri_circle_create(pt1, pt2);
		NS.push_back(Tri_cir);
		outer.push_back(Tri_cir);
		Vec4 p1t = { p1.x, 2 * m_shaft_h - p1.y,0 };
		Vec4 b1t = { b1.x, 2 * m_shaft_h - b1.y,0 };
		for (auto& x : Tri_cir.m_CtrlPts) {
			x.y = 2 * m_shaft_h - x.y;
		}
		NS.push_back(Tri_cir);
		outer.push_back(Tri_cir);
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		c2 = { x1 + m_shaft_r1 , y2,0 };
		c3 = { x2 - m_shaft_r2, y2, 0 };
		Vec4 d2 = Model_Solution::NurbsCirle_mid_pt(c2, b9, x1, y1, m_shaft_r1);
		Vec4 d3 = Model_Solution::NurbsCirle_mid_pt(c3, b6, x2, y2, m_shaft_r2);
		pt1.clear(), pt2.clear();
		pt1.push_back(c2);
		pt1.push_back(d2);
		pt1.push_back(b9);
		pt2.push_back(c3);
		pt2.push_back(d3);
		pt2.push_back(b6);
		SplineSurface Rect_Cir = Rect_Cir_2create(pt1, pt2);
		NS.clear();
		tempVols.clear();
		NS.push_back(Rect_Cir);
		mid.push_back(Rect_Cir);
		for (auto& x : Rect_Cir.m_CtrlPts) {
			x.y = 2 * m_shaft_h - x.y;
		}
		NS.push_back(Rect_Cir);
		mid.push_back(Rect_Cir);
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		c4 = { x2 + m_shaft_r2, y2, 0 };
		c5 = { x3 - m_shaft_r3 , y3, 0 };
		Vec4 d4 = Model_Solution::NurbsCirle_mid_pt(c4, b7, x2, y2, m_shaft_r2);
		Vec4 d5 = Model_Solution::NurbsCirle_mid_pt(c5, b11, x3, y3, m_shaft_r3);
		NS.clear();
		tempVols.clear();
		pt1.clear(), pt2.clear();
		pt1.push_back(c4);
		pt1.push_back(d4);
		pt1.push_back(b7);
		pt2.push_back(c5);
		pt2.push_back(d5);
		pt2.push_back(b11);
		Rect_Cir = Rect_Cir_2create(pt1, pt2);
		NS.push_back(Rect_Cir);
		mid.push_back(Rect_Cir);
		for (auto& x : Rect_Cir.m_CtrlPts) {
			x.y = 2 * m_shaft_h - x.y;
		}
		NS.push_back(Rect_Cir);
		mid.push_back(Rect_Cir);

		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		c6 = { x3 + m_shaft_r3, y3, 0 };
		c7 = { m_box_w, y3, 0 };
		Vec4 d6 = Model_Solution::NurbsCirle_mid_pt(c6, b5, x3, y3, m_shaft_r3);
		pt1.clear(), pt2.clear();
		pt1.push_back(c6);
		pt1.push_back(d6);
		pt1.push_back(b5);
		pt2.push_back(c7);
		pt2.push_back(p4);

		Rect_Cir = Tri_circle_create(pt1, pt2);
		NS.push_back(Rect_Cir);
		outer.push_back(Rect_Cir);
		for (auto& x : Rect_Cir.m_CtrlPts) {
			x.y = 2 * m_shaft_h - x.y;
		}
		NS.push_back(Rect_Cir);
		outer.push_back(Rect_Cir);
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		Vec4 e1 = { x1 - box_rib_t / 2 , y1 - sqrt(m_shaft_r1*m_shaft_r1 - (box_rib_t / 2)*(box_rib_t / 2)),0 };
		Vec4 e5 = Model_Solution::NurbsCirle_mid_pt(b1t, e1, x1, y1, m_shaft_r1);
		Vec4 e3 = { b1t.x,0,0 };
		Vec4 e4 = { e1.x,0,0 };

		pt1.clear(), pt2.clear();
		pt1.push_back(b1t);
		pt1.push_back(e5);
		pt1.push_back(e1);
		pt2.push_back(e3);
		pt2.push_back(e4);
		NS.clear();
		tempVols.clear();
		Rect_Cir = Tri_circle_create(pt1, pt2);
		NS.push_back(Rect_Cir);
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		Vec4 e2 = { 2 * x1 - e1.x, e1.y , 0 };
		Vec4 b9t = { b9.x, 2 * y1 - b9.y,0 };
		Vec4 e6 = Model_Solution::NurbsCirle_mid_pt(e2, b9t, x1, y1, m_shaft_r1);
		Vec4 e7 = { 2 * x1 - e4.x , e4.y ,0 };
		Vec4 e8 = { 2 * x1 - e3.x, e3.y,0 };
		pt1.clear(), pt2.clear();
		pt1.push_back(e2);
		pt1.push_back(e6);
		pt1.push_back(b9t);
		pt2.push_back(e7);
		pt2.push_back(e8);
		NS.clear();
		tempVols.clear();
		Rect_Cir = Tri_circle_create(pt1, pt2);
		NS.push_back(Rect_Cir);
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		Vec4 e9 = Model_Solution::NurbsCirle_mid_pt(e1, e2, x1, y1, m_shaft_r1);
		pt1.clear(), pt2.clear();
		pt1.push_back(e1);
		pt1.push_back(e9);
		pt1.push_back(e2);
		pt2.push_back(e4);
		pt2.push_back(e7);
		Rect_Cir = Tri_circle_create(pt1, pt2);
		NS.clear();
		NS.push_back(Rect_Cir);
		ribs.push_back(Rect_Cir);
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		Vec4 b6t = { b6.x, 2 * y2 - b6.y, 0 };
		Vec4 b7t = { b7.x , 2 * y2 - b7.y,0 };
		Vec4 b11t = { b11.x, 2 * y2 - b11.y,0 };
		Vec4 b5t = { b5.x, 2 * y2 - b5.y ,0 };
		Vec4 p4t = { p4.x , 2 * y2 - p4.y,0 };
		Vec4 e10 = { x2 - box_rib_t / 2, y2 - sqrt(m_shaft_r2* m_shaft_r2 - box_rib_t * box_rib_t / 4),0 };
		Vec4 e11 = { x2 + box_rib_t / 2, y2 - sqrt(m_shaft_r2* m_shaft_r2 - box_rib_t * box_rib_t / 4),0 };
		Vec4 e12 = { b6t.x, 0,0 };
		Vec4 e13 = { e10.x, 0,0 };
		Vec4 e14 = { e11.x, 0,0 };
		Vec4 e15 = { b7t.x, 0,0 };
		Vec4 e16 = Model_Solution::NurbsCirle_mid_pt(b6t, e10, x2, y2, m_shaft_r2);
		Vec4 e17 = Model_Solution::NurbsCirle_mid_pt(e11, b7t, x2, y2, m_shaft_r2);
		Vec4 e18 = Model_Solution::NurbsCirle_mid_pt(e10, e11, x2, y2, m_shaft_r2);
		pt1.clear(), pt2.clear();
		pt1.push_back(b6t);
		pt1.push_back(e16);
		pt1.push_back(e10);
		pt2.push_back(e12);
		pt2.push_back(e13);
		Rect_Cir = Tri_circle_create(pt1, pt2);
		NS.clear();
		NS.push_back(Rect_Cir);
		pt1.clear(), pt2.clear();
		pt1.push_back(e11);
		pt1.push_back(e17);
		pt1.push_back(b7t);
		pt2.push_back(e14);
		pt2.push_back(e15);
		Rect_Cir = Tri_circle_create(pt1, pt2);
		NS.push_back(Rect_Cir);
		pt1.clear(), pt2.clear();
		pt1.push_back(e10);
		pt1.push_back(e18);
		pt1.push_back(e11);
		pt2.push_back(e13);
		pt2.push_back(e14);
		Rect_Cir = Tri_circle_create(pt1, pt2);
		NS.push_back(Rect_Cir);
		ribs.push_back(Rect_Cir);
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		Vec4 e19 = { x3 - box_rib_t / 2 , y3 - sqrt(m_shaft_r3*m_shaft_r3 - box_rib_t * box_rib_t / 4) , 0 };
		Vec4 e20 = { 2 * x3 - e19.x, e19.y, 0 };
		Vec4 e21 = Model_Solution::NurbsCirle_mid_pt(e19, e20, x3, y3, m_shaft_r3);
		Vec4 e22 = Model_Solution::NurbsCirle_mid_pt(b11t, e19, x3, y3, m_shaft_r3);
		Vec4 e23 = Model_Solution::NurbsCirle_mid_pt(e21, b5t, x3, y3, m_shaft_r3);
		Vec4 e24 = { b11t.x,0,0 };
		Vec4 e25 = { e19.x,0,0 };
		Vec4 e26 = { e20.x,0,0 };
		Vec4 e27 = { b5t.x,0,0 };
		pt1.clear(), pt2.clear();
		NS.clear();
		pt1.push_back(b11t);
		pt1.push_back(e22);
		pt1.push_back(e19);
		pt2.push_back(e24);
		pt2.push_back(e25);
		Rect_Cir = Tri_circle_create(pt1, pt2);
		NS.push_back(Rect_Cir);
		pt1.clear(), pt2.clear();
		pt1.push_back(e19);
		pt1.push_back(e21);
		pt1.push_back(e20);
		pt2.push_back(e25);
		pt2.push_back(e26);
		Rect_Cir = Tri_circle_create(pt1, pt2);
		NS.push_back(Rect_Cir);
		ribs.push_back(Rect_Cir);
		pt1.clear(), pt2.clear();
		pt1.push_back(e20);
		pt1.push_back(e23);
		pt1.push_back(b5t);
		pt2.push_back(e26);
		pt2.push_back(e27);
		Rect_Cir = Tri_circle_create(pt1, pt2);
		NS.push_back(Rect_Cir);
		tempVols = M.CreatSweepVol(NS, 1, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		//Բ��
		NS.clear();
		pt1.clear();
		pt2.clear();
		pt1.push_back(b1t);
		pt1.push_back(e5);
		pt1.push_back(e1);
		SplineSurface cir_arc = scaling_cir_creat(pt1, x1, y1, m_shaft_r1, m_shaft_t, 1);
		NS.push_back(cir_arc);
		pt1.clear();
		pt1.push_back(b1);
		pt1.push_back(b0);
		pt1.push_back(b2);
		cir_arc = scaling_cir_creat(pt1, x1, y1, m_shaft_r1, m_shaft_t, 1);
		NS.push_back(cir_arc);
		pt1.clear();
		pt1.push_back(b1t);
		pt1.push_back(c1);
		cir_arc = scaling_cir_creat(pt1, x1, y1, m_shaft_r1, m_shaft_t, 1);
		NS.push_back(cir_arc);
		pt1.clear();
		pt1.push_back(c1);
		pt1.push_back(b1);
		cir_arc = scaling_cir_creat(pt1, x1, y1, m_shaft_r1, m_shaft_t, 1);
		NS.push_back(cir_arc);
		pt1.clear();
		pt1.push_back(e1);
		pt1.push_back(e2);
		cir_arc = scaling_cir_creat(pt1, x1, y1, m_shaft_r1, m_shaft_t, 1);
		NS.push_back(cir_arc);
		pt1.clear();
		pt1.push_back(e2);
		pt1.push_back(b9t);
		cir_arc = scaling_cir_creat(pt1, x1, y1, m_shaft_r1, m_shaft_t, 1);
		NS.push_back(cir_arc);
		pt1.clear();
		pt1.push_back(b9t);
		pt1.push_back(c2);
		cir_arc = scaling_cir_creat(pt1, x1, y1, m_shaft_r1, m_shaft_t, 1);
		NS.push_back(cir_arc);
		pt1.clear();
		pt1.push_back(c2);
		pt1.push_back(b9);
		cir_arc = scaling_cir_creat(pt1, x1, y1, m_shaft_r1, m_shaft_t, 1);
		NS.push_back(cir_arc);
		pt1.clear();
		pt1.push_back(b9);
		pt1.push_back(b2);
		cir_arc = scaling_cir_creat(pt1, x1, y1, m_shaft_r1, m_shaft_t, 1);
		NS.push_back(cir_arc);
		pt1.clear();
		pt1.push_back(b6t);
		pt1.push_back(c3);
		cir_arc = scaling_cir_creat(pt1, x2, y2, m_shaft_r2, m_shaft_t, 1);
		NS.push_back(cir_arc);
		pt1.clear();
		pt1.push_back(c3);
		pt1.push_back(b6);
		cir_arc = scaling_cir_creat(pt1, x2, y2, m_shaft_r2, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(b6);
		pt1.push_back(b7);
		cir_arc = scaling_cir_creat(pt1, x2, y2, m_shaft_r2, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(b7t);
		pt1.push_back(c4);
		cir_arc = scaling_cir_creat(pt1, x2, y2, m_shaft_r2, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(c4);
		pt1.push_back(b7);
		cir_arc = scaling_cir_creat(pt1, x2, y2, m_shaft_r2, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(e10);
		pt1.push_back(b6t);
		cir_arc = scaling_cir_creat(pt1, x2, y2, m_shaft_r2, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(e10);
		pt1.push_back(e11);
		cir_arc = scaling_cir_creat(pt1, x2, y2, m_shaft_r2, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(e11);
		pt1.push_back(b7t);
		cir_arc = scaling_cir_creat(pt1, x2, y2, m_shaft_r2, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(e19);
		pt1.push_back(e20);
		cir_arc = scaling_cir_creat(pt1, x3, y3, m_shaft_r3, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(e19);
		pt1.push_back(b11t);
		cir_arc = scaling_cir_creat(pt1, x3, y3, m_shaft_r3, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(e20);
		pt1.push_back(b5t);
		cir_arc = scaling_cir_creat(pt1, x3, y3, m_shaft_r3, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(b5t);
		pt1.push_back(c6);
		cir_arc = scaling_cir_creat(pt1, x3, y3, m_shaft_r3, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(c6);
		pt1.push_back(b5);
		cir_arc = scaling_cir_creat(pt1, x3, y3, m_shaft_r3, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(b5);
		pt1.push_back(b4);
		cir_arc = scaling_cir_creat(pt1, x3, y3, m_shaft_r3, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(b11);
		pt1.push_back(b4);
		cir_arc = scaling_cir_creat(pt1, x3, y3, m_shaft_r3, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(b11t);
		pt1.push_back(c5);
		cir_arc = scaling_cir_creat(pt1, x3, y3, m_shaft_r3, m_shaft_t, 1);
		NS.push_back(cir_arc);

		pt1.clear();
		pt1.push_back(c5);
		pt1.push_back(b11);
		cir_arc = scaling_cir_creat(pt1, x3, y3, m_shaft_r3, m_shaft_t, 1);
		NS.push_back(cir_arc);

		for (auto&x : NS) { ribs.push_back(x); }

		tempVols = M.CreatSweepVol(NS, m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		tempVols = M.CreatSweepVol(ribs, m_box_t, z_forward);
		M.Trans(tempVols, m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		Vec4 T0 = { 0,0,0 };
		pt1.clear(), pt2.clear();
		pt1.push_back(T0);
		pt1.push_back(e3);
		pt2.push_back(p1t);
		pt2.push_back(b1t);
		SplineSurface Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);

		pt1.clear(), pt2.clear();
		pt1.push_back(e8);
		pt1.push_back(e12);
		pt2.push_back(b9t);
		pt2.push_back(b6t);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);

		pt1.clear(), pt2.clear();
		pt1.push_back(e15);
		pt1.push_back(e24);
		pt2.push_back(b7t);
		pt2.push_back(b11t);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);

		Vec4 T1 = { m_box_w,0,0 };
		pt1.clear(), pt2.clear();
		pt1.push_back(e27);
		pt1.push_back(T1);
		pt2.push_back(b5t);
		pt2.push_back(p4t);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		tempVols = M.CreatSweepVol(NS, m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		tempVols = M.CreatSweepVol(outer, m_box_t, z_forward);
		M.Trans(tempVols, m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		tempVols = M.CreatSweepVol(mid, m_box_t, z_forward);
		M.Trans(tempVols, m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		Vec4 bT0 = T0.Trans(-m_box_t, 0, 0);
		Vec4 bp1t = p1t.Trans(-m_box_t, 0, 0);
		pt1.clear(), pt2.clear();
		pt1.push_back(bT0);
		pt1.push_back(T0);
		pt2.push_back(bp1t);
		pt2.push_back(p1t);
		Rect = Rect_Create(pt1, pt2);
		NS.clear();
		NS.push_back(Rect);
		outer_boundary.push_back(Rect);
		Vec4 bc0 = c0.Trans(-m_box_t, 0, 0);
		Vec4 bp1 = p1.Trans(-m_box_t, 0, 0);
		pt1.clear(), pt2.clear();
		pt1.push_back(bp1t);
		pt1.push_back(p1t);
		pt2.push_back(bc0);
		pt2.push_back(c0);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		outer_boundary.push_back(Rect);
		pt1.clear(), pt2.clear();
		pt1.push_back(bc0);
		pt1.push_back(c0);
		pt2.push_back(bp1);
		pt2.push_back(p1);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		outer_boundary.push_back(Rect);
		varray<SplineSurface> nss = outer_boundary;
		M.Trans(nss, m_box_w + m_box_t, x_forward);
		for (auto& n : nss) {
			outer_boundary.push_back(n);
		}
		tempVols = M.CreatSweepVol(NS, m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, m_box_w + m_box_t, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		NS.clear();
		double k = (m_gear_r1 + m_box_t) / m_gear_r1;
		Vec4 bp2 = { x1 + k * (p2.x - x1) , y1 + k * (p2.y - y1) , p2.z , p2.w };
		Vec4 bp0 = Model_Solution::NurbsCirle_mid_pt(bp1, bp2, x1, y1, m_shaft_r1);
		pt1.clear(), pt2.clear();
		pt1.push_back(bp1);
		pt1.push_back(bp0);
		pt1.push_back(bp2);
		pt2.push_back(p1);
		pt2.push_back(p0);
		pt2.push_back(p2);
		cir_arc = Rect_Cir_2create(pt1, pt2);
		NS.push_back(cir_arc);
		outer_boundary.push_back(cir_arc);
		tempVols = M.CreatSweepVol(NS, m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		k = (m_gear_r2 + m_box_t) / m_gear_r2;
		Vec4 bp3 = { x3 + k * (p3.x - x3) , y3 + k * (p3.y - y3) , p3.z , p3.w };
		Vec4 bp4 = p4.Trans(m_box_t, 0, 0);
		Vec4 bp5 = Model_Solution::NurbsCirle_mid_pt(bp3, bp4, x3, y3, m_box_t);
		pt1.clear(), pt2.clear();
		NS.clear();
		pt1.push_back(p4);
		pt1.push_back(p5);
		pt1.push_back(p3);
		pt2.push_back(bp4);
		pt2.push_back(bp5);
		pt2.push_back(bp3);
		cir_arc = Rect_Cir_2create(pt1, pt2);
		NS.push_back(cir_arc);
		outer_boundary.push_back(cir_arc);
		tempVols = M.CreatSweepVol(NS, m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		pt1.clear(), pt2.clear();
		pt1.push_back(p2);
		pt1.push_back(p3);
		pt2.push_back(bp2);
		pt2.push_back(bp3);
		Rect = Rect_Create(pt1, pt2);
		NS.clear();
		NS.push_back(Rect);
		tempVols = M.CreatSweepVol(NS, m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		front = reducer;
		for (auto& v : front) {
			for (auto& pt : v.m_CtrlPts) {
				pt.z = -m_box_L - pt.z;
			}
		}
		for (auto& v : front) reducer.push_back(v);

		tempVols = M.CreatSweepVol(outer_boundary, -m_box_L, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		pt1.clear(), pt2.clear();
		pt1.push_back(bc0);
		pt1.push_back(c0);
		pt2.push_back(bp1);
		pt2.push_back(p1);
		Rect = Rect_Create(pt1, pt2);
		NS.clear();
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		pt1.push_back(bp1t);
		pt1.push_back(p1t);
		pt2.push_back(bc0);
		pt2.push_back(c0);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		tempVols = M.CreatSweepVol(NS, m_box_t, z_forward);
		M.Trans(tempVols, m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, m_box_w + m_box_t, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, -m_box_L - 3 * m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, -(m_box_w + m_box_t), x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		double bL = sqrt((p2.x - p3.x)*(p2.x - p3.x) + (p2.y - p3.y)*(p2.y - p3.y));
		double sL = 0.8*bL;
		double bW = m_box_L;
		double sW = 0.8*bW;
		NS = Trape4_Create(bL, sL, bW, sW);
		double dalph = atan((p2.y - p3.y) / (p3.x - p2.x));
		M.Rolate(NS, -PI / 2, x_forward);
		tempVols = M.CreatSweepVol(NS, m_box_t, y_forward);
		M.Rolate(tempVols, -dalph, z_forward);
		M.Trans(tempVols, (p2.x + p3.x) / 2, x_forward);
		M.Trans(tempVols, (p2.y + p3.y) / 2, y_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		Vec4 btT0 = T0.Trans(0, -m_box_t, 0);
		Vec4 be3 = e3.Trans(0, -m_box_t, 0);
		pt1.clear(), pt2.clear();
		pt1.push_back(btT0);
		pt1.push_back(be3);
		pt2.push_back(T0);
		pt2.push_back(e3);
		NS.clear();
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		Vec4 be4 = e4.Trans(0, -m_box_t, 0);
		pt1.clear(), pt2.clear();
		pt1.push_back(be3);
		pt1.push_back(be4);
		pt2.push_back(e3);
		pt2.push_back(e4);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		Vec4 be7 = e7.Trans(0, -m_box_t, 0);
		pt1.push_back(be4);
		pt1.push_back(be7);
		pt2.push_back(e4);
		pt2.push_back(e7);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		Vec4 be8 = e8.Trans(0, -m_box_t, 0);
		pt1.clear(), pt2.clear();
		pt1.push_back(be7);
		pt1.push_back(be8);
		pt2.push_back(e7);
		pt2.push_back(e8);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		Vec4 be12 = e12.Trans(0, -m_box_t, 0);
		pt1.push_back(be8);
		pt1.push_back(be12);
		pt2.push_back(e8);
		pt2.push_back(e12);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		Vec4 be13 = e13.Trans(0, -m_box_t, 0);
		pt1.push_back(be12);
		pt1.push_back(be13);
		pt2.push_back(e12);
		pt2.push_back(e13);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		Vec4 be14 = e14.Trans(0, -m_box_t, 0);
		pt1.push_back(be13);
		pt1.push_back(be14);
		pt2.push_back(e13);
		pt2.push_back(e14);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		Vec4 be15 = e15.Trans(0, -m_box_t, 0);
		pt1.push_back(be14);
		pt1.push_back(be15);
		pt2.push_back(e14);
		pt2.push_back(e15);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		Vec4 be24 = e24.Trans(0, -m_box_t, 0);
		pt1.clear(), pt2.clear();
		pt1.push_back(be15);
		pt1.push_back(be24);
		pt2.push_back(e15);
		pt2.push_back(e24);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		Vec4 be25 = e25.Trans(0, -m_box_t, 0);
		pt1.push_back(be24);
		pt1.push_back(be25);
		pt2.push_back(e24);
		pt2.push_back(e25);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		Vec4 be26 = e26.Trans(0, -m_box_t, 0);
		pt1.push_back(be25);
		pt1.push_back(be26);
		pt2.push_back(e25);
		pt2.push_back(e26);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		Vec4 be27 = e27.Trans(0, -m_box_t, 0);
		pt1.push_back(be26);
		pt1.push_back(be27);
		pt2.push_back(e26);
		pt2.push_back(e27);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		Vec4 e28 = { p4t.x,0,0 };
		Vec4 be28 = e28.Trans(0, -m_box_t, 0);
		pt1.push_back(be27);
		pt1.push_back(be28);
		pt2.push_back(e27);
		pt2.push_back(e28);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		Vec4 ebT0 = bT0.Trans(0, -m_box_t, 0);
		pt1.push_back(ebT0);
		pt1.push_back(btT0);
		pt2.push_back(bT0);
		pt2.push_back(T0);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		M.Trans(Rect, m_box_t + m_box_w, x_forward);
		NS.push_back(Rect);
		tempVols = M.CreatSweepVol(NS, m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		varray<SplineVolume> v1 = tempVols;
		for (auto& v : v1) {
			for (auto& pt : v.m_CtrlPts) {
				pt.z = -m_box_L - pt.z;
			}
		}
		for (auto&x : v1) {
			reducer.push_back(x);
		}
		M.Trans(v1, -m_box_t, z_forward);
		for (auto&x : v1) {
			reducer.push_back(x);
		}

		M.Trans(tempVols, m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		NS.clear();
		pt1.clear(), pt2.clear();
		pt1.push_back(ebT0);
		pt1.push_back(btT0);
		pt2.push_back(bT0);
		pt2.push_back(T0);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		tempVols = M.CreatSweepVol(NS, -m_box_L, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		M.Trans(tempVols, -m_box_t, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, m_box_w + 2 * m_box_t, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, m_box_t, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		tempVols = M.CreatSweepVol(NS, m_box_t, z_forward);
		M.Trans(tempVols, -m_box_t, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, m_box_w + 3 * m_box_t, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, -(m_box_L + m_box_t), z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, -(m_box_w + 3 * m_box_t), x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		pt1.clear(), pt2.clear();
		pt1.push_back(bp1t);
		pt1.push_back(p1t);
		pt2.push_back(bc0);
		pt2.push_back(c0);
		Rect = Rect_Create(pt1, pt2);
		NS.clear();
		NS.push_back(Rect);
		pt1.clear(), pt2.clear();
		pt1.push_back(bc0);
		pt1.push_back(c0);
		pt2.push_back(bp1);
		pt2.push_back(p1);
		Rect = Rect_Create(pt1, pt2);
		NS.push_back(Rect);
		tempVols = M.CreatSweepVol(NS, m_box_t, z_forward);
		M.Trans(tempVols, -m_box_t, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, m_box_w + 3 * m_box_t, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, -(m_box_L + m_box_t), z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, -(m_box_w + 3 * m_box_t), x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		tempVols = M.CreatSweepVol(NS, -m_box_L, z_forward);
		M.Trans(tempVols, -m_box_t, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, m_box_w + 3 * m_box_t, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		NS = Tri_cir4_create(m_box_t, m_box_t, 0.4*m_box_t);
		M.Rolate(NS, -PI / 2, x_forward);
		tempVols = M.CreatSweepVol(NS, -m_box_t, y_forward);
		M.Trans(tempVols, -3 * m_box_t / 2, x_forward);
		//M.Trans(tempVols, -m_box_t, y_forward);
		M.Trans(tempVols, 2 * m_box_t, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, -(3 * m_box_t + m_box_L), z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, 3 * m_box_t + m_box_w, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, 3 * m_box_t + m_box_L, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
	
		tempVols = M.CreatSweepVol(NS, box_outer_t / 2, y_forward);
		M.Trans(tempVols, -3 * m_box_t / 2, x_forward);
		M.Trans(tempVols, 2 * m_box_t, z_forward);
		M.Trans(tempVols, y1 - box_outer_t / 2, y_forward);
		v1 = tempVols;
		M.Trans(v1, box_outer_t / 2, y_forward);
		for (auto& v : v1) {
			tempVols.push_back(v);
		}

		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

		M.Trans(tempVols, -(3 * m_box_t + m_box_L), z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, 3 * m_box_t + m_box_w, x_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}
		M.Trans(tempVols, 3 * m_box_t + m_box_L, z_forward);
		for (auto& x : tempVols) {
			reducer.push_back(x);
		}

	}
	varray<SplineVolume> getVols() { return reducer; }

private:
	varray<SplineSurface> Tri_cir4_create(double L, double w, double r) {
		varray<SplineSurface> NS;
		Cirle_Triangle* ct = new Cirle_Triangle(w / 2, L / 2, L / 2, 0, w / 2, r, 0, 0, 0, PI / 2);
		NS.push_back(ct->getSuf());
		ct = new Cirle_Triangle(L / 2, w / 2, w / 2, 0, L / 2, r, L/2, w / 2, 0, -PI);
		NS.push_back(ct->getSuf());
		ct = new Cirle_Triangle(w / 2, L / 2, L / 2, 0, w / 2, r, 0, w, 0, 3*PI / 2);
		NS.push_back(ct->getSuf());
		ct = new Cirle_Triangle(L / 2, w / 2, w / 2, 0, L / 2, r, -L/2, w / 2, 0, 0);
		NS.push_back(ct->getSuf());
		return NS;
	}
	varray<SplineSurface> Trape4_Create(double bL, double sL, double bW, double sW) {
		varray<SplineSurface> trapes;
		trapezoid* t = new trapezoid(sL, bL, (bW - sW) / 2, 0, 0, 0, PI / 2);
		trapes.push_back(t->getSuf());
		t = new trapezoid(sW, bW, (bL - sL) / 2, bL / 2, bW / 2, 0, 0);
		trapes.push_back(t->getSuf());
		t = new trapezoid(sL, bL, (bW - sW) / 2, 0, bW, 0, 3 * PI / 2);
		trapes.push_back(t->getSuf());
		t = new trapezoid(sW, bW, (bL - sL) / 2, -bL / 2, bW / 2, 0, -PI);
		trapes.push_back(t->getSuf());
		return trapes;
	}

	SplineSurface Rect_Create(varray<Vec4>& pt1, varray<Vec4>& pt2) {
		varray<Spline> Ls;
		Spline L;
		L.m_Degree = m_degree;
		L.m_Knots = m_knots;
		SplineSurface NS;
		L.m_CtrlPts.push_back(pt1[0]);
		L.m_CtrlPts.push_back(pt1[0] / 2 + pt1[1] / 2);
		L.m_CtrlPts.push_back(pt1[1]);
		Ls.push_back(L);

		L.m_CtrlPts.clear();
		L.m_CtrlPts.push_back(pt1[0]);
		L.m_CtrlPts.push_back((pt1[0] + pt2[0]) / 2);
		L.m_CtrlPts.push_back(pt2[0]);
		Ls.push_back(L);
		
		L.m_CtrlPts.push_back(pt2[0]);
		L.m_CtrlPts.push_back(pt2[0] / 2 + pt2[1] / 2);
		L.m_CtrlPts.push_back(pt2[1]);
		Ls.push_back(L);
		
		L.m_CtrlPts.clear();
		L.m_CtrlPts.push_back(pt1[1]);
		L.m_CtrlPts.push_back((pt1[1] + pt2[1]) / 2);
		L.m_CtrlPts.push_back(pt2[1]);
		Ls.push_back(L);
		NS.CoonsInterpolate(Ls);
		return NS;
		
	


	}

	//����Բ���Ŵ�mode0����Сmode1
	SplineSurface scaling_cir_creat(varray<Vec4>& pt1,double xt, double yt,double r, double t, int mode) {
		if (pt1.size() == 2) {
			auto p = pt1[1];
			pt1.push_back(p);
			pt1[1] = Model_Solution::NurbsCirle_mid_pt(pt1[0], pt1[2], xt, yt, r);
		
		}
		varray<Vec4> pt2;
		if (mode == 1) {
			double k = (r - t) / r;
			for (auto& p : pt1) {
	
				Vec4 t = { xt + k * (p.x - xt) , yt + k * (p.y - yt) , p.z , p.w };
				pt2.push_back(t);
			}
		}
		else if (mode == 0) {
			double k = (r + t) / r;
			for (auto& p : pt1) {
				Vec4 t = { xt + k * (p.x - xt) , yt + k * (p.y - yt), p.z,p.w };
			}
		}
		SplineSurface NS = Rect_Cir_2create(pt1, pt2);
		return NS;
	}
	void Clear_all() {

	}
	SplineSurface Rect_Cir_2create(varray<Vec4>& pt1, varray<Vec4>&pt2) {
		varray<Spline> Ls;
		Spline L;
		L.m_Degree = m_degree;
		L.m_Knots = m_knots;
		SplineSurface gear;
		L.m_CtrlPts.push_back(pt1[0]);
		L.m_CtrlPts.push_back(pt1[0] / 2 + pt2[0] / 2);
		L.m_CtrlPts.push_back(pt2[0]);
		Ls.push_back(L);
		L.m_CtrlPts.clear();
		L.m_CtrlPts = pt1;
		Ls.push_back(L);
		L.m_CtrlPts.clear();
		L.m_CtrlPts.push_back(pt1[2]);
		L.m_CtrlPts.push_back(pt1[2]/2+pt2[2]/2);
		L.m_CtrlPts.push_back(pt2[2]);
		Ls.push_back(L);
		L.m_CtrlPts.clear();
		L.m_CtrlPts = pt2;
		Ls.push_back(L);
		gear.CoonsInterpolate(Ls);
		return gear;
	}

	SplineSurface Tri_circle_create(varray<Vec4>& pt1, varray<Vec4>& pt2) {
		varray<Spline> Ls;
		Spline L;
		L.m_Degree = m_degree;
		L.m_Knots = m_knots;
		SplineSurface gear;
		L.m_CtrlPts = pt1;
		Ls.push_back(L);
		L.m_CtrlPts.clear();
		L.m_CtrlPts.push_back(pt1[0]);
		L.m_CtrlPts.push_back((pt2[0] + pt2[0])/2);
		L.m_CtrlPts.push_back(pt2[0]);
		Ls.push_back(L);
		L.m_CtrlPts.clear();
		L.m_CtrlPts.push_back(pt1[0]);
		L.m_CtrlPts.push_back((pt1[0] + pt2[0]) / 2);
		L.m_CtrlPts.push_back(pt2[0]);
		Ls.push_back(L);
		L.m_CtrlPts.clear();
		L.m_CtrlPts.push_back(pt1[2]);
		L.m_CtrlPts.push_back((pt1[2] + pt2[1]) / 2);
		L.m_CtrlPts.push_back(pt2[1]);
		Ls.push_back(L);
		gear.CoonsInterpolate(Ls);
		return gear;
	}
	
	//pts1Ϊ�����Ƶ㣬pts2Ϊ�Ҳ���Ƶ㣬��Ϊ���µ���
	SplineSurface gear_create(varray<Vec4>& pts1, varray<Vec4>& pts2) {
		varray<Spline> Ls;
		Spline L;
		L.m_Degree = m_degree;
		L.m_Knots = m_knots;
		L.m_CtrlPts.push_back(pts1[0]);
		L.m_CtrlPts.push_back((pts1[0] + pts2[0]) / 2);
		L.m_CtrlPts.push_back(pts2[0]);
		Ls.push_back(L);
		L.m_CtrlPts.clear();
		L.m_CtrlPts = pts1;
		Ls.push_back(L);
		L.m_CtrlPts.clear();
		L.m_CtrlPts.push_back(pts1[2]);
		L.m_CtrlPts.push_back((pts2[2] + pts1[2]) / 2);
		L.m_CtrlPts.push_back(pts2[2]);
		Ls.push_back(L);
		L.m_CtrlPts.clear();
		L.m_CtrlPts = pts2;
		Ls.push_back(L);
		SplineSurface gear;
		gear.CoonsInterpolate(Ls);
		return gear;
	}

};


//ֱ�����е�ֱ����ĳ�������
//�˴�z<=41,����Բֱ�����ڳݸ�Բֱ��
class Gear_Straight {
private:
	//�߲����
	double m = 3;				//ģ��
	int z = 25;					//����
	double alph = PI / 9;		//ѹ����
	double hax = 1;				//�ݶ���ϵ��
	double cx = 0.25;			//��ϵϵ��
	double B = 10;				//���ֿ��(������)
	double x = 0;				//��λϵ��
	double Dk = 34;				//�����м�׵�ֱ��

	//�м����
	double ha=0;				//�ݶ���
	double hf=0;				//�ݸ���
	double Da=0;				//�ݶ�Բֱ��
	double Db=0;				//��Բֱ��
	double Df=0;				//�ݸ�Բֱ��
	double D=0;					//�ֶ�Բֱ��
	double S = 0;				//�ֶ�Բ�ݺ�
	double r = 0;				//�ݸ���Բ�ǰ뾶

	//�ײ����
	Vec4 Pa1;				//��������ݶ�Բ���㣬���
	Vec4 Pa2;				//ͬ�ϣ��Ҳ�
	Vec4 P1;					//��������ֶ�Բ���㣬���
	Vec4 P2;					//ͬ�ϣ��Ҳ�
	Vec4 Pb1;				//���������Բ���㣬���
	Vec4 Pb2;				//ͬ���Ҳ�
	Vec4 Cirpt;				//������Բ��Բ��
	Vec4 ChamferPt;			//��Բ��Բ��
	double rc;					//������Բ�İ뾶
	Vec4 Pr1;				//Բ����ݸ�Բ����
	Vec4 Pr2;				//ͬ��
	double betak;				//һ���ݵ�Բ�ĽǶȵ�1/2
	Vec4 Pmr1;
	Vec4 Pmr2;	
public:

	//������Ƭ,��һƬΪ������ߣ��ڶ�ƬΪ�����ұ�
	varray<SplineSurface> Gear;

	//������ս��
	varray<SplineVolume> GearVols;

	//��ʼ������
	//mΪģ����zΪ������alphѹ���ǣ�hax�ݶ���ϵ����cx��϶ϵ����B���ֿ�ȣ�xΪ��λϵ����Dk�����м��ֱ��
	void Init(double m , int z , double alph ,double hax ,double  cx ,
				double B,double x,double Dk )
	{
		this->m = m;
		this->z = z;
		this->alph = alph;
		this->hax = hax;
		this->cx = cx;
		this->B = B;
		this->x = x;
		this->Dk = Dk;
	}

	//�����߷���
	double invak(double a)
	{
		return tan(a) - a;
	}

	//��ѹ����
	double calpressangle(double r)
	{
		double rb = Db / 2;
		return acos(rb / r);
	}

	//����Բ���ϵĵ�,�������ڼ������
	Vec4 GetCirclrPts(Vec3 p1,Vec3 p2,Vec3 p3)
	{
		double x1 = p1.x, y1 = p1.y;
		double x2 = p2.x, y2 = p2.y;
		double x3 = p3.x, y3 = p3.y;
		double a = 2 * (x2 - x1);
		double b = 2 * (y2 - y1);
		double c = x2 * x2 + y2 * y2 - x1 * x1 - y1 * y1;
		double d = 2 * (x3 - x2);
		double e = 2 * (y3 - y2);
		double f = x3 * x3 + y3 * y3 - x2 * x2 - y2 * y2;
		double x = (b*f - e * c) / (b*d - e * a);
		double y = (d*c - a * f) / (b*d - e * a);
		double r = sqrt((x - x1)*(x - x1) + (y - y1)*(y - y1));
		Cirpt = { x,y,0 };
		rc = r;
		double k1 =  -1*(x1 - x)/ (y1 - y);
		double b1 = y1 - k1 * x1;
		double k2 = -1 * (x3 - x) / (y3 - y);
		double b2 = y3 - k2 * x3;
		double Xres = (b2 - b1) / (k1 - k2);
		double Yres = k1 * Xres + b1;
		double l = sqrt((Xres - x1)*(Xres - x1) + (Yres - y1)*(Yres - y1));
		double beta0 = atan(l / r);
		return Vec4(Xres, Yres, 0, cos(beta0));
	}


	//��������ϵĵ�
	Vec4 GetGearPts(double R)
	{
		double ai = calpressangle(R);
		double ri = Db / 2 / cos(ai);
		double ceta = invak(ai);
		return Vec4(ri*cos(ceta), ri*sin(ceta), 0);
	}

	//��ȡԲ���������ߵĽ���
	Vec4 GetTangentPts(Vec4 p1,Vec4 p2,Vec4 Cpts,double R)
	{
		double x1 = p1.x, y1 = p1.y;
		double x2 = p2.x, y2 = p2.y;
		double a = Cpts.x, b = Cpts.y;
		double k1 = -(x1 - a) / (y1 - b), k2 = -(x2 - a) / (y2 - b);
		double b1 = y1 - k1 * x1, b2 = y2 - k2 * x2;
		double x = -(b1 - b2) / (k1 - k2);
		double y = k1 * x + b1;
		double r0 = R;
		double l = sqrt((x - x1)*(x - x1) + (y - y1)*(y - y1));
		double delt0 = asin(l/r0);
		return Vec4(x, y, 0, cos(delt0));
	}

	//��ȡֱ����Բ�Ľ���
	Vec4 GetCirTangPts()
	{
		double Sb = S * Db / D + Db * invak(alph);
		double k = Sb / D;
		double rk = Dk / 2;
		return Vec4(rk / (sqrt(1 + k * k)), k*rk / (sqrt(1 + k * k)), 0, 1);
	}

	Vec4 GetSystemMetricsPts(double a, double b, double c, Vec4 p)
	{
		double w = p.w;
		double x0 = p.x, y0 = p.y;
		double x1 = (((b*b - a * a)*x0 - 2 * a*b*y0 - 2 * a*c) / (a*a + b * b));
		double y1 = (((a*a - b * b)*y0 - 2 * a*b*x0 - 2 * b*c) / (a*a + b * b));
		return Vec4(x1, y1, 0,w);
	}

	void GetMParameter()
	{
		D = m * z;
		ha = (hax + x)*m;
		hf = (hax + cx - x)*m;
		Da = D + 2 * ha;
		Db = D * cos(alph);
		Df = D - 2 * hf;
		S = PI * m / 2;
		r = 0.38*m;
		betak = PI / z;
	}


	void GetLParameter()
	{
		double Sb = S * Db / D + Db * invak(alph);
		double k = Sb / D;
		Pb2 = { Db / 2,0,0 };
		P2 = GetGearPts(D/2);
		Pa2 = GetGearPts(Da / 2);
		Pb1 = GetSystemMetricsPts(k, -1, 0, Pb2);
		P1 = GetSystemMetricsPts(k, -1, 0, P2);
		Pa1 = GetSystemMetricsPts(k, -1, 0, Pa2);
	}

	void AdjustPatches()
	{
		double Sb = S * Db / D + Db * invak(alph);
		double delt = atan(Sb / D);
		delt = PI / 2 - delt;
		for (int i = 0; i < Gear.size(); i++)
		{
			for (int j = 0; j < Gear[i].m_CtrlPts.size(); j++)
			{
				double w = Gear[i].m_CtrlPts[j].w;
				Gear[i].m_CtrlPts[j] = Gear[i].m_CtrlPts[j].RotateZ(delt);
				Gear[i].m_CtrlPts[j].w = w;
			}
		}
	}
	
	//���㵹Բ�ǵ��е���м��
	void GetChamferPts()
	{
		double x1 = Cirpt.x;
		double y1 = Cirpt.y;
		double r1 = rc;
		double r2 = Df / 2;
		double a = (r + r1)*(r + r1);
		double b = (r + r2)*(r + r2);
		double c = (a - b-x1*x1-y1*y1)/(-2);
		double a1 = 1+x1*x1/(y1*y1);
		double b1 = -1*2 * c*x1/(y1*y1) ;
		double c1 = c*c/(y1*y1)-b;
		double Xres1 = (-1*b1+sqrt(b1*b1-4*a1*c1)) / (2 * a1);
		double Xres2 = (-1*b1 - sqrt(b1*b1 - 4 * a1*c1)) / (2 * a1);
		double Xres = Xres1 > Xres2 ? Xres1 : Xres2;
		double Yres = (c - x1*Xres)/y1;
		ChamferPt = { Xres,Yres,0 };
		double delt1 = atan(abs(Xres - x1) / abs(Yres - y1));
		Pb2 = { x1 + r1 * sin(delt1),y1 - r1 * cos(delt1),0 };
		double Sb = S * Db / D + Db * invak(alph);
		double k = Sb / D;
		Pb1 = GetSystemMetricsPts(k, -1, 0, Pb2);
		double delt2 = atan(abs(Yres) / abs(Xres));
		Pr2 = { r2*cos(delt2),-r2*sin(delt2),0 };
		Pmr2 = GetTangentPts(Pb2, Pr2, ChamferPt, r);
		Pmr1 = GetSystemMetricsPts(k, -1, 0, Pmr2);
		Pr1 = GetSystemMetricsPts(k, -1, 0, Pr2);
	}
	
	Gear_Straight(double m = 2, int z = 30, double alph = PI / 9, double hax = 1, double  cx = 0.25,
		double B = 10, double x = 0, double Dk = 34) :m(m), z(z), alph(alph), hax(hax), cx(cx),
		B(B), x(x), Dk(Dk) {
		varray<Spline> NL1;
		varray<Spline> NL2;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		GetMParameter();
		GetLParameter();
		GetCirclrPts(Pb2, P2, Pa2);
		GetChamferPts();
		double Sb = S * Db / D + Db * invak(alph);
		double k = Sb / D;
		Vec4 new_p2 = GetCirclrPts(Pb2, P2, Pa2);
		Vec4 new_p1 = GetSystemMetricsPts(k, -1, 0, new_p2);
		P1 = new_p1;
		P2 = new_p2;
		Vec4 Pam1 = (Pa1 + Pa2) / 2;
		Vec4 Pm1 = (P1 + P2) / 2;
		Pm1.w = 1;
		Vec4 Pbm1 = (Pb1 + Pb2) / 2;
		NL1.resize(4);
		NL2.resize(4);
		for (int i = 0; i < 4; i++)
		{
			NL1[i].m_Degree = 2;
			NL1[i].m_Knots = knots;
			NL2[i].m_Degree = 2;
			NL2[i].m_Knots = knots;
		}
		NL1[0].m_CtrlPts.push_back(Pb1);
		NL1[0].m_CtrlPts.push_back(Vec3(Pb1 + Pbm1) / 2);
		NL1[0].m_CtrlPts[1].w = 1;
		NL1[0].m_CtrlPts.push_back(Pbm1);
		NL1[1].m_CtrlPts.push_back(Pb1);
		NL1[1].m_CtrlPts.push_back(P1);
		NL1[1].m_CtrlPts.push_back(Pa1);
		NL1[2].m_CtrlPts.push_back(Pa1);
		NL1[2].m_CtrlPts.push_back(Vec3(Pa1 + Pam1) / 2);
		NL1[2].m_CtrlPts[1].w = 1;
		NL1[2].m_CtrlPts.push_back(Pam1);
		NL1[3].m_CtrlPts.push_back(Pbm1);
		NL1[3].m_CtrlPts.push_back(Pm1);
		NL1[3].m_CtrlPts.push_back(Pam1);
		SplineSurface S1;
		S1.CoonsInterpolate(NL1);
		Gear.push_back(S1);
		NL2[0].m_CtrlPts.push_back(Pbm1);
		NL2[0].m_CtrlPts.push_back(Vec3(Pbm1 + Pb2) / 2);
		NL2[0].m_CtrlPts[1].w = 1;
		NL2[0].m_CtrlPts.push_back(Pb2);
		NL2[1].m_CtrlPts.push_back(Pbm1);
		NL2[1].m_CtrlPts.push_back(Pm1);
		NL2[1].m_CtrlPts.push_back(Pam1);
		NL2[2].m_CtrlPts.push_back(Pam1);
		NL2[2].m_CtrlPts.push_back(Vec3(Pam1 + Pa2) / 2);
		NL2[2].m_CtrlPts[1].w = 1;
		NL2[2].m_CtrlPts.push_back(Pa2);
		NL2[3].m_CtrlPts.push_back(Pb2);
		NL2[3].m_CtrlPts.push_back((P2));
		NL2[3].m_CtrlPts.push_back(Pa2);
		SplineSurface S2;
		S2.CoonsInterpolate(NL2);
		Gear.push_back(S2);
		varray<Spline> NL3, NL4;
		NL3.resize(4);
		NL4.resize(4);
		for (int i = 0; i < 4; i++)
		{
			NL3[i].m_Degree = 2;
			NL3[i].m_Knots = knots;
			NL4[i].m_Degree = 2;
			NL4[i].m_Knots = knots;
		}

		double rk = Dk / 2;
		Vec4 Prm1 = Vec4(rk / (sqrt(1 + k * k)), k*rk / (sqrt(1 + k * k)), 0, 1);

		NL3[0].m_CtrlPts.push_back(Prm1);
		NL3[0].m_CtrlPts.push_back(Vec3(Prm1 + Pr2) / 2);
		NL3[0].m_CtrlPts[1].w = 1;
		NL3[0].m_CtrlPts.push_back(Pr2);
		NL3[1].m_CtrlPts.push_back(Prm1);
		NL3[1].m_CtrlPts.push_back(Vec3(Prm1 + Pbm1) / 2);
		NL3[1].m_CtrlPts[1].w = 1;
		NL3[1].m_CtrlPts.push_back(Pbm1);
		NL3[2].m_CtrlPts.push_back(Pbm1);
		NL3[2].m_CtrlPts.push_back(Vec3(Pbm1 + Pb2) / 2);
		NL3[2].m_CtrlPts[1].w = 1;
		NL3[2].m_CtrlPts.push_back(Pb2);
		NL3[3].m_CtrlPts.push_back(Pr2);
		NL3[3].m_CtrlPts.push_back(Pmr2);
		NL3[3].m_CtrlPts.push_back(Pb2);
		SplineSurface S3;
		S3.CoonsInterpolate(NL3);
		Gear.push_back(S3);

		NL4[0].m_CtrlPts.push_back(Pr1);
		NL4[0].m_CtrlPts.push_back(Vec3(Pr1 + Prm1) / 2);
		NL4[0].m_CtrlPts[1].w = 1;
		NL4[0].m_CtrlPts.push_back(Prm1);
		NL4[1].m_CtrlPts.push_back(Pr1);
		NL4[1].m_CtrlPts.push_back(Pmr1);
		NL4[1].m_CtrlPts.push_back(Pb1);
		NL4[2].m_CtrlPts.push_back(Pb1);
		NL4[2].m_CtrlPts.push_back(Vec3(Pb1 + Pbm1) / 2);
		NL4[2].m_CtrlPts[1].w = 1;
		NL4[2].m_CtrlPts.push_back(Pbm1);
		NL4[3].m_CtrlPts.push_back(Prm1);
		NL4[3].m_CtrlPts.push_back(Vec3(Prm1 + Pbm1) / 2);
		NL4[3].m_CtrlPts[1].w = 1;
		NL4[3].m_CtrlPts.push_back(Pbm1);
		SplineSurface S4;
		S4.CoonsInterpolate(NL4);
		Gear.push_back(S4);
		double rf = Df / 2;
		double ceta1 = betak - atan(k);
		double nk = -1 * tan(ceta1);
		Vec4 Pf2 = { rf / (sqrt(1 + nk * nk)),nk*rf / (sqrt(1 + nk * nk)) ,0 };
		Vec4 Pf1 = GetSystemMetricsPts(k, -1, 0, Pf2);
		Vec4 Pk2 = { rk / (sqrt(1 + nk * nk)),nk*rk / (sqrt(1 + nk * nk)) ,0 };
		Vec4 Pk1 = GetSystemMetricsPts(k, -1, 0, Pk2);
		Vec4 Cpts = { 0,0,0 };
		Vec4 Pmf2 = GetTangentPts(Pr2, Pf2, Cpts, rf);
		Vec4 Pmf1 = GetSystemMetricsPts(k, -1, 0, Pmf2);
		Vec4 Pmk2 = GetTangentPts(Prm1, Pk2, Cpts, rk);
		Vec4 Pmk1 = GetSystemMetricsPts(k, -1, 0, Pmk2);


		varray<Spline> NL6, NL5;
		NL5.resize(4);
		NL6.resize(4);
		for (int i = 0; i < 4; i++)
		{
			NL5[i].m_Degree = 2;
			NL5[i].m_Knots = knots;
			NL6[i].m_Degree = 2;
			NL6[i].m_Knots = knots;
		}
		NL5[0].m_CtrlPts.push_back(Prm1);
		NL5[0].m_CtrlPts.push_back(Pmk2);
		NL5[0].m_CtrlPts.push_back(Pk2);
		NL5[1].m_CtrlPts.push_back(Prm1);
		NL5[1].m_CtrlPts.push_back(Vec3(Prm1 + Pr2) / 2);
		NL5[1].m_CtrlPts[1].w = 1;
		NL5[1].m_CtrlPts.push_back(Pr2);
		NL5[2].m_CtrlPts.push_back(Pr2);
		NL5[2].m_CtrlPts.push_back(Pmf2);
		NL5[2].m_CtrlPts.push_back(Pf2);
		NL5[3].m_CtrlPts.push_back(Pk2);
		NL5[3].m_CtrlPts.push_back(Vec3(Pk2 + Pf2) / 2);
		NL5[3].m_CtrlPts[1].w = 1;
		NL5[3].m_CtrlPts.push_back(Pf2);

		SplineSurface S5;
		S5.CoonsInterpolate(NL5);
		Gear.push_back(S5);

		NL6[0].m_CtrlPts.push_back(Pk1);
		NL6[0].m_CtrlPts.push_back(Pmk1);
		NL6[0].m_CtrlPts.push_back(Prm1);
		NL6[1].m_CtrlPts.push_back(Pk1);
		NL6[1].m_CtrlPts.push_back(Vec3(Pk1 + Pf1) / 2);
		NL6[1].m_CtrlPts[1].w = 1;

		NL6[1].m_CtrlPts.push_back(Pf1);
		NL6[2].m_CtrlPts.push_back(Pr1);
		NL6[2].m_CtrlPts.push_back(Pmf1);
		NL6[2].m_CtrlPts.push_back(Pf1);
		NL6[3].m_CtrlPts.push_back(Prm1);
		NL6[3].m_CtrlPts.push_back(Vec3(Prm1 + Pr1) / 2);
		NL6[3].m_CtrlPts[1].w = 1;

		NL6[3].m_CtrlPts.push_back(Pr1);

		SplineSurface S6;
		S6.CoonsInterpolate(NL6);
		Gear.push_back(S6);

		AdjustPatches();							//������׼��λ��

		varray<SplineSurface> temp = Gear;
		for (int i = 0; i < z - 1; i++)
		{
			for (int j = 0; j < temp.size(); j++)
			{
				for (int k = 0; k < temp[j].m_CtrlPts.size(); k++)
				{
					double w = temp[j].m_CtrlPts[k].w;
					temp[j].m_CtrlPts[k] = temp[j].m_CtrlPts[k].RotateZ(2 * betak);
					temp[j].m_CtrlPts[k].w = w;
				}
			}

			for (int j = 0; j < temp.size(); j++)
			{
				Gear.push_back(temp[j]);
			}
		}
		Model_Solution m0;
		GearVols = m0.CreatSweepVol(Gear, B, z_forward);
	}

	void GeneratePathes(double m = 2, int z = 30, double alph = PI / 9, double hax = 1, double  cx = 0.25,
		double B = 10, double x = 0, double Dk = 34)
	{
		Init(m, z, alph, hax, cx,B, x, Dk);
		varray<Spline> NL1;
		varray<Spline> NL2;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		GetMParameter();
		GetLParameter();
		GetCirclrPts(Pb2, P2, Pa2);
		GetChamferPts();
		double Sb = S * Db / D + Db * invak(alph);
		double k = Sb / D;
		Vec4 new_p2 = GetCirclrPts(Pb2, P2, Pa2);
		Vec4 new_p1 = GetSystemMetricsPts(k,-1,0,new_p2);
		P1 = new_p1;
		P2 = new_p2;
		Vec4 Pam1 = (Pa1 + Pa2) / 2;
		Vec4 Pm1 = (P1 + P2) / 2;
		Pm1.w = 1;
		Vec4 Pbm1 = (Pb1 + Pb2) / 2;
		NL1.resize(4);
		NL2.resize(4);
		for (int i = 0; i < 4; i++)
		{
			NL1[i].m_Degree = 2;
			NL1[i].m_Knots = knots;
			NL2[i].m_Degree = 2;
			NL2[i].m_Knots = knots;
		}
		NL1[0].m_CtrlPts.push_back(Pb1);
		NL1[0].m_CtrlPts.push_back(Vec3(Pb1 + Pbm1) / 2);
		NL1[0].m_CtrlPts[1].w = 1;
		NL1[0].m_CtrlPts.push_back(Pbm1);
		NL1[1].m_CtrlPts.push_back(Pb1);
		NL1[1].m_CtrlPts.push_back(P1);
		NL1[1].m_CtrlPts.push_back(Pa1);
		NL1[2].m_CtrlPts.push_back(Pa1);
		NL1[2].m_CtrlPts.push_back(Vec3(Pa1 + Pam1) / 2);
		NL1[2].m_CtrlPts[1].w = 1;
		NL1[2].m_CtrlPts.push_back(Pam1);
		NL1[3].m_CtrlPts.push_back(Pbm1);
		NL1[3].m_CtrlPts.push_back(Pm1);
		NL1[3].m_CtrlPts.push_back(Pam1);
		SplineSurface S1; 
		S1.CoonsInterpolate(NL1);
		Gear.push_back(S1);
		NL2[0].m_CtrlPts.push_back(Pbm1);
		NL2[0].m_CtrlPts.push_back(Vec3(Pbm1 + Pb2) / 2);
		NL2[0].m_CtrlPts[1].w = 1;
		NL2[0].m_CtrlPts.push_back(Pb2);
		NL2[1].m_CtrlPts.push_back(Pbm1);
		NL2[1].m_CtrlPts.push_back(Pm1);
		NL2[1].m_CtrlPts.push_back(Pam1);
		NL2[2].m_CtrlPts.push_back(Pam1);
		NL2[2].m_CtrlPts.push_back(Vec3(Pam1 + Pa2) / 2);
		NL2[2].m_CtrlPts[1].w = 1;
		NL2[2].m_CtrlPts.push_back(Pa2);
		NL2[3].m_CtrlPts.push_back(Pb2);
		NL2[3].m_CtrlPts.push_back((P2));
		NL2[3].m_CtrlPts.push_back(Pa2);
		SplineSurface S2;
		S2.CoonsInterpolate(NL2);
		Gear.push_back(S2);
		varray<Spline> NL3, NL4;
		NL3.resize(4);
		NL4.resize(4);
		for (int i = 0; i < 4; i++)
		{
			NL3[i].m_Degree = 2;
			NL3[i].m_Knots = knots;
			NL4[i].m_Degree = 2;
			NL4[i].m_Knots = knots;
		}
		
		double rk = Dk / 2;
		Vec4 Prm1 = Vec4(rk / (sqrt(1 + k * k)), k*rk / (sqrt(1 + k * k)), 0, 1);

		NL3[0].m_CtrlPts.push_back(Prm1);
		NL3[0].m_CtrlPts.push_back(Vec3(Prm1 + Pr2) / 2);
		NL3[0].m_CtrlPts[1].w = 1;
		NL3[0].m_CtrlPts.push_back(Pr2);
		NL3[1].m_CtrlPts.push_back(Prm1);
		NL3[1].m_CtrlPts.push_back(Vec3(Prm1 + Pbm1) / 2);
		NL3[1].m_CtrlPts[1].w = 1;
		NL3[1].m_CtrlPts.push_back(Pbm1);
		NL3[2].m_CtrlPts.push_back(Pbm1);
		NL3[2].m_CtrlPts.push_back(Vec3(Pbm1 + Pb2) / 2);
		NL3[2].m_CtrlPts[1].w = 1;
		NL3[2].m_CtrlPts.push_back(Pb2);
		NL3[3].m_CtrlPts.push_back(Pr2);
		NL3[3].m_CtrlPts.push_back(Pmr2);
		NL3[3].m_CtrlPts.push_back(Pb2);
		SplineSurface S3;
		S3.CoonsInterpolate(NL3);
		Gear.push_back(S3);
		
		NL4[0].m_CtrlPts.push_back(Pr1);
		NL4[0].m_CtrlPts.push_back(Vec3(Pr1 + Prm1) / 2);
		NL4[0].m_CtrlPts[1].w = 1;
		NL4[0].m_CtrlPts.push_back(Prm1);
		NL4[1].m_CtrlPts.push_back(Pr1);
		NL4[1].m_CtrlPts.push_back(Pmr1);
		NL4[1].m_CtrlPts.push_back(Pb1);
		NL4[2].m_CtrlPts.push_back(Pb1);
		NL4[2].m_CtrlPts.push_back(Vec3(Pb1 + Pbm1) / 2);
		NL4[2].m_CtrlPts[1].w = 1;
		NL4[2].m_CtrlPts.push_back(Pbm1);
		NL4[3].m_CtrlPts.push_back(Prm1);
		NL4[3].m_CtrlPts.push_back(Vec3(Prm1 + Pbm1) / 2);
		NL4[3].m_CtrlPts[1].w = 1;
		NL4[3].m_CtrlPts.push_back(Pbm1);
		SplineSurface S4;
		S4.CoonsInterpolate(NL4);
		Gear.push_back(S4);
		double rf = Df / 2;
		double ceta1 = betak - atan(k);
		double nk = -1 * tan(ceta1);
		Vec4 Pf2 = { rf / (sqrt(1 + nk * nk)),nk*rf / (sqrt(1 + nk * nk)) ,0 };
		Vec4 Pf1 = GetSystemMetricsPts(k, -1, 0, Pf2);
		Vec4 Pk2 = { rk / (sqrt(1 + nk * nk)),nk*rk / (sqrt(1 + nk * nk)) ,0 };
		Vec4 Pk1 = GetSystemMetricsPts(k, -1, 0, Pk2);
		Vec4 Cpts = { 0,0,0 };
		Vec4 Pmf2 = GetTangentPts(Pr2, Pf2, Cpts, rf);
		Vec4 Pmf1 = GetSystemMetricsPts(k, -1, 0, Pmf2);
		Vec4 Pmk2 = GetTangentPts(Prm1, Pk2, Cpts, rk);
		Vec4 Pmk1 = GetSystemMetricsPts(k, -1, 0, Pmk2);


		varray<Spline> NL6, NL5;
		NL5.resize(4);
		NL6.resize(4);
		for (int i = 0; i < 4; i++)
		{
			NL5[i].m_Degree = 2;
			NL5[i].m_Knots = knots;
			NL6[i].m_Degree = 2;
			NL6[i].m_Knots = knots;
		}
		NL5[0].m_CtrlPts.push_back(Prm1);
		NL5[0].m_CtrlPts.push_back(Pmk2);
		NL5[0].m_CtrlPts.push_back(Pk2);
		NL5[1].m_CtrlPts.push_back(Prm1);
		NL5[1].m_CtrlPts.push_back(Vec3(Prm1+Pr2)/2);
		NL5[1].m_CtrlPts[1].w = 1;
		NL5[1].m_CtrlPts.push_back(Pr2);
		NL5[2].m_CtrlPts.push_back(Pr2);
		NL5[2].m_CtrlPts.push_back(Pmf2);
		NL5[2].m_CtrlPts.push_back(Pf2);
		NL5[3].m_CtrlPts.push_back(Pk2);
		NL5[3].m_CtrlPts.push_back(Vec3(Pk2+Pf2)/2);
		NL5[3].m_CtrlPts[1].w = 1;
		NL5[3].m_CtrlPts.push_back(Pf2);

		SplineSurface S5;
		S5.CoonsInterpolate(NL5);
		Gear.push_back(S5);

		NL6[0].m_CtrlPts.push_back(Pk1);
		NL6[0].m_CtrlPts.push_back(Pmk1);
		NL6[0].m_CtrlPts.push_back(Prm1);
		NL6[1].m_CtrlPts.push_back(Pk1);
		NL6[1].m_CtrlPts.push_back(Vec3(Pk1+Pf1)/2);
		NL6[1].m_CtrlPts[1].w = 1;

		NL6[1].m_CtrlPts.push_back(Pf1);
		NL6[2].m_CtrlPts.push_back(Pr1);
		NL6[2].m_CtrlPts.push_back(Pmf1);
		NL6[2].m_CtrlPts.push_back(Pf1);
		NL6[3].m_CtrlPts.push_back(Prm1);
		NL6[3].m_CtrlPts.push_back(Vec3(Prm1+Pr1)/2);
		NL6[3].m_CtrlPts[1].w = 1;

		NL6[3].m_CtrlPts.push_back(Pr1);

		SplineSurface S6;
		S6.CoonsInterpolate(NL6);
		Gear.push_back(S6);

		AdjustPatches();							//������׼��λ��

		varray<SplineSurface> temp = Gear;
		for (int i = 0; i < z-1; i++)
		{
			for (int j = 0; j < temp.size(); j++)
			{
				for (int k = 0; k < temp[j].m_CtrlPts.size(); k++)
				{
					double w = temp[j].m_CtrlPts[k].w;
					temp[j].m_CtrlPts[k] = temp[j].m_CtrlPts[k].RotateZ(2 * betak);
					temp[j].m_CtrlPts[k].w = w;
				}
			}

			for (int j = 0; j < temp.size(); j++)
			{
				Gear.push_back(temp[j]);
			}
		}
		Model_Solution m0;
		GearVols = m0.CreatSweepVol(Gear, B, z_forward);
	}
	varray<SplineVolume> getGear() { return GearVols; }
	void show()
	{
		cout << Cirpt.x << "," << Cirpt.y << endl;
		cout << rc << endl;
	}

	
};


class BearingBlock2
{
public:
	varray<SplineVolume> BearingBlock_Vols;
	varray<SplineVolume> getVol() { return BearingBlock_Vols; }
private:
	//����
	double L1;
	double L2;
	double H1;
	double H2;
	double r1;
	double r2;
	double t;

	//��ɿ�
	
	varray<Cirle_Triangle> Cts;
	varray<RecTangle> Rts;

	Model_Solution M;

	void InitCts()
	{
		//���沿��
		varray<SplineSurface> BearingBlock_Surf;
		varray<SplineVolume> Vols;
		Cirle_Triangle C0;
		C0 = Cirle_Triangle(H2 / 2, L2 / 2, L2 / 2, 0, H2 / 2, r2, 0, 0, 0, PI / 2);
		//M.Rolate(C0.CircuSquare, PI / 2, x_forward);
		Cts.push_back(C0);
		C0 = Cirle_Triangle(H2 / 2, L2 / 2, L2 / 2, 0, H2 / 2, r2, -L2/2, H2/2, 0, 0);
		Cts.push_back(C0);
		C0 = Cirle_Triangle(H2 / 2, L2 / 2, L2 / 2, 0, H2 / 2, r2, 0, H2, 0, -PI/2);
		Cts.push_back(C0);
		C0 = Cirle_Triangle(H2 / 2, L2 / 2, L2 / 2, 0, H2 / 2, r2, L2/2, H2/2, 0, -PI);
		Cts.push_back(C0);

		for (auto block : Cts)
		{
			BearingBlock_Surf.push_back(block.CircuSquare);
		}

		Vols = M.CreatSweepVol(BearingBlock_Surf, t, z_forward);

		for (auto Vol : Vols)
		{
			BearingBlock_Vols.push_back(Vol);
		}
		
		//���Ҳ���
		BearingBlock_Surf.clear();
		Vols.clear();
		Cts.clear();

		C0 = Cirle_Triangle(t / 2, L1 / 2, L1 / 2, 0, t / 2, r1, -(L1 + L2) / 2, 0, 0, PI / 2);
		Cts.push_back(C0);

		C0 = Cirle_Triangle(t / 2, L1 / 2, L1 / 2, 0, t / 2, r1, -(L1 + L2) / 2 - L1/2, t/2, 0, 0);
		Cts.push_back(C0);

		C0 = Cirle_Triangle(t / 2, L1 / 2, L1 / 2, 0, t / 2, r1, -(L1 + L2) / 2, t, 0, -PI / 2);
		Cts.push_back(C0);

		C0 = Cirle_Triangle(t / 2, L1 / 2, L1 / 2, 0, t / 2, r1, -(L1 + L2) / 2 + L1/2, t/2, 0,-PI);
		Cts.push_back(C0);

		for (auto block : Cts)
		{
			BearingBlock_Surf.push_back(block.CircuSquare);
		}

		M.Rolate(BearingBlock_Surf, PI / 2, 1);
		Vols = M.CreatSweepVol(BearingBlock_Surf, H1, y_backward);

		for (auto Vol : Vols)
		{
			BearingBlock_Vols.push_back(Vol);
		}

		BearingBlock_Surf.clear();
		Vols.clear();
		Cts.clear();

		C0 = Cirle_Triangle(t / 2, L1 / 2, L1 / 2, 0, t / 2, r1, (L1 + L2) / 2, 0, 0, PI / 2);
		Cts.push_back(C0);

		C0 = Cirle_Triangle(t / 2, L1 / 2, L1 / 2, 0, t / 2, r1, (L1 + L2) / 2 - L1 / 2, t / 2, 0, 0);
		Cts.push_back(C0);

		C0 = Cirle_Triangle(t / 2, L1 / 2, L1 / 2, 0, t / 2, r1, (L1 + L2) / 2, t, 0, -PI / 2);
		Cts.push_back(C0);

		C0 = Cirle_Triangle(t / 2, L1 / 2, L1 / 2, 0, t / 2, r1, (L1 + L2) / 2 + L1 / 2, t / 2, 0, -PI);
		Cts.push_back(C0);

		for (auto block : Cts)
		{
			BearingBlock_Surf.push_back(block.CircuSquare);
		}

		M.Rolate(BearingBlock_Surf, PI / 2, 1);
		Vols = M.CreatSweepVol(BearingBlock_Surf, H1, y_backward);

		for (auto Vol : Vols)
		{
			BearingBlock_Vols.push_back(Vol);
		}

		
	}
	void InitRts()
	{
		varray<SplineSurface> BearingBlock_Surf;
		varray<SplineVolume> Vols;
		RecTangle Rect;
		Rect = RecTangle(L2, H1, 0, -H1, 0, PI / 2);
		BearingBlock_Surf.push_back(Rect.Rect);
		Vols = M.CreatSweepVol(BearingBlock_Surf, t, z_forward);
		for (auto Vol : Vols)
			BearingBlock_Vols.push_back(Vol);

	}
public:
	void BearingBlock(double L1 = 5, double L2 = 10, double H1 = 5, double H2 = 10, double r1 = 2, double r2 = 4, double t = 5)
	{
		this->H1 = H1;
		this->H2 = H2;
		this->L1 = L1;
		this->L2 = L2;
		this->r1 = r1;
		this->r2 = r2;
		this->t = t;
		InitRts();
		InitCts();
	}

};


//���������ģ�飬�������ÿ���Ӵ����Ƿ����
class TestBolcks
{
public:
	static double distance(Vec4 p1, Vec4 p2)
	{
		return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
	}

	static void get_Align(varray<SplineVolume>& Vols) {
		varray<map<Vec4, int>> Vmaps;
		varray<varray<map<Vec4, int>>> Smaps;
		varray<varray<map<Vec4, int>>> Bmaps;
		Model_Solution M;
		//��װVmaps
		int NVidx = 0;
		for (auto& NV : Vols) {
			map<Vec4, int> Vmap;
			for (int i = 0; i < NV.m_CtrlPts.size(); i++) {
				if(Vmap.find(NV.m_CtrlPts[i]) != Vmap.end())
				{
					cout << "idx"<< NVidx << "ERROR" << endl;
				}
				Vmap[NV.m_CtrlPts[i]] = i;
			}
			Vmaps.push_back(Vmap);
			NVidx++;
		}
		varray<varray<SplineSurface>> NS = M.GetSurfaces(Vols);
		//��װBmaps��Smaps
		for (int i = 0; i < NS.size(); i++) {
			varray<map<Vec4, int>> Bmap;
			varray<map<Vec4, int>>Smap;
		
			for (int j = 0; j < NS[i].size(); j++) {
				map<Vec4, int>temp;
				temp[NS[i][j].GetSurFacePoint(0, 0)] = 0;
				temp[NS[i][j].GetSurFacePoint(1, 0)] = 1;
				temp[NS[i][j].GetSurFacePoint(0, 1)] = 2;
				temp[NS[i][j].GetSurFacePoint(1, 1)] = 3;
				Bmap.push_back(temp);
				temp.clear();
				for (int k = 0; k < NS[i][j].m_CtrlPts.size(); k++) {
					temp[NS[i][j].m_CtrlPts[k]] = k;
				}
				Smap.push_back(temp);

			}
			Bmaps.push_back(Bmap);
			Smaps.push_back(Smap);
		}
		
		for (int Vidx0 = 0; Vidx0 < Bmaps.size() - 1; Vidx0++)
		{	
			for (int Sufidx0 = 0; Sufidx0 < 6; Sufidx0++) {

				for (int Vidx1 = Vidx0 + 1; Vidx1 < Bmaps.size(); Vidx1++) {
					for (int Sufidx1 = 0; Sufidx1 < 6; Sufidx1++) {
						int count = 0;
						for (auto& n : Bmaps[Vidx1][Sufidx1]) {
							if (Bmaps[Vidx0][Sufidx0].find(n.first) == Bmaps[Vidx0][Sufidx0].end())
							{
								break;
							}
							count++;
						}
						if (count != 4)
							continue;
						else {
						/*	cout << "patch0" << Vidx0 << "," << Sufidx0 << endl;
							cout << "patch1" << Vidx1 << "," << Sufidx1 << endl;
							cout << endl;*/
							//��������Ӵ���
							for (auto& n : Smaps[Vidx1][Sufidx1]) {
								if (Smaps[Vidx0][Sufidx0].find(n.first) == Smaps[Vidx0][Sufidx0].end()) {
									double min_dis = INT_MAX;
									double cur_dis = 0;
									Vec4 ans;
									for (auto m : Smaps[Vidx0][Sufidx0])
									{
										cur_dis = distance(n.first, m.first);
										if (cur_dis > min_dis)
											continue;
										min_dis = cur_dis;
										ans = m.first;
									}
				/*					cout << endl;
									cout << Vols[Vidx0].m_CtrlPts[Vmaps[Vidx0][ans]].x << " " << Vols[Vidx0].m_CtrlPts[Vmaps[Vidx0][ans]].y << " " << Vols[Vidx0].m_CtrlPts[Vmaps[Vidx0][ans]].z << endl;
									cout << n.first.x << " " << n.first.y << " " << n.first.z << endl;
				*/					Vols[Vidx0].m_CtrlPts[Vmaps[Vidx0][ans]] = n.first;
								}
							/*	else {
									Smaps[Vidx0][Sufidx0].erase(n.first);
									Smaps[Vidx1][Sufidx1].erase(n.first);
								}*/
							}
						}
					}
				}
			}
		}
		
	/*	for (const auto& n : Bmaps) {
			for (const auto& m : n) {
				for (const auto& p : m) {
					cout << p.first.x << "," << p.first.y << "," << p.first.z << "\t";
				}
				cout << endl;
			}
			cout << endl;
		}*/
	}
	

	struct ListNode {
		Vec4 val;
		ListNode *next;
		ListNode() { next = nullptr; };
	};

	class pList {
	public:
		//ListNode* pList;
		map<Vec4, int> ptMap;
		varray<SplineVolume> Vols;
		varray<SplineSurface> SS;
		varray<int> ptidx;
		varray<varray<int>> CtrlIDs;
	public:
		
		void writeCtrlptID(const varray<SplineVolume>& Vols) {
			
		}
		void OutputParaVolumeDataTxt(const varray<SplineVolume>& vol, const string& path)
		{
			this->Vols = vol;
			string filename = path + "ctrlpts.txt";
			string fileIDname = path + "ctrlptsIdx.txt";
			ofstream ofs(filename);
			ofs << "PN " << " " << Vols.size() << "\n";  //patch number
			for (int i = 0; i < Vols.size(); i++)
			{
				ofs << "PI" << " " << i << "\n";   //��ǰƬ��id�š�

				SplineVolume& sv = Vols.at(i);
				ofs << "OD" << " " << sv.m_uDegree + 1 << " "
					<< sv.m_vDegree + 1 << " " << sv.m_wDegree + 1 << "\n"; //Order

				ofs << "UK" << " " << sv.m_uKnots.size() << "\n";
				for (int i = 0; i < sv.m_uKnots.size(); i++)
					ofs << sv.m_uKnots.at(i) << " ";   //�ڵ�����
				ofs << "\n";
				ofs << "VK" << " " << sv.m_vKnots.size() << "\n";
				for (int i = 0; i < sv.m_vKnots.size(); i++)
					ofs << sv.m_vKnots.at(i) << " ";   //�ڵ�����
				ofs << "\n";
				ofs << "WK" << " " << sv.m_wKnots.size() << "\n";
				for (int i = 0; i < sv.m_wKnots.size(); i++)
					ofs << sv.m_wKnots.at(i) << " ";   //�ڵ�����
				ofs << "\n";

				ofs << "CP" << " " << sv.m_uNum << " " << sv.m_vNum << " " << sv.m_wNum << "\n"; //control point number
				for (int i = 0; i < sv.m_CtrlPts.size(); i++)
					ofs << sv.m_CtrlPts.at(i).x << " " << sv.m_CtrlPts.at(i).y << " " << sv.m_CtrlPts.at(i).z << " " << sv.m_CtrlPts.at(i).w << "\n";
			}
			ofs.close();

			ofstream idofs(fileIDname);
			getVolidx(Vols);
			ofs << "PN " << " " << Vols.size() << "\n";  //patch number
			int idx = 0;
			for (int i = 0; i < Vols.size(); i++)
			{
				int idxnum = Vols[i].m_CtrlPts.size();
				idofs << "PN" << " " << i << "\n";
				idofs << "PI" << " " << idxnum << "\n";
				for (int j = 0; j < idxnum; j++)
				{
					idofs << ptidx[idx++] << " ";
				}
				idofs << "\n";
			}
			idofs.close();
		}
		
		void OutputParaSurfaceDataTxt(const varray<SplineSurface>& ss, const string& path)
		{
			this->SS = ss;
			string filename = path + "ctrlpts.txt";
			string fileIDname = path + "ctrlptsIdx.txt";
			ofstream ofs(filename);
			ofs << "PN " << " " << SS.size() << "\n";  //patch number
			for (int i = 0; i < SS.size(); i++)
			{
				ofs << "PI" << " " << i << "\n";   //��ǰƬ��id�š�

				SplineSurface& s = SS.at(i);
				ofs << "OD" << " " << s.m_uDegree + 1 << " "
					<< s.m_vDegree + 1 << "\n"; //Order

				ofs << "UK" << " " << s.m_uKnots.size() << "\n";
				for (int i = 0; i < s.m_uKnots.size(); i++)
					ofs << s.m_uKnots.at(i) << " ";   //�ڵ�����
				ofs << "\n";
				ofs << "VK" << " " << s.m_vKnots.size() << "\n";
				for (int i = 0; i < s.m_vKnots.size(); i++)
					ofs << s.m_vKnots.at(i) << " ";   //�ڵ�����
				ofs << "\n";

				ofs << "CP" << " " << s.m_uNum << " " << s.m_vNum << "\n"; //control point number
				for (int i = 0; i < s.m_CtrlPts.size(); i++)
					ofs << s.m_CtrlPts.at(i).x << " " << s.m_CtrlPts.at(i).y << " " << s.m_CtrlPts.at(i).z << " " << s.m_CtrlPts.at(i).w << "\n";
			}
			ofs.close();

			ofstream idofs(fileIDname);
			
			getSurfidx(SS);
			ofs << "PN " << " " << SS.size() << "\n";  //patch number
			int idx = 0;
			for (int i = 0; i < SS.size(); i++)
			{
				int idxnum = SS[i].m_CtrlPts.size();
				idofs << "PN" << " " << i << "\n";
				idofs << "PI" << " " << idxnum << "\n";
				for (int j = 0; j < idxnum; j++)
				{
					idofs << ptidx[idx++] << " ";
				}
				idofs << "\n";
			}
			idofs.close();
		}

		void getVolidx(const varray<SplineVolume>& Vols)
		{
			ptidx.clear();
			this->Vols = Vols;
			int idx = 0;
			for (auto vol : Vols)
			{
				for (auto pt : vol.m_CtrlPts)
				{
					if (ptMap.find(pt) == ptMap.end())
					{
						ptMap[pt] = idx++;
					}
				}
			}
			for (auto vol : Vols)
			{
				for (auto pt : vol.m_CtrlPts)
				{
					ptidx.push_back(ptMap[pt]);
				}
			}
		}

		void getSurfidx(const varray<SplineSurface>& Surfs)
		{
			ptidx.clear();
			this->SS = Surfs;
			int idx = 0;
			for (auto s : SS)
			{
				for (auto pt : s.m_CtrlPts)
				{
					if (ptMap.find(pt) == ptMap.end())
					{
						ptMap[pt] = idx++;
					}
				}
			}
			for (auto s : SS)
			{
				for (auto pt : s.m_CtrlPts)
				{
					ptidx.push_back(ptMap[pt]);
				}
			}
		}

		void showdata()
		{
			for (int i = 0; i < Vols.size(); i++)
			{
				for (int j = 0; j < Vols[i].m_CtrlPts.size(); j++)
				{
					//cout << Vols[i].m_CtrlPts[j].x << " " << Vols[i].m_CtrlPts[j].y << " " << Vols[i].m_CtrlPts[j].z << "\t ";
					//cout << ptidx[i*Vols[i].m_CtrlPts.size()+j] << "\t";
					cout << ptMap[Vols[i].m_CtrlPts[j]]<< "\t";
				}
				cout << endl;
			}
		}

		varray<varray<int>> getfaceidx(const varray<SplineVolume>& Vols,const varray<SplineSurface>& NSF)
		{
			getVolidx(Vols);
			varray<varray<int>> res;
			varray<int>cur_res;
			map<Vec4, int> tidx = ptMap;
			int num = 0;
			for (auto NS : NSF) {
				for (auto NSpt : NS.m_CtrlPts)
				{
					if (tidx.find(NSpt) != tidx.end())
					{
						num = tidx[NSpt];
						cur_res.push_back(num);
						tidx.erase(NSpt);
					}
				}
				res.push_back(cur_res);
				cur_res.clear();
			}
			return res;
		}

		varray<int> getfaceidx(const varray<SplineVolume>& Vols, const SplineSurface& NSF)
		{
			varray<int> res;
			getVolidx(Vols);
			res = getfaceidx(NSF);
			return res;

		}

		varray<int> getfaceidx(const SplineSurface& NS)
		{
			map<Vec4, int> tidx = ptMap;
			varray<int> res;
			int num = 0;
			for (auto NSpt : NS.m_CtrlPts)
			{
				if (tidx.find(NSpt) != tidx.end())
				{
					num = tidx[NSpt];
					res.push_back(num);
					tidx.erase(NSpt);
				}
			}
			return res;
		}

		//varray<SplineVolume>& operator = ()

	};
	

};


//����
class Mspring
{
public:
	Mspring() {}

	//���캯��
	Mspring(double D, double d, double t, double n);

	//��ʼ���߲����
	void InitHighParameter(double D, double d, double t, int n);

	//���NURBS��
	varray<SplineVolume> GetMspring();


private:
	varray<Feature_Vol> m_fvs;			//�����弯��
	Model_Solution ms;

private:
	//�߲����
	double h_D;							//�����о�(ֱ��)
	double h_d;							//����Բֱ��
	double h_t;							//�ھ�
	int h_n;							//��ЧȦ��

	//�в����(�ϼ�,����ӳ�乹��)

	//�Ͳ����(�ϼ򵥣�����ӳ�乹��)
};

//����(�ʷְ汾)
class Gear_Wheel
{
public:
	//Gear_Wheel() {}

	Gear_Wheel(double m = 2, int z = 30, double alph = PI / 9, double hax = 1, double  cx = 0.25,
		double B = 8, double x = 0, double Dk = 32);

	//Gear_Straight(double m = 2, int z = 30, double alph = PI / 9, double hax = 1, double  cx = 0.25,
	//	double B = 10, double x = 0, double Dk = 34) :m(m), z(z), alph(alph), hax(hax), cx(cx),
	//	B(B), x(x), Dk(Dk)

	//�������NURBS��
	varray<SplineVolume> GetGearVols();

	
private:
	varray<Feature_Vol> m_fvs;			//�����弯��
	Model_Solution ms;

private:
	double m = 3;				//ģ��
	int z = 25;					//����
	double alph = PI / 9;		//ѹ����
	double hax = 1;				//�ݶ���ϵ��
	double cx = 0.25;			//��ϵϵ��
	double B = 8;				//���ֿ��(������)
	double x = 0;				//��λϵ��
	double Dk = 32;				//�����м�׵�ֱ��

};

class Box_Quad
{
public:

	Box_Quad(double L = 9, double L1 = 8, double L2 = 2, double L3 = 2, double W = 6, double H = 5, double H1 = 4, double SL1 = 1.5, double SL2 = 0.5,
		double r1 = 1.5, double r2 = 2.5, double r3 = 1, double r4 = 0.15, double t1 = 0.5, double t2 = 0.6, double tk1 = 0.5, double tk2 = 0.5, int n = 4);

	//����ӳ��
	void Map1();
	void Map2();

	//�������
	varray<SplineVolume> GetBox();

private:
	varray<Feature_Vol> m_fvs;			//�����弯��
	Model_Solution ms;

private:
	//�߲����
	double h_L;
	double h_L1;
	double h_L2;
	double h_L3;
	double h_W;
	double h_H;
	double h_H1;
	double h_SL1;
	double h_SL2;
	double h_r1;
	double h_r2;
	double h_r3;
	double h_r4;
	double h_t1;
	double h_t2;
	double h_tk1;
	double h_tk2;
	int h_n;

	//�в����
	double m_L;
	double m_L1;
	double m_L2;
	double m_L3;
	double m_L4;//ԭ�㵽r3Բ��X�������
	double m_W;
	double m_W1;//�װ��
	double m_H;
	double m_H1;
	double m_H2;//�����κ�r3Բ�ĸ߶�
	double m_SL1;
	double m_SL2;
	double m_r1;
	double m_r2;
	double m_r3;
	double m_r4;
	double m_t1;
	double m_t2;
	double m_tk1;
	double m_tk2;

	//�Ͳ����
	//�ݿ�����λ������
	Vec3 l_p1;
	Vec3 l_p2;
	Vec3 l_p3;
	Vec3 l_p4;
	Vec3 l_p5;
	Vec3 l_p6;
	Vec3 l_p7;
	Vec3 l_p8;

	//�װ�y����ֱ������
	Vec3 l_p9;
	Vec3 l_p10;
};


//�����
class Bearing_Chock
{
public:
	Bearing_Chock(double L = 80, double L1 = 12, double L2 = 40, double L3 = 8, double W = 30, double W1 = 15, double W2 = 22, double W3 = 8, double W4 = 30, double W5 = 17,
		double H1 = 45, double H2 = 5, double d1 = 6, double r1 = /*7.5*/5.5, double r2 = /*15*/12, double sita = 10.0 / 180 * PI);

	//����ӳ��
	void Map1();
	void Map2();

	//��������
	varray<SplineVolume> GetBearing();
private:
	varray<Feature_Vol> m_fvs;			//�����弯��
	Model_Solution ms;

private:
	//�߲����
	double h_L;
	double h_L1;
	double h_L2;
	double h_L3;
	double h_W;
	double h_W1;
	double h_W2;
	double h_W3;
	double h_W4;
	double h_W5;
	double h_H1;
	double h_H2;
	double h_d1;
	double h_r1;
	double h_r2;
	double h_sita;		//Բ����ת�Ƕ�

	//�в����
	double m_L;
	double m_L1;
	double m_L2;
	double m_L3;
	double m_W;
	double m_W1;
	double m_W2;
	double m_W3;
	double m_W4;
	double m_W5;
	double m_H1;
	double m_H2;
	double m_d1;
	double m_r1;
	double m_r2;

	double m_L4;		//��Բ�����������
	double m_L5;		//����γ�
	double m_W6;		//ԲͲ���������
	double m_H3;		//�װ������ԲͲԲ�ľ���
	double m_sita_1;	//�ײ�Բ���νǶ�
	double m_sita_2;	//��Բ���νǶ�
	double m_sita_3;	//����Բ���νǶ�

	//�Ͳ����
	Vec3 l_p1;		//�м����(��)���ĵ�����
	Vec3 l_p2;		//�м����(��)���ĵ�����
	Vec3 l_p3;		//��������ĵ�����
	Vec3 l_p4;		//��Բ���ĵ�����
	Vec3 l_p5;		//ԲͲ���ĵ�����

};
