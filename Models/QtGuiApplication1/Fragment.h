//////////////////��ģ�͵���Ƭ��////////////////////////
///////////////2020-3-1 193780714///////////////////////
#pragma once

#include"stdafx.h"

using std::min;
using std::max;
using std::vector;

class Fragment
{
public:
	Fragment();
	Fragment(float a, float b, float c,float begin,float end,float degree);//a,b,c��Ӧƽ���һ��������̣�begin��end��Ӧƽ��Ŀ�ʼ����ʼλ�ã�degree��Ӧϸ�ֶ�
	
	~Fragment();
	//vector<float> ReturnCurve();//���صȲ��ߵĵ�
	vector<float> ReturnSupportPoints();//����֧�ŵ�
	//vector<float> ReturnFourBoundary();//���ؽ���߽�
	//vector<float> ReturnSurface();//�����ͽ���

	bool intersectionLinePlane(point3d p1, point3d p2, vector<float> & temp,  float D);  //���㽻��
	vector<vector<float>> getPoints(varray < varray<varray<point3d>>> points); //�����Ƭ��
	vector<float> getSupportPoints(); //���֧�ŵ�
	bool SetInOrOut(point3d p,int n);//���Ƿ�����ƽ����
	point3d getShadow(point3d p, float D);
protected:
	//�������㺯��
	point3d alf_mult_point(double alf, point3d p);

	//����
	//void QuadsAndHexs();                    //������Ƭ�������ģ��
	//void SolvingSection();                  //�Ȳε�Ԫ����
	//void SolvingSupport();                  //֧�ŵ����

	//
	//void ReturnSupport();                                             //����֧�ŵ�
	//
	//void GetSupporPoints();                                           //���֧�ŵ�
	//
	//void FindBoundaryCurve();                                         //���߽�����
	//void ClearHex();                                                  //��յ�Ԫ
	//void FindBoundaryPoints();                                         //���ұ߽�
	//void CoordinateTransformation();                                  //�ռ�����ת��
	//void ChangePoints(point3d & p1, point3d & p2);                     //��������
	//double Angle(point3d pt1, point3d pt2, point3d pt3);                //���������н�
	//
	//void SurfaceFitting();                                            //ƽ�����
	//
	//void CIntersecPoints();                                           //ֱ����ƽ���ཻ�����
	//void GetAllCell(int Unum, int Vnum, int Wnum, int Uorder, int Vorder, int Worder, varray<double> &uknots, varray<double> &vknots, varray<double> &wknots, varray<point3d> &bpts);//����Ȳ��ߵ�Ԫ
	//point3d GetCubicPoint(int Unum, int Vnum, int Wnum, int uid, int vid, int wid, int uorder, int vorder, int worder, varray<double> Nu, varray<double> Nv, varray<double> Nw, varray<point3d> &bpt, double u, double v, double w);//�������ϵĵ�
	//int Findnum(int left, int right, double x, varray<double> kont);
private:
	//��Ƭ�������ģ�ʹ洢
	varray<double>  u, v, w;                           //u,v,w������Ľڵ�ʸ��
	int u_num, v_num, w_num;                           //u,v,w������Ľڵ����
	int u_order, v_order, w_order;                     //u,v,w��������Ľ���
	int Patchnum;                                    //�������ģ�͵�Ƭ��
	varray<point4d> m_allpoints;                     //��Ƭ���Ƶ�����
	varray<varray<point3d>> PatchControlPoints;      //���п��Ƶ�����

	varray<varray<varray<point3d>>> BoundaryCurve;   //�߽��

	//����ƽ����Ƭ
	float _A, _B, _C, _D,_begin,_end,_each;                                  //��������ƽ��
	int _degree;   //�ܹ��з��˶��ٸ�ƽ��
	vector<vector<float>> eachSurfacePoints;  //ÿһ����ƽ������е�
	vector<vector<point3d>> eachSurfacePoints_point3d;

	//�ռ�����ת��
	double sinrx, cosrx, sinry, cosry;                  //ת���н�

	//��Ƭ����洢
	varray<point3d> m_allIntersetPoints;
	varray<varray<point3d>> m_allSupportPoints;
	varray<point3d> SupportPointsCube;
	varray<varray<double>> Tran;

	//���������Ϣ
	double Wdegree, Vdegree, Udegree;//ϸ�־���
	varray<point3d> m_allcubepoints;
	int Cubenum;
};

