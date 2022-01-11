#pragma once

#include <qstring.h>
#include <fstream>
#include <QOpenGLFunctions_4_5_Core>
#include "stdafx.h"
#include "RWGeometric.h"
#include "PolyIGA.h"

class Option:protected QOpenGLFunctions_4_5_Core
{
public:
	//��һ�λ�ͼ���ǵڶ��λ�ͼ��0�ǵ�һ�Σ�����ڶ���
	int m_again = 0;
	Option();
	~Option();
	//��������������
	void setLockOn();//����
	void setLockClose();//����
	bool isLcokOn();//�ж���״̬

	//�ļ���ȡ֮��Ĳ���
	void OnBnClickedCreateModel(QOpenGLContext * p, QWidget * father);//������ʾNurbsVol
	void OnBnClickedCreateModel(QOpenGLContext * p, QWidget * father, CPolyParaVolume polyv);//��ȡһ��CPolyParaVolume����ʾ
	void OnBnClickedViewre();
	void OnBnClicked3MF();//����3mf�ļ�

	//ģ�͵���Ⱦ����
	void createData(QOpenGLContext * p, varray<varray<varray<point3d>>> varrayPoint, varray<varray<varray<point3d>>> varrayLine);//��Ⱦ����ģ��
	void createData(QOpenGLContext * p, varray<varray<varray<Vec3>>> varrayPoint, varray<varray<varray<Vec3>>> varrayLine);//��Ⱦ����ģ��
	void createAllData(QOpenGLContext * p, QWidget * father, CPolyParaVolume polyv);//�Ǿ���ģ�͵���Ⱦ
	void createData(QOpenGLContext * p, vector<float> data);//��Ⱦһ���Ǿ���ģ��Ƭ
	void createMode();
	void createCP(QOpenGLContext * p);//���Ƶ������

	//ת������
	void TransCnurbsvolToPolyIGA(varray<NurbsVol> &nurbsVols, const CPolyParaVolume poly);//Curbsת��ΪCPolyParaVolume

	//����֮�е���ʹ�õ���Ⱦ����
	void drawNurbsVol(QOpenGLContext * p, QWidget * father, varray<NurbsVol> nvs);//��ʾ����NurbsVol
	void drawNurbsSurface(QOpenGLContext * p, QWidget * father, varray<NurbsSurface> nsfs);//��ʾ����NurbsSurface
	void drawNurbsLine(QOpenGLContext * p, QWidget * father, varray<NurbsLine> nls);//��ʾ����NurbsLine
	void drawBezierLine(QOpenGLContext * p, QWidget * father, varray<BezierLine> bzls);//��ʾBezierLine
	void drawpoint4d(QOpenGLContext * p, QWidget * father, varray<varray<point4d>> p4d2);//��ʾpoint4d
	void drawpoint3d(QOpenGLContext * p, QWidget * father, varray<varray<point3d>> p3d2);//��ʾpoint3d
	void reDraw(QOpenGLContext * p, QWidget * father);//���»���

	void drawSplineVolume(QOpenGLContext * p, QWidget * father, varray<SplineVolume> nvs);//��ʾ����SplineVolume
	void drawSplineSurface(QOpenGLContext * p, QWidget * father, varray<SplineSurface> nsfs);//��ʾ����SplineSurface
	void drawSpline(QOpenGLContext * p, QWidget * father, varray<Spline> nls);//��ʾ����Spline
	void drawVec4(QOpenGLContext * p, QWidget * father, varray<varray<Vec4>> p4d2);//��ʾVec4
	void drawVec3(QOpenGLContext * p, QWidget * father, varray<varray<Vec3>> p3d2);//��ʾVec3
private:
	//��ǰ����ͼ�δ��ڵľ��
	QOpenGLContext * _openglContext;
	//�ؼ���ֵ����
	QString m_modelPath;
	//���ݸ��˵ĵ��Կ������ò�ͬ������̸߳�����
	//ע�⣬�������Ǽ����ܼ��ͳ��򣬲�����Խ����߳�Խ��
	//�߳����෴������ɳ���ļ��������½�
	int m_threadNum = 6;
	//�����޸�������ͼ����ʾ֮�����
	bool _modifyLock = false;
protected:
	////���������Setting
	//void SaveSetting();
	////ˢ�²���
	//void ReflashPara();
	//���㷨����
	point3d LegalVector(point3d a, point3d b, point3d c);
	//string�ָ�
	void Split(const std::string& s, varray<std::string>& v, const std::string& c);
};

