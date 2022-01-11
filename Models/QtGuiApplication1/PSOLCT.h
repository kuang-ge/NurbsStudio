
#pragma once
#include "PSO.h"
#include "pointXd.h"
#include "CToric.h"
#include "CNurbs.h"

//PSO����ֱ����toric����߽罻��
class psoLCT :public PSO
{
public:
	psoLCT(const toric& mtoric,const int Eidx);

	virtual ~psoLCT() {}

	//PSO����ֱ����toric���潻��Ĳ�����ֵ��������Ӧ��
	void CalPSOLCT(point2d& para,double& dis);

private:
	point2d m_p0;//�߽�ֱ�����
	point2d m_v0;//�߽�ֱ�߷�������
	toric m_toric;
	double m_len;

	//��дĿ�꺯��
	double CalDis(const varray<double>& P);
};

//������������Ľ�ƽ����
class psoCenLine :public PSO
{
public:
	psoCenLine();

	virtual ~psoCenLine(){}

	//������������Ľ�ƽ����
	void CalPsoCenLine(const varray<point3d>& lineDu, point3d& Cenline, double& dis);

private:
	varray<point3d> m_lineDu;

	//��дĿ�꺯��
	double CalDis(const varray<double>& P);
};

class psoPtsDis :public PSO
{
public:
	psoPtsDis();
	virtual ~psoPtsDis() {}

	//���������P0����Ϊdis�ĵ�
	void FindPtsWithDis(const NurbsLine line, const point3d P0, const double dis, double & u,double& err);

private:
	NurbsLine m_line;
	point3d m_P0;
	double m_dis;
	//��дĿ�꺯��
	double CalDis(const varray<double>& P);
};