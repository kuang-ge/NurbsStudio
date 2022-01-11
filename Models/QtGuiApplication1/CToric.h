#pragma once
//toric����
//USST ���ı� 2018/11

#include "varray.h"
#include "pointXd.h"

class toric
{
public:
	varray<point2d> m_ConerPts;//��ά������Ķ���,���ж�������Ӧ���ڵ���0,����ʱ������
	varray<varray<double>> m_cm;//����Ӧ�Ļ�����ϵ��������0�������ҪΪ���������棬��ȡ����ʽϵ��
	varray<varray<point4d>> m_CtrlPts;//���Ƶ�

	//���³�Ա�ɹ��캯�������
	varray<varray<point2d>> m_LatticePts;//��ά��������
	int m_LineN;//������ߵ�����
	varray<varray<int>> m_BoudaryLine;//������߽�ֱ�߷��̲���h(u,v)=au+bv+c=0
	int m_xmin, m_xmax, m_ymin, m_ymax;//������ֵ

public:
	//Ĭ�Ϲ��캯���������ڷ���ռ�
	toric();

	//���캯��
	toric(const varray<point2d>& ConerPts);

	//���캯��
	//ConerPts����ά������Ķ��㣬���ж�������Ӧ���ڵ���0,����ʱ������
	//CtrlPts�����Ƶ�
	//cm������Ӧ�Ļ�����ϵ���������ҪΪ���������棬��ȡ����ʽϵ��
	toric(const varray<point2d>& ConerPts, const varray<varray<point4d>>& CtrlPts, const varray<varray<double>>& cm);

	//��������
	~toric();

	//����Toric�����
	point3d GetToricPts(const double u, const double v)const;

	//����Ȳ���
	void CalSurfacePts(const int UptsNum, const int VptsNum, varray<varray<point3d>>& uPts, varray<varray<point3d>>& vPts);

	//����߽�����
	void CalBondaryPts(const int ptsNum, varray<varray<point3d>>& BondaryPts);

	//�ж��Ƿ�Ϊtoric�߽���, ���򷵻ر߽�ţ����򷵻�-1
	int InEdge(const point2d &p);

	//�жϲ����Ƿ�Ϸ�
	bool IsInSF(const double u, const double v);

	//COONS��ֵ���۽��������ֵ
	//BondaryCtrlPts���߽����߿��Ƶ㣬��toric��������˳������
	void CoonsToric(const varray<varray<point4d>>& BondaryCtrlPts);

protected:
	varray<varray<point4d>> m_BondaryCtrlPts;//�߽����߿��Ƶ㣬��toric��������˳������

	//Coonsֱ����,���������Ƶ�
	void CoonsS12TCtrlPts(varray<varray<point4d>>& UCtrlPts, varray<varray<point4d>>& VCtrlPts, varray<varray<point4d>>& TCtrlPts);

	//�����׵㷨ʸ
	//center���������
	//point3d FirstDir(const varray<point4d>& cur, const point3d center);

	//Ѱ�ҵ�ֵ����߽罻��
	//val:��ֵ�ߵ�ֵ
	//UorV:trueΪU��ֵ�ߣ�falseΪV��ֵ��
	//startVal����ֵ���������߽�Ľ���
	//endVal����ֵ�������һ���߽�Ľ���
	//void FindStartAndEndVal(const double val, const bool UorV, double& startVal, double& endVal);


	//���ݸ���������������߽緽��
	void CalBoudaryLine();

	//�����������
	void CalParaLattice();

	//����(u,v)��Ӧ�Ķ�Ԫ�߽緽��ֵ
	double CalhiVal(const double u, const double v, const int lineIdx)const;
	double CalhiVal(const point2d& m, const int lineIdx)const;

	//����(u,v)��Ӧ�Ļ�����ֵ
	/*m����������
	cm��������ϵ��
	*/
	double ToricBasic(const double u, const double v, const point2d& m, const double cm)const;

	//������߽緽�̷Ǹ���У��
	void LineNonNeg(varray<int>& line, int idx);

	//(u,v)�����m��Ӧ�Ļ�����һ�׵���
	//dir��true��ʾ��U�󵼣�false��ʾ��V��
	double CalhiDer1(const double u, const double v, const point2d& m, const double cm, const bool dir);

	//(u,v)�����m��Ӧ�Ļ��������׵���
	//dir����0��1��2��3��=��uu,vv,uv,vu��
	double CalhiDer2(const double u, const double v, const point2d& m, const double cm, const int dir);

public:
	//������U��ֵ�ߵĶ˵�,���ض˵���m_LatticePts������
	void FindEndpointV(const point2d& pts, point2d& stIdx, point2d& edIdx);

	//toric����ʽ
	point3d CalM(const double u, const double v);

	//toric����ʽƫ��
	//dir��true��ʾ��U�󵼣�false��ʾ��V��
	point3d CalMDir1(const double u, const double v, const bool dir);

	//(u,v)��һ�׵���
	//dir��true��ʾ��U�󵼣�false��ʾ��V��
	point3d CalPtsDir1(const double u, const double v, const bool dir);

	//(u,v)�����׵���
	//dir����0��1��2��3��=��uu,vv,uv,vu��
	point3d CalPtsDir2(const double u, const double v, const int dir);
};
