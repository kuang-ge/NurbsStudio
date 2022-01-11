#pragma once
#include "varray.h"
#include "pointXd.h"
#include "CToric.h"
#include "CNurbs.h"

class toricSegment
{
private:
	//���
	point3d m_CenPts;
	varray<point3d> m_linesDu;//�Ǽ����׵���ʸ
	varray<NurbsLine> m_lines;//�Ǽ���
	double m_r;

	//������
	varray<point3d> m_CenLine;
public:
	//���ʸ��Ϣ
	struct outsideVec
	{
		point3d pts;//�ǵ�����
		point3d vec;//���ʸ����(��λ��)
		point3d twinPts;//���ʸ����������Ľ���
		bool isEdge = false;//�Ƿ�Ϊ�߽��
		bool hasTwinPts = false;
	};

	//������Ϣ
	struct lineinf
	{
		int idx0;//���߶˵��Ӧ��cube
		int idx1;
		varray<point4d> ctrlPts;
	};
	varray<lineinf> m_lineInf;
	

	struct VecTop
	{
		int volIdx = -1;
		varray<int> coners;//[4]���ǵ���
		varray<varray<int>> linePara;//[4][3],ÿ���ߵ�uvw�����2��ʾΪ������������3Ϊ����
	};
	varray<VecTop> m_OutsideVecsTop;//���ʸ����
	varray<outsideVec> m_OutsideVecs;//���ʸ
	varray<NurbsSurface> m_OutsideSurface;
	NurbsLine sweepPath;
	//�ָ��ĹǼ�����
	//varray<varray<NurbsLine>> m_SegLines;

	//varray<varray<point4d>> m_BondaryCtrlPts;//�߽����߿��Ƶ㣬��toric��������˳������
	//varray<double> m_knots;
	varray<varray<point3d>> m_samPts; //�ڲ��ָ���
	//varray<varray<varray<point4d>>> m_segSFBoundary; //�ָ�������߹�����Ƭ�Ŀ��Ƶ�
	//int m_lineNum;

	varray<toric>& m_torics;
	varray<NurbsVol> m_Vols;

	//toricSegment() {};
	toricSegment(varray<toric>& mtoric);
	~toricSegment();

	/*
	*/
	void ToricSegmentation(const point3d CenPts, const varray<NurbsLine>& lines, const double r);

private:

	//ֱ����toric���潻��Ĳ�����ֵ,zֵ��ʾm_toric�еı��
	//p0:ֱ�ߵ�
	//v0:ֱ�߷�������
	double LineCrossTorics(const point3d& p0, const point3d& v0, point3d& res);

	//��һ���ڲ���
	void FirstVols(const double la, const double lb);

	//��ƽ����
	void CalBis(varray<varray<point3d>>& bis);

	//�������ʸ
	void CalOutsideVec();

	//���ʸ����
	void FindOutsideVecTop();

	int FindIndex(const point3d pts);

	//ȥ���غ���
	void DelSurface();

	//������Ԫ���Ƿ���ͬ(����)
	bool IsSameEle(const varray<int>& a, const varray<int>& b);

	//����
	void Sampling(const int topIdx, const int lineIdx, const int num, varray<point3d>& sampts);

	//�߽�������
	void CreateBounVol();

	//������toric���潻��Ĳ�����ֵ
	//p0:ֱ�ߵ�
	//v0:ֱ�߷�������
	//��������������(-1,-1)
	double LineCrossToric(const point3d& p0, const point3d& v0, const int ToricIdx, point2d& res);

	//���߲���
	int FindLineInf(const int idx0, const int idx1, bool& isRev);

	//����top��ϢѰ�����ߵ�
	//topIdx��m_OutsideVecsTop����
	//lineIdx��VecsTop�бߵ�����
	//x:���������ֵ
	point3d FindPrismPts(const int topIdx, const int lineIdx, const double x);

	//�߽������
	void FittingBounSF();

	//�߽��ж�
	int InEdge(toric & tor, const point2d & p0);

	//����VecTop���cube���ϵĿ��Ƶ�
	void GetSFCtrlPtsWithVecTop(const VecTop & vectop, varray<point4d>& ctrlPts);

	//��ţ�ٵ�����
	double NewtonItera(const point3d& p0, const point3d& v0, const int ToricIdx, point2d& res);

	//test
	void WriteData();
	void testData();
};