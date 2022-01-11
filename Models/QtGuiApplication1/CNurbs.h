//Nurbs��صı�׼����
//�㷨���ԡ��Ǿ�������B����(��2��)�� Les Piegl,Wayne Tiller��; �����; �廪��ѧ������
//@ USST ���ı� 2018

#pragma once

#include "stdafx.h"
#include "XVec.h"
class BezierLine
{
	//����u������n��Bernstein����ʽ��ֵ
	void AllBernstein(const int n, const double u, varray<double>& B)const;
public:
	int m_Degree;//���ߴ���
	varray<point4d> m_CtrlPts;//���Ƶ�
	//����Bezier���ߵ�
	point3d GetLinePoint(const double u)const;
	//����Bezier����
	void CalLinePoint(const int Unum, varray<point3d>& linePts)const;
};

	class NurbsBase
	{
	public:
		/*����ڵ��±�
		x���ڵ�ֵ
		degree������
		CtrlPtsNum�����Ƶ�����
		knots���ڵ�ʸ��*/
		int FindSpan(const double x, const int degree, const int CtrlPtsNum, const varray<double>& knots)const;

		/*���ݲ���ֵ���ڵ��±꣬���������
		//u������ֵ
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

	//NURBS����
	struct NurbsLine :public NurbsBase
	{
		int m_Degree;//���ߴ���
		varray<double> m_Knots;//���߽ڵ�ʸ��
		varray<point4d> m_CtrlPts;//���Ƶ�

		//����ֱ����ʽ����
		void CreatLineWithTwoPoints(const point3d &p1, const point3d &p2, int degree = 2);

		/*����u�ڵ��Ӧ�����������
		u���ڵ����*/
		point3d GetLinePoint(const double u)const;

		/*��ȡ���ϵĵ�
		Unum�����ߵ�����
		linePts�����ߵ�	*/
		void CalLinePoint(const int Unum, varray<point3d>& linePts)const;

		/*����ʸֵ����A(u)����p�׵���ʸֵ����A(u)��NURBS���߷���ʽ
		u������
		n����ʸ����
		Der��A(u)����p�׵�*/
		void PtsDerivsAu(const double u, const int n, varray<point4d>& Der)const;

		/*������u������n�׵�
		u������
		n����ʸ����
		Der������n�׵� */
		void PtsDerivs(const double u, const int n, varray<point4d>& Der)const;

		/*��������n�׵�
		step��������ȷֲ�����
		Der�����ߵ��Ӧ��n�׵�*/
		void CurveDerivs(const int n, const double step, varray<varray<point3d>>& Der);

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
		void Segmentation(const double u, varray<NurbsLine>& lines);

		/*NURBS���߷ֽ�ΪBEZIER
		Qw�������BEZIER���߶ο��Ƶ�*/
		void Decompose(varray<varray<point4d>>& Qw);

		//��ת���߷���
		void CruveReverse();
	};

	//NURBS����
	struct NurbsSurface :public NurbsBase
	{
		int m_uDegree;//v�������
		int m_vDegree;//v�������
		int m_uNum;//v������Ƶ�����
		int m_vNum;//v������Ƶ�����
		varray<double> m_uKnots;//v����ڵ�ʸ��
		varray<double> m_vKnots;//v����ڵ�ʸ��
		varray<point4d> m_CtrlPts;//���Ƶ�,��v-u����洢

		void SetSurface(const int uDegree, const int vDegree, const int uNum, const int vNum, 
			const varray<double>& uKnots, const varray<double>& vKnots);



		//���㣨u��v����Ӧ�������ϵĵ�
		point3d GetSurFacePoint(const double u, const double v)const;

		//�����ı�����Ƭ��ʾ����
		//num:�÷������������
		//quads:��Ƭ����
		//lines:�Ȳ�������
		threadParam CalQuads(const int Unum, const int Vnum, varray<varray<point3d>>& quads, varray<varray<point3d>>& lines)const;

		/*Coons��ֵ
		EndgCtrlPts:�߽���Ƶ㣨e0,e1,e2,e3˳��*/
		void CoonsInterpolate(const varray<varray<point4d>>& EdgeCtrlPts);

		/*Coons��ֵ
		EdgeLines:�߽����ߣ�e0,e1,e2,e3˳��*/
		void CoonsInterpolate(const varray<NurbsLine>& EdgeLines);

        //��������
		//Udegree,Vdegree:���׺����
		void DegreeElevate(const int Udegree, const int Vdegree);

		//����ڵ����
		//Uknot,Vknot:��Ҫ����Ľڵ�
		void KnotsRefine(const varray<double>& Uknot, const varray<double>& Vknot);

		//����ֶ�ΪBezier����
		//dir:0=U����1=V����
		//QW��Bezier������Ƶ�
		void Decompose(const bool dir, varray<varray<varray<point4d>>>& QW);  

		//���ݿ��Ƶ��ά��ż���һά���
		int CtrlPtsIdx(const int uIdx, const int vIdx);

		//���Ƶ�����ת��ΪU-V
		void OrderCtrlPts();

		//���Ƶ�����ת��ΪU-V
		void OrderCtrlPts(NurbsSurface& sf);

		//��ȡ�߽���,�����Ѿ���U-V�����
		void GetEdgeLines(varray<NurbsLine>& EdgeLines);

	};

	//NURBS��
	struct NurbsVol :public NurbsBase
	{
		int m_uDegree;//u�������
		int m_vDegree;//v�������
		int m_wDegree;//w�������
		int m_uNum;//v������Ƶ�����
		int m_vNum;//v������Ƶ�����
		int m_wNum;//w������Ƶ�����
		varray<double> m_uKnots;//v����ڵ�ʸ��
		varray<double> m_vKnots;//v����ڵ�ʸ��
		varray<double> m_wKnots;//w����ڵ�ʸ��
		varray<point4d> m_CtrlPts;//���Ƶ㣬��v-u-w����洢

		void SetVol(const int uDegree, const int vDegree, const int wDegree, const int uNum, const int vNum, const int wNum,
			const varray<double>& uKnots, const varray<double>& vKnots, const varray<double>& wKnots);

		//(u,v,w)��Ӧ�����ϵĵ�
		point3d GetVolPoint(const double u, const double v, const double w)const;
		

		//�����ı�����Ƭ��ʾ����
		//num:�÷������������
		//quads:��Ƭ����
		//lines:�Ȳ�������
		threadParamVOL CalQuads(const int Unum, const int Vnum, const int Wnum, 
			varray<varray<varray<point3d>>>& quads, varray<varray<varray<point3d>>>& lines)const;

		//���ݿ��Ƶ���ά��ż���һά���
		int CtrlPtsIdx(const int uIdx, const int vIdx, const int wIdx);

		//���Ƶ�����ת��ΪU-V-W
		void OrderCtrlPts();

		//���Ƶ�����ת��ΪU-V-W
		void OrderCtrlPts(NurbsVol& vol);

		//�ڵ����
		//Uknot,Vknot,Wknot:��Ҫ����Ľڵ�
		void KnotsRefine(varray<double>& Uknot, varray<double>& Vknot, varray<double>& Wknot);
		
		//����
		//Udegree,Vdegree,Wdegree:���׺����
		void DegreeElevate(const int Udegree, const int Vdegree, const int Wdegree);

		/*ɨ������Nurbs��ģ�ͣ����洹ֱ��·��
		pathT��ɨ��·��
		nurbsSF����ʼ����
		K������ʵ������,һ��ȡ·�����Ƶ�������1
		*/
		void CreateSweepNurbsVol(const NurbsLine& pathT, const NurbsSurface& nurbsSF, const int K);

		/*ƽ��ɨ������Nurbs��ģ�ͣ����治��ֱ��·��
		pathT��ɨ��·��
		nurbsSF����ʼ����*/
		void CreateTransSweepNurbsVol(const NurbsLine& pathT, const NurbsSurface& nurbsSF);

		//����
		//path:·��
		//surfaces:��������
		void LoftingNurbsVol(const NurbsLine& path, const varray<NurbsSurface>& surfaces);

		void LoftingNurbsVol2(const NurbsLine& path, const varray<NurbsSurface>& surfaces);

		//����
		//path:·��
		//surfaces0:·������������
		//surfaces1:·���յ��������
		void LoftingNurbsVol(const NurbsLine& path, const NurbsSurface& surfaces0,const NurbsSurface& surfaces1);

		//��ȡ�����߽���
		//dir:1->6��ʾ����
		NurbsSurface GetSingleSurfaces(int dir) const;

		//��COONS��ֵ
		void VolCoonsInterpolate(const varray<NurbsLine>& EdgeLines);

		//refinenum:ϸ������
		void Knots_Refine_Num(int refinenum)
		{
			//����²����Ľڵ㣬��Ϊϸ�������Ĳ���
			varray <double> new_uknots;
			varray <double> new_vknots;
			varray <double> new_wknots;
			int i = 1;
			while (i <= refinenum)
			{
				//�ҳ��м�ڵ� u����
				for (int i = 1; i < m_uKnots.size(); i++)
				{
					if (m_uKnots[i - 1] != m_uKnots[i])
					{
						new_uknots.push_back((m_uKnots[i - 1] + m_uKnots[i]) / 2);
					}
				}

				//�ҳ��м�ڵ� v����
				for (int i = 1; i < m_vKnots.size(); i++)
				{
					if (m_vKnots[i - 1] != m_vKnots[i])
					{
						new_vknots.push_back((m_vKnots[i - 1] + m_vKnots[i]) / 2);
					}
				}
				//�ҳ��м�ڵ� w����
				for (int i = 1; i < m_wKnots.size(); i++)
				{
					if (m_wKnots[i - 1] != m_wKnots[i])
					{
						new_wknots.push_back((m_wKnots[i - 1] + m_wKnots[i]) / 2);
					}
				}
				//ϸ��
				KnotsRefine(new_uknots, new_vknots, new_wknots);

				new_uknots.clear();
				new_vknots.clear();
				new_wknots.clear();

				i++;
			}
		}
	private:
		//����Ȳ���
		//uvw:0=u,1=v,2=w
		//t:����
		//num:��u-v-w˳��
		//L:�Ȳ����ı�����Ƭ��
		void CalIsoSurface(const int uvw, const double t, const int num1, const int num2,
			varray<varray<point3d>>& quads, varray<varray<point3d>>& lines)const;

		/*ɨ��ʱ�Ľ���ʵ��λ��
		pathT��ɨ��·��
		K������ʵ������
		pos������ʵ��λ��
		NewKnots���½ڵ�ʸ��*/
		void InsLocation(const NurbsLine & pathT, int K, varray<double> & pos);

		/*����ɨ��·���ϵľֲ�����ϵ
		pathT��ɨ��·��
		pos������ʵ��λ��
		TranMat��pos���ľֲ�����ϵ*/
		void LocalCoordinates(const NurbsLine & pathT, const varray<double> & pos,
			varray<varray<point4d>> & TranMat);

		/*��ɨ��·���Խ�����б任
		nurbsSF��������Ƶ�
		TranMat���任����
		OriCoordinates:����ľֲ�����ϵ
		allNurbsSF���õ������н���ʵ�����Ƶ�
		*/
		void MatrixTran(const varray<varray<point4d>> & nurbsSF, const varray<varray<point4d>>& TranMat,
			varray<varray<varray<point4d>>>& allNurbsSF);

		//ɨ������Nurbs��
		void SweepSurface(const varray<varray<varray<point4d>>>& allNurbsSF, 
			const varray<varray<varray<double>>>& SFw, const varray<double>& pos);

		//ȡ��ߴ�
		void MaxDegree(const varray<NurbsSurface>& surfaces, int& uDegree, int& vDegree);

		//�ڵ�ʸ������
		void KnotsUnify(const varray<NurbsSurface>& surfaces, varray<double>& NewUKnots, varray<double>& NewVKnots);
	};