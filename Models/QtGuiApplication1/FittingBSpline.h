//����B��������������С���˷����

#pragma once
#ifndef CB_spline
#define CB_spline
#include "XVec.h"


//#include "varray.h"
#include "XVec.h"
//#include "D:/study/project/Visual Studio/Library/Math/Eigen/dense"
#include "MathUtil.h"
#include "eigen/Eigen/dense"

#include <string>
#include <cmath>
#include <iostream>

namespace bsl
{
	using namespace base;
	using namespace Eigen;
	class FitBSpline
	{
	public:
		//Ĭ�Ϲ��캯��
		FitBSpline();

		//Ĭ����������
		~FitBSpline();
	public:
		double m_FittingErr;//�������׼�
		varray<Vec4> m_CtrlPts;//���Ƶ㼯
		varray<double> m_Knots;//�ڵ�ʸ��
		int m_Degree;//���ߴ���
		
		varray<double> m_oriPtsPara;//��ɢ������������
	private:
		varray<Vec4> m_oriPtsVec4s;//��ɢ�㼯
		int m_CtrlPtsNum;//���Ƶ�����
		//�����ҳ�����������ɢ��
		void ParaPts();

		//�����ڵ�ʸ��
		void CreateKnots();

		//��С���˷�����Ƶ�
		void CalCtrlPts();

		//������Ͼ���
		void CalFittingErr();

	public:
		//B��������
		Vec4 Bspl(double t);

		/*�������
		inputVecs����ɢ�㼯
		sdegree�����ߴ���
		ctrlNum:����������Ƶ�����
		*/
		void FittingBspl(const varray<Vec4>& inputVecs, int sdegree = 3, int ctrlNum = 4);

		/*�������
		inputVecs����ɢ�㼯
		Knots:�ڵ�ʸ��
		sdegree�����ߴ���
		ctrlNum:����������Ƶ�����
		*/
		void FittingBspl(const varray<Vec4>& inputVecs, varray<double> Knots, int sdegree = 3, int ctrlNum = 4);

		//��PtsVec��ֵ��ʹ������������ptsNum
		void Interpolate(varray<Vec4>& PtsVec, int ptsNum);
	};



//////////////////////////////FitBSplineSurface////////////B-����������ɢ�����///////////////////////////////////

	class FitBSplineSurface
	{
	public:
		varray<varray<Vec4>> m_uvCtrlPts;//���Ƶ㼯
		varray<double> m_uKnots;//u����ڵ�ʸ��
		varray<double> m_vKnots;//v����ڵ�ʸ��
		double m_FittingErr;

		FitBSplineSurface();
		~FitBSplineSurface();

		/*�ܼ�ɢ�ҵ��������
		allInputPts:��ɢ�㼯
		all4Edge:��ɢ��߽�
		degree:���ߴ���
		uCtrlPnum:u������Ƶ���
		vCtrlPnum:v������Ƶ���
		iterations:��������,��ʡ�Դ˲�����ʹ������ɢ�ҵ�������
		*/
		void FittingBsurfaceNew(varray<Vec4>& allInputPts, varray<varray<Vec4>>& all4Edge,
			int udegree, int vdegree, int uCtrlPnum, int vCtrlPnum, int iterations);


		/*������ɢ������������
		allInputdata����ɢ�㼯
		all4edge���߽�㼯
		degree�����ߴ���
		uCtrlNum:u������Ƶ�����
		vCtrlNum:v������Ƶ�����
		*/
		void FittingBsurface(const varray<Vec4>& allInputdata, const varray<varray<Vec4>>& all4edge,
			int udegree = 3, int vdegree = 3, int uCtrlNum = 4, int vCtrlNum = 4);

		//private:
		varray<varray<Vec4>> m_uPtsVec;//�̶�v��u����ĵ㼯

	private:
		struct Bsurface//����ṹ��
		{
			int udegree;//U�������ߴ���
			int vdegree;//V�������ߴ���
			varray<varray<Vec4>> uvCtrlPts;//���Ƶ㼯
			varray<double> uKnots;//u����ڵ�ʸ��
			varray<double> vKnots;//v����ڵ�ʸ��

			void clear()
			{
				uvCtrlPts.clear();
				uKnots.clear();
				vKnots.clear();
			}
		};

		//�����ԭʼ����
		varray<Vec4> m_allInputPts;//��ɢ�㼯
		varray<varray<Vec4>> m_all4Edge;//�߽���ɢ�㼯
		int m_udegree;//U�������ߴ���
		int m_vdegree;//V�������ߴ���
		int m_uCtrlPtsNum;//u������Ƶ�����
		int m_vCtrlPtsNum;//v������Ƶ�����
		
		//���̱���
        Bsurface m_baseSurface;//����

		Vec4 m_uvOriginPts;//u,vԭ�㣨m_all4Edge[0]��m_all4Edge[1]���㣩
		Vec4 m_uDir;//u��������(m_all4Edge[0]Ϊu����m_all4Edge[1]Ϊv����
		Vec4 m_vDir;//v��������(m_all4Edge[0]Ϊu����m_all4Edge[1]Ϊv����

		varray<int> m_uTx; //uvw����ϵ��Ӧ��xyz����ϵ, 0 = x, 1 = y, 2 = z

		varray<varray<Vec4>> m_EdgeCtrlPts;//�߽���Ƶ�
		varray<varray<double>> m_EdgePara;//�߽���ɢ�����
		
		varray<varray<Vec4>> m_uCtrlPts;//���u�����Ŀ��Ƶ㼯
		varray<varray<double>> m_uPara;
		varray<varray<double>> m_vPara;

		//B�������溯��
		Vec4 BSF(Bsurface surface, double u, double v);

		//ת�ö�άPtsVec
		template<typename _T>
		void TransPose(varray<varray<_T>>& PtsVec);

		/*Coons��ֵ
		EndgCtrlPts:�߽���Ƶ㣨e0,e1,e2,e3˳��
		coonsPatchCtrlpts:���ص�������Ƶ�
		*/
		void CoonsInterpolate(const varray<varray<Vec4>>& EndgCtrlPts, varray<varray<Vec4>>& coonsPatchCtrlpts);

		//��ԭ��
		void FindOriginPts();

		//��ת�㼯�ڵ��˳��
		void Reverse(varray<Vec4>& pnts);

		//����߽�
		void AligningEdge();

		//����u,v��������
		void CalUVDir();

		//ȷ��uvw����ϵ��Ӧ��xyz����ϵ,0=x,1=y,2=z
		void UVWtoXYZ();

		//������dir�ķ����С��������dir=0,1,2
		void Sorted(varray<Vec4>& ptsArr, int dir);

		//���u����
		void FittingUdir();

		//���v����
		void FittingVdir();

//////////////////////////////////////////�ܼ�ɢ�ҵ����/////////////////////////////////////

		//�߽�Coons��ֵ����
		void InterpolateBaseSurface();

		//��ȡ�㼯ptsArr���Ե�ptsΪ���ģ���a*b���ο��ڵĵ�
		void GetNearPts(const varray<Vec4>& ptsArr, Vec4 pts, double a, double b, varray<Vec4>& nearPts);

		//�������ռ��Ӧ��
		void Segmentation(const varray<double>& uPara, const varray<double>& vPara, varray<varray<Vec4>>& ptsArr);

		//���hardy'sϵ������C
		//ptsArr:���������ݵ�
		//b:hardy��ֵ���ڳ���
		void CalHardysC(const varray<Vec4>& ptsArr, double b, MatrixXd& C);

		//Hardy's˫���β�ֵ
		Vec4 HardysInterpolate(const Vec4 pts, const varray<Vec4>& ptsArr, const MatrixXd& C, const double b);

		//���㷶Χ
		void CalExtreme(double extrm[4]);

		//��ֵ�����
		void InterpolatGridPts(int uNum, int vNum);


/////////////////////////////����ɢ�ҵ�������///////////////////////////////////////////////

		//v������8����߽���ж�
		bool IsVEdgePts(Vec4 pts, const varray<Vec4>& ptsArr);

		//�ж�pnt�Ƿ��ڵ㼯ptsArr��,��������,-1Ϊ���ڵ㼯��
		int InArr(Vec4 pnt, const varray<Vec4>& ptsArr);

		//����̶�v��u����ĵ㼯������m_uPtsVec
		void CalUPtsVec();

		//v���������
		void VParaPts();

		//����������
		void CalSurFaceFittingErr();

	};

}

#endif

