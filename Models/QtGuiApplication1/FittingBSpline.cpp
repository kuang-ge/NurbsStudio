
#include "FittingBSpline.h"
#include "qmatrix.h"

using namespace std;
namespace bsl
{

	/////////////////FitBSpline///////////////////////////////FitBSpline///////////////////


	//Ĭ�Ϲ��캯��
	FitBSpline::FitBSpline()
	{
		m_FittingErr = 10000.0;
	}

	//Ĭ����������
	FitBSpline::~FitBSpline()
	{
	}

	//�����ҳ�����������ɢ��
	void FitBSpline::ParaPts()
	{
		double ltemp, sumL = 0.0;
		//�ۼ��ҳ�
		for (int i = 1; i < m_oriPtsVec4s.size(); ++i)
		{
			ltemp = (m_oriPtsVec4s[i] - m_oriPtsVec4s[i - 1]).Magnitude();
			sumL += ltemp;
		}
		//�淶���ҳ��������ݵ�Ĳ�����
		m_oriPtsPara.clear();
		m_oriPtsPara.push_back(0.0);
		double utemp = 0.0;
		for (int i = 1; i < m_oriPtsVec4s.size() - 1; ++i)
		{
			utemp = m_oriPtsPara[i - 1] + (m_oriPtsVec4s[i] - m_oriPtsVec4s[i - 1]).Magnitude() / sumL;
			m_oriPtsPara.push_back(utemp);
		}
		m_oriPtsPara.push_back(1.0);
	}

	//�����ڵ�ʸ��
	void FitBSpline::CreateKnots()
	{
		m_Knots.clear();
		for (int i = 0; i < m_CtrlPtsNum + m_Degree + 1; i++)
		{
			double temp = 0.0;
			if (i > m_Degree && i < m_CtrlPtsNum)
			{
				double c = 1.0* m_oriPtsVec4s.size() / (m_CtrlPtsNum - m_Degree);
				int j = i - m_Degree;
				int k = (int)(j*c);
				double a = j*c - k;
				temp = (1 - a)*m_oriPtsPara[k - 1] + a*m_oriPtsPara[k];
			}
			else if (i > m_CtrlPtsNum - 1)
				temp = 1.0;
			else;
			m_Knots.push_back(temp);
		}
	}

	//��С���˷���ǰ���Ƶ������µĿ��Ƶ�
	void FitBSpline::CalCtrlPts()
	{
		varray<Vec4> r;
		Vec4 temp = Vec4(0.0, 0.0, 0.0);//��ʱ������������㼯����Ԫ��
		r.push_back(temp);
		int oriptNum = m_oriPtsVec4s.size();
		for (int i = 1; i < oriptNum - 1; ++i)
		{
			temp = m_oriPtsVec4s[i]
				- m_oriPtsVec4s[0] * OneBasisFun(m_Degree, m_Knots.size() - 1, m_Knots, 0, m_oriPtsPara[i])
				- m_oriPtsVec4s[oriptNum - 1] * OneBasisFun(m_Degree, m_Knots.size() - 1, m_Knots, m_CtrlPtsNum - 1, m_oriPtsPara[i]);
			r.push_back(temp);
		}
		//������С���˷�����������
		MatrixXd matN(oriptNum - 2, m_CtrlPtsNum - 2), matNT(m_CtrlPtsNum - 2, oriptNum - 2), matR(m_CtrlPtsNum - 2, 3), matD(m_CtrlPtsNum - 2, 3);
		//matN
		for (int i = 1; i < oriptNum - 1; ++i)
		{
			for (int j = 1; j < m_CtrlPtsNum - 1; ++j)
				matN(i - 1, j - 1) = OneBasisFun(m_Degree, m_Knots.size() - 1, m_Knots, j, m_oriPtsPara[i]);
		}
		//matNT
		matNT = matN.transpose();
		//matR
		for (int i = 1; i < m_CtrlPtsNum - 1; ++i)
		{
			temp = Vec4(0.0, 0.0, 0.0);
			for (int j = 1; j < m_oriPtsVec4s.size() - 1; ++j)
			{
				temp = temp + r[j] * OneBasisFun(m_Degree, m_Knots.size() - 1, m_Knots, i, m_oriPtsPara[j]);
			}
			matR(i - 1, 0) = temp.x;
			matR(i - 1, 1) = temp.y;
			matR(i - 1, 2) = temp.z;
		}
		//NT*N*D=R
		MatrixXd A = matNT*matN;
		matD = A.colPivHouseholderQr().solve(matR);
		//���������Ƶ㼯
		m_CtrlPts.clear();
		m_CtrlPts.push_back(m_oriPtsVec4s[0]);
		for (int i = 0; i < m_CtrlPtsNum - 2; ++i)
		{
			temp = Vec4(matD(i, 0), matD(i, 1), matD(i, 2));
			m_CtrlPts.push_back(temp);
		}
		m_CtrlPts.push_back(m_oriPtsVec4s[oriptNum - 1]);
	}

	//B��������
	Vec4 FitBSpline::Bspl(double t)
	{
		Vec4 res(0.0, 0.0, 0.0);
		int crlpnt_size = m_CtrlPts.size();
		for (int i = 0; i < crlpnt_size; ++i)
		{
			res = res + m_CtrlPts[i] * OneBasisFun(m_Degree, m_Knots.size() - 1, m_Knots, i, t);
		}
		return res;
	}

	//������Ͼ���
	void FitBSpline::CalFittingErr()
	{
		double sumeps = 0.0;
		for (int i = 0; i < m_oriPtsVec4s.size(); ++i)
		{
			sumeps += pow((m_oriPtsVec4s[i] - Bspl(m_oriPtsPara[i])).Magnitude(), 2);
		}
		m_FittingErr = sqrt(sumeps / m_oriPtsVec4s.size());
	}

	/*�������
	inputVecs����ɢ�㼯
	sdegree�����ߴ���
	ctrlNum:����������Ƶ�����
	*/
	void FitBSpline::FittingBspl(const varray<Vec4>& inputVecs, int sdegree/* = 3 */, int ctrlNum/* = 4 */)
	{

		using namespace std;
		if (inputVecs.empty())
		{
			cout << "������ɢ�㼯Ϊ�գ���0�˳�" << endl;
			while (cin)
			{
				exit(0);
			}
		}
		else
		{
			m_oriPtsVec4s.clear();
			m_oriPtsVec4s = inputVecs;
		}
		m_Degree = sdegree;
		m_CtrlPtsNum = ctrlNum;

		ParaPts();

		if (m_CtrlPtsNum < m_Degree + 1)
		{
			cout << "���Ƶ���������ȷ����0�˳�" << endl;
			int c_in;
			while (cin >> c_in)
			{
				if (c_in == 0)
					exit(0);
			}
		}
		else
		{
			CreateKnots();
			CalCtrlPts();
			CalFittingErr();
		}
	}

	/*�������
	inputVecs����ɢ�㼯
	Knots:�ڵ�ʸ��
	sdegree�����ߴ���
	ctrlNum:����������Ƶ�����
	*/
	void FitBSpline::FittingBspl(const varray<Vec4>& inputVecs, varray<double> Knots, int sdegree /*= 3*/, int ctrlNum /*= 4*/)
	{
		using namespace std;
		if (inputVecs.empty())
		{
			cout << "������ɢ�㼯Ϊ�գ���0�˳�" << endl;
			while (cin)
			{
				exit(0);
			}
		}
		else
		{
			m_oriPtsVec4s.clear();
			m_oriPtsVec4s = inputVecs;
		}
		m_Degree = sdegree;
		m_CtrlPtsNum = ctrlNum;

		ParaPts();

		if (m_CtrlPtsNum < m_Degree + 1)
		{
			cout << "���Ƶ���������ȷ����0�˳�" << endl;
			int c_in;
			while (cin >> c_in)
			{
				if (c_in == 0)
					exit(0);
			}
		}
		else
		{
			m_Knots = Knots;
			CalCtrlPts();
			CalFittingErr();
		}

	}

	//��PtsVec��ֵ��ʹ������������ptsNum
	//void FitBSpline::Interpolate(varray<Vec4>& PtsVec, int ptsNum)
	//{
	//	int i = 0;
	//	while (PtsVec.size() < ptsNum)
	//	{
	//		//��ֵ��
	//		Vec4 dp = PtsVec[i] + (PtsVec[i + 1] - PtsVec[i]) / 2;
	//		//�����ֵ��
	//		varray<Vec4>::iterator inp = PtsVec.begin() + 2 * i + 1;
	//		PtsVec.push_back(dp);
	//		++i;
	//	}
	//}
	//���º���Ϊ����ģ���޸ģ���������ģ��û���Թ���
	//��PtsVec��ֵ��ʹ������������ptsNum
	void FitBSpline::Interpolate(varray<Vec4>& PtsVec, int ptsNum)
	{
		while (PtsVec.size() < ptsNum)
		{
			varray<Vec4> newPtsVec;
			int num = PtsVec.size() - 1;
			for (int i = 0; i < num; ++i)
			{
				newPtsVec.push_back(PtsVec[i]);
				//��ֵ��
				Vec4 dp = (PtsVec[i + 1] + PtsVec[i]) / 2;
				newPtsVec.push_back(dp);
			}
			newPtsVec.push_back(PtsVec[num]);
			PtsVec = newPtsVec;
		}
	}
//////////////////////////////FitBSplineSurface////////////B-����������ɢ�����///////////////////////////////////

	FitBSplineSurface::FitBSplineSurface()
	{
		m_FittingErr = -1;
	}

	FitBSplineSurface::~FitBSplineSurface()
	{
	}

	/*�ܼ�ɢ�ҵ��������
	allInputPts:��ɢ�㼯
	all4Edge:��ɢ��߽�
	degree:���ߴ���
	uCtrlPnum:u������Ƶ���
	vCtrlPnum:v������Ƶ���
	iterations:��������,��ʡ�Դ˲�����ʹ������ɢ�ҵ�������
	*/
	void FitBSplineSurface::FittingBsurfaceNew(varray<Vec4>& allInputPts, varray<varray<Vec4>>& all4Edge,
		int udegree, int vdegree, int uCtrlPnum, int vCtrlPnum, int iterations)
	{
		//����
		m_allInputPts = allInputPts;
		m_all4Edge = all4Edge;
		m_udegree = udegree;
		m_vdegree = vdegree;
		m_uCtrlPtsNum = uCtrlPnum;
		m_vCtrlPtsNum = vCtrlPnum;
		
		//����߽磬��ԭ�߽����˹�ϵ
		AligningEdge();
		CalUVDir();
		UVWtoXYZ();
		//�߽�Coons��ֵ����
		InterpolateBaseSurface();
		int uNum = m_uCtrlPtsNum;
		int vNum = m_vCtrlPtsNum;

		for (int i = 0; i < iterations; ++i)
		{
			uNum = std::pow(2.0, i)*uNum;
			vNum = std::pow(2.0, i)*vNum;
			bool isuover = uNum > m_all4Edge[0].size();
			bool isvover = vNum > m_all4Edge[1].size();
			if (isuover || isvover)
				break;
			else
			{
				//Hardy��ֵ
				InterpolatGridPts(uNum, vNum);
				//���U
				FittingUdir();
				//���V
				FittingVdir();
				//�������
				m_baseSurface.uvCtrlPts = m_uvCtrlPts;
			}
		}
	}

	//B�������溯��
	Vec4 FitBSplineSurface::BSF(Bsurface surface, double u, double v)
	{
		Vec4 res(0.0, 0.0, 0.0);
		for (int i = 0; i < surface.uvCtrlPts.size(); ++i)
		{
			for (int j = 0; j < surface.uvCtrlPts[i].size(); j++)
				res = res + surface.uvCtrlPts[i][j] * OneBasisFun(surface.vdegree, surface.vKnots.size() - 1, surface.vKnots, j, v)
				* OneBasisFun(surface.udegree, surface.uKnots.size() - 1, surface.uKnots, i, u);
		}
		return res;
	}

	//ת�ö�άPtsVec
	template<typename _T>
	void FitBSplineSurface::TransPose(varray<varray<_T>>& PtsVec)
	{
		varray<varray<_T>> PtsVec_T;
		varray<_T> PtsVec_Ti;
		int PtsVec_sz = PtsVec.size();
		int PtsVec0_sz = PtsVec[0].size();
		for (int i = 0; i < PtsVec0_sz; ++i)
		{
			PtsVec_Ti.clear();
			for (int j = 0; j < PtsVec_sz; ++j)
				PtsVec_Ti.push_back(PtsVec[j][i]);
			PtsVec_T.push_back(PtsVec_Ti);
		}
		PtsVec.clear();
		PtsVec = PtsVec_T;
	}

	/*Coons��ֵ
	EndgCtrlPts:�߽���Ƶ㣨e0,e1,e2,e3˳��
	coonsPatchCtrlpts:���ص�������Ƶ�
	*/
	void FitBSplineSurface::CoonsInterpolate(const varray<varray<Vec4>>& EdgeCtrlPts, varray<varray<Vec4>>& coonsPatchCtrlpts)
	{
		coonsPatchCtrlpts.clear();
		int upnum = EdgeCtrlPts[0].size();//u��������
		int vpnum = EdgeCtrlPts[1].size();//v��������
		Vec4 P00, P10, P01, P11;//�ǵ�
		P00 = EdgeCtrlPts[0][0];
		P10 = EdgeCtrlPts[0][upnum - 1];
		P01 = EdgeCtrlPts[2][0];
		P11 = EdgeCtrlPts[2][upnum - 1];
		for (int i = 0; i < upnum; ++i)
		{
			double tui = 1.0*i / (upnum-1);
			varray<Vec4> Ctrlpts_i;
			for (int j = 0; j < vpnum; ++j)
			{
				double tvj = 1.0*j / (vpnum-1);
				Vec4 Pij = (1 - tui)*EdgeCtrlPts[1][j] + tui * EdgeCtrlPts[3][j]
					+ (1 - tvj)*EdgeCtrlPts[0][i] + tvj * EdgeCtrlPts[2][i]
					- (1 - tvj)*((1 - tui)*P00 + tui * P10)
					- tvj * ((1 - tui)*P01 + tui * P11);
				Ctrlpts_i.push_back(Pij);
			}
			coonsPatchCtrlpts.push_back(Ctrlpts_i);
		}
	}

	//��ԭ��
	void FitBSplineSurface::FindOriginPts()
	{
		Vec4 de[4];
		int e0_s = m_all4Edge[0].size();
		int e1_s = m_all4Edge[1].size();
		Vec4 e0[4] = { m_all4Edge[0][0],m_all4Edge[0][0] ,m_all4Edge[0][e0_s - 1] ,m_all4Edge[0][e0_s - 1] };
		Vec4 e1[4] = { m_all4Edge[1][0],m_all4Edge[1][e1_s - 1] ,m_all4Edge[1][0] ,m_all4Edge[1][e1_s - 1] };
		for (int i = 0; i < 4; ++i)
		{
			de[i] = e0[i] - e1[i];
			if (de[i].x == 0.0 && de[i].y == 0.0 && de[i].z == 0.0)
			{
				m_uvOriginPts = e0[i];
				break;
			}
		}
	}

	//��ת�㼯�ڵ��˳��
	void FitBSplineSurface::Reverse(varray<Vec4>& pnts)
	{
		int p_s = pnts.size();
		Vec4 temp;
		for (int i = 0; i < p_s / 2; ++i)
		{
			temp = pnts[p_s - 1 - i];
			pnts[p_s - 1 - i] = pnts.at(i);
			pnts.at(i) = temp;
		}
	}

	//����߽�
	void FitBSplineSurface::AligningEdge()
	{
		FindOriginPts();
		//��ԭ�����e0��ʹԭ������ǰ
		int e_s = m_all4Edge[0].size();
		Vec4 dtemp = m_all4Edge[0][e_s - 1] - m_uvOriginPts;
		if (dtemp.x == 0.0 && dtemp.y == 0.0 && dtemp.z == 0.0)
			Reverse(m_all4Edge[0]);
		//��ԭ�����e1��ʹԭ������ǰ
		e_s = m_all4Edge[1].size();
		dtemp = m_all4Edge[1][e_s - 1] - m_uvOriginPts;
		if (dtemp.x == 0.0 && dtemp.y == 0.0 && dtemp.z == 0.0)
			Reverse(m_all4Edge[1]);
		//��e1���һ�����e2
		Vec4 elp = m_all4Edge[1][e_s - 1];
		e_s = m_all4Edge[2].size();
		dtemp = m_all4Edge[2][e_s - 1] - elp;
		if (dtemp.x == 0.0 && dtemp.y == 0.0 && dtemp.z == 0.0)
			Reverse(m_all4Edge[2]);
		//��e2���һ�����e3(�������)
		elp = m_all4Edge[2][e_s - 1];
		e_s = m_all4Edge[3].size();
		dtemp = m_all4Edge[3][0] - elp;
		if (dtemp.x == 0.0 && dtemp.y == 0.0 && dtemp.z == 0.0)
			Reverse(m_all4Edge[3]);
	}

	//�߽�Coons��ֵ����
	void FitBSplineSurface::InterpolateBaseSurface()
	{
		m_baseSurface.clear();
		m_EdgePara.clear();
		//�߽����		
		FitBSpline fl;
		varray<Vec4> ecpi;
		varray<double> parai;
		m_baseSurface.udegree = m_udegree;
		m_baseSurface.vdegree = m_vdegree;
		//e0
		fl.FittingBspl(m_all4Edge[0], m_udegree, m_uCtrlPtsNum);
		ecpi = fl.m_CtrlPts;
		m_EdgeCtrlPts.push_back(ecpi);
		m_baseSurface.uKnots = fl.m_Knots;
		parai = fl.m_oriPtsPara;
		m_EdgePara.push_back(parai);

		//e1
		fl.FittingBspl(m_all4Edge[1], m_vdegree, m_vCtrlPtsNum);
		ecpi.clear();
		ecpi = fl.m_CtrlPts;
		m_EdgeCtrlPts.push_back(ecpi);
		m_baseSurface.vKnots = fl.m_Knots;
		parai = fl.m_oriPtsPara;
		m_EdgePara.push_back(parai);

		//e2
		fl.FittingBspl(m_all4Edge[2], m_baseSurface.uKnots, m_udegree, m_uCtrlPtsNum);
		ecpi.clear();
		ecpi = fl.m_CtrlPts;
		m_EdgeCtrlPts.push_back(ecpi);
		parai = fl.m_oriPtsPara;
		m_EdgePara.push_back(parai);

		//e3
		fl.FittingBspl(m_all4Edge[3], m_baseSurface.vKnots, m_vdegree, m_vCtrlPtsNum);
		ecpi.clear();
		ecpi = fl.m_CtrlPts;
		m_EdgeCtrlPts.push_back(ecpi);
		parai = fl.m_oriPtsPara;
		m_EdgePara.push_back(parai);

		//Coons��ֵ
		CoonsInterpolate(m_EdgeCtrlPts, m_baseSurface.uvCtrlPts);
	}

	//����u,v��������
	void FitBSplineSurface::CalUVDir()
	{
		int e0_s = m_all4Edge[0].size();
		m_uDir = m_all4Edge[0][e0_s - 1] - m_uvOriginPts;
		int e1_s = m_all4Edge[1].size();
		m_vDir = m_all4Edge[1][e1_s - 1] - m_uvOriginPts;
	};

	////ȷ��uvw����ϵ��Ӧ��xyz����ϵ,0=x,1=y,2=z
	//void FitBSplineSurface::UVWtoXYZ()
	//{
	//	int wd = -1;
	//	//u��v�����Ӧ���ᣬΪ���������ǵĻ��б�w������1��2��3��Ӧx,y,z
	//	int ud = 1;
	//	int vd = 1;
	//	double umax = fabs(m_uDir.x), vmax = fabs(m_vDir.x);
	//	//u��
	//	if (fabs(m_uDir.y) > umax)
	//	{
	//		umax = fabs(m_uDir.y);
	//		ud = 2;
	//	}
	//	if (fabs(m_uDir.z) > umax)
	//	{
	//		umax = fabs(m_uDir.z);
	//		ud = 3;
	//	}
	//	//v��
	//	if (fabs(m_vDir.y) > vmax)
	//	{
	//		vmax = fabs(m_vDir.y);
	//		vd = 2;
	//	}
	//	if (fabs(m_vDir.z) > vmax)
	//	{
	//		vmax = fabs(m_vDir.z);
	//		vd = 3;
	//	}
	//	//�����غ�
	//	if (ud == vd)
	//	{
	//		int vdm = 1;
	//		double vmin = fabs(m_vDir.x);
	//		if (fabs(m_vDir.y) < vmin)
	//		{
	//			vmin = fabs(m_vDir.y);
	//			vdm = 2;
	//		}
	//		if (fabs(m_vDir.z) < vmin)
	//		{
	//			vmin = fabs(m_vDir.z);
	//			vdm = 3;
	//		}
	//		vd = abs(5 - vd*vdm);
	//	}
	//	//ȷ��ͶӰ����
	//	switch (ud*vd)
	//	{
	//	case 2:
	//	{
	//		wd = 3;
	//		break;
	//	}
	//	case 3:
	//	{
	//		wd = 2;
	//		break;
	//	}
	//	case 6:
	//	{
	//		wd = 1;
	//		break;
	//	}
	//	default:
	//		break;
	//	}
	//	m_uTx.clear();
	//	m_uTx.push_back(ud - 1);
	//	m_uTx.push_back(vd - 1);
	//	m_uTx.push_back(wd - 1);
	//}
	//���º���������ģ��ʵ�ã�����ģ��û���Թ���
	//ȷ��uvw����ϵ��Ӧ��xyz����ϵ,0=x,1=y,2=z
	void FitBSplineSurface::UVWtoXYZ()
	{
		int wd = -1;
		//u��v�����Ӧ���ᣬΪ���������ǵĻ��б�w������1��2��3��Ӧx,y,z
		int ud = 1;
		int vd = 1;
		double umax = fabs(m_uDir.x), vmax = fabs(m_vDir.x);
		//u��
		if (fabs(m_uDir.y) > umax)
		{
			umax = fabs(m_uDir.y);
			ud = 2;
		}
		if (fabs(m_uDir.z) > umax)
		{
			umax = fabs(m_uDir.z);
			ud = 3;
		}
		//v��
		if (fabs(m_vDir.y) > vmax)
		{
			vmax = fabs(m_vDir.y);
			vd = 2;
		}
		if (fabs(m_vDir.z) > vmax)
		{
			vmax = fabs(m_vDir.z);
			vd = 3;
		}
		//�����غ�
		if (ud == vd)
		{
			/*int vdm = 1;
			double vmin = fabs(m_vDir.x);
			if (fabs(m_vDir.y) < vmin)
			{
				vmin = fabs(m_vDir.y);
				vdm = 2;
			}
			if (fabs(m_vDir.z) < vmin)
			{
				vmin = fabs(m_vDir.z);
				vdm = 3;
			}
			vd = abs(5 - vd*vdm);*/
			if (umax > vmax)
			{
				//vѰ�ҵڶ���ķ���
				int vdm = 1;
				double vmin = fabs(m_vDir.x);
				if (fabs(m_vDir.y) < vmin)
				{
				vmin = fabs(m_vDir.y);
				vdm = 2;
				}
				if (fabs(m_vDir.z) < vmin)
				{
				vmin = fabs(m_vDir.z);
				vdm = 3;
				}
				vd = abs(5 - vd*vdm);
			}
			else
			{
				//uѰ�ҵڶ���ķ���
				int udm = 1;
				double umin = fabs(m_uDir.x);
				if (fabs(m_uDir.y) < umin)
				{
					umin = fabs(m_uDir.y);
					udm = 2;
				}
				if (fabs(m_uDir.z) < umin)
				{
					umin = fabs(m_uDir.z);
					udm = 3;
				}
				ud = abs(5 - ud*udm);
			}
		}
		//ȷ��ͶӰ����
		switch (ud*vd)
		{
		case 2:
		{
			wd = 3;
			break;
		}
		case 3:
		{
			wd = 2;
			break;
		}
		case 6:
		{
			wd = 1;
			break;
		}
		default:
			break;
		}
		m_uTx.clear();
		m_uTx.push_back(ud - 1);
		m_uTx.push_back(vd - 1);
		m_uTx.push_back(wd - 1);
	}
	//��dir�ķ����С��������dir=0,1,2
	void FitBSplineSurface::Sorted(varray<Vec4>& ptsArr, int dir)
	{
		if(dir!=0 && dir!=1 && dir!=2)
		{
			using namespace std;
			cout << "Sorted����������ֻ��Ϊx,y��z����0�˳�" << endl;
			int c_in;
			while (cin >> c_in)
			{
				if (c_in == 0)
					exit(0);
			}
		}
		//��С��������
		int i, j;
		for (i = 0; i < ptsArr.size() - 1; i++)
		{
			for (j = i + 1; j < ptsArr.size(); j++)
			{
				if (ptsArr[i][dir] > ptsArr[j][dir])//��ǰһ�㵽ƽ�����Ⱥ�һ��Ĵ��򽻻���
				{
					Vec4 temp = ptsArr[i];
					ptsArr[i] = ptsArr[j];
					ptsArr[j] = temp;
				}
			}
		}
	}

	//��ȡ�㼯ptsArr���Ե�ptsΪ���ģ���a*b���ο��ڵĵ�
	void FitBSplineSurface::GetNearPts(const varray<Vec4>& ptsArr, Vec4 pts, double a, double b, varray<Vec4>& nearPts)
	{
		nearPts.clear();
		double umin = pts[m_uTx[0]] - a / 2;
		double umax = pts[m_uTx[0]] + a / 2;
		double vmin = pts[m_uTx[1]] - b / 2;
		double vmax = pts[m_uTx[1]] + b / 2;
		for (int i = 0; i < ptsArr.size(); ++i)
		{
			double pu = ptsArr[i][m_uTx[0]];
			double pv = ptsArr[i][m_uTx[1]];
			if (pu >= umin && pu <= umax && pv >= vmin && pv <= vmax)
				nearPts.push_back(ptsArr[i]);
		}
	}

	//���������ռ��Ӧ��
	void FitBSplineSurface::Segmentation(const varray<double>& uPara, const varray<double>& vPara, varray<varray<Vec4>>& ptsArr)
	{
		ptsArr.clear();
		for (int j = 0; j < vPara.size(); ++j)
		{
			varray<Vec4> ptsArri;
			for (int i = 0; i < uPara.size(); ++i)
			{
				Vec4 pts = BSF(m_baseSurface, uPara[i], vPara[j]);
				ptsArri.push_back(pts);
			}
			ptsArr.push_back(ptsArri);
		}
	}

	//���hardy'sϵ��C
	//ptsArr:���������ݵ�
	//b:hardy��ֵ���ڳ���
	void FitBSplineSurface::CalHardysC(const varray<Vec4>& ptsArr, double b, MatrixXd& C)
	{
		int L = ptsArr.size();
		MatrixXd A(L, L), Z(L, 1);
		//A
		for (int r = 0; r < L; ++r)
		{
			double xr = ptsArr[r][m_uTx[0]];
			double yr = ptsArr[r][m_uTx[1]];
			for (int l = 0; l < L; ++l)
			{
				double xl = ptsArr[l][m_uTx[0]];
				double yl = ptsArr[l][m_uTx[1]];
				double arl = sqrt((xr - xl)*(xr - xl) + (yr - yl)*(yr - yl) + b*b);
				A(r, l) = arl;
			}
			
		}
		//Z
		for (int i = 0; i < L; ++i)
		{
			double zr = ptsArr[i][m_uTx[2]];
			Z(i,0) = zr;
		}
		//C
		C = A.colPivHouseholderQr().solve(Z);
	}

	//Hardy's˫���β�ֵ
	Vec4 FitBSplineSurface::HardysInterpolate(const Vec4 pts, const varray<Vec4>& ptsArr, const MatrixXd& C, const double b)
	{
		double x = pts[m_uTx[0]];
		double y = pts[m_uTx[1]];
		int L = ptsArr.size();
		MatrixXd A(1, L);
		for (int i = 0; i < L; ++i)
		{
			double xl = ptsArr[i][m_uTx[0]];
			double yl = ptsArr[i][m_uTx[1]];
			A(0,i) = sqrt((x - xl)*(x - xl) + (y - yl)*(y - yl) + b*b);
		}
		Vec4 res(0,0,0);
		res[m_uTx[0]] = x;
		res[m_uTx[1]] = y;
		res[m_uTx[2]] = (A*C)(0,0);
		return res;
	}

	//���㷶Χ
	void FitBSplineSurface::CalExtreme(double extrm[4])
	{
		double ud = m_uDir[m_uTx[0]];
		double vd = m_vDir[m_uTx[1]];
		varray<Vec4> ptsArr1 = m_all4Edge[1];
		varray<Vec4> ptsArr3 = m_all4Edge[3];
		Sorted(ptsArr1, m_uTx[0]);
		Sorted(ptsArr3, m_uTx[0]);

		varray<Vec4> ptsArr0 = m_all4Edge[0];
		varray<Vec4> ptsArr2 = m_all4Edge[2];
		Sorted(ptsArr0, m_uTx[1]);
		Sorted(ptsArr2, m_uTx[1]);
		if (ud > 0)
		{
			extrm[0] = ptsArr1[0][m_uTx[0]];
			extrm[1] = ptsArr3[ptsArr3.size() - 1][m_uTx[0]];
		}
		else
		{
			extrm[0] = ptsArr3[0][m_uTx[0]];
			extrm[1] = ptsArr1[ptsArr1.size() - 1][m_uTx[0]];
		}
		if (vd > 0)
		{
			extrm[2] = ptsArr0[0][m_uTx[1]];
			extrm[3] = ptsArr2[ptsArr2.size() - 1][m_uTx[1]];
		}
		else
		{
			extrm[2] = ptsArr2[0][m_uTx[1]];
			extrm[3] = ptsArr0[ptsArr0.size() - 1][m_uTx[1]];
		}
	}

	//��ֵ�����
	void FitBSplineSurface::InterpolatGridPts(int uNum, int vNum)
	{
		//������ָ�
		varray<double> uPara, vPara;
		for (int i = 0; i <= uNum; ++i)
			uPara.push_back(1.0*i / uNum);
		for (int i = 0; i <= vNum; ++i)
			vPara.push_back(1.0*i / vNum);
		//��ֵ���������
		varray<varray<Vec4>> uvGrid;
		Segmentation(uPara, vPara, uvGrid);
		//��ÿһ��������ֵ
		//��Χ
		double extrm[4];
		CalExtreme(extrm);
		double a = 2 * (extrm[1] - extrm[0]) / uNum;
		double b = 2 * (extrm[3] - extrm[2]) / vNum;
		FitBSpline fl1, fl3;
		fl1.m_Degree = m_vdegree;
		fl1.m_CtrlPts = m_EdgeCtrlPts[1];
		fl1.m_Knots = m_baseSurface.vKnots;
		fl3.m_Degree = m_vdegree;
	    fl3.m_CtrlPts = m_EdgeCtrlPts[3];
		fl3.m_Knots = m_baseSurface.vKnots;

		m_uPtsVec.clear();
		for (int i = 0; i < uvGrid.size(); ++i)
		{
			varray<Vec4> UPVi;
			//��λe1��Ӧvֵ��
			Vec4 hp = fl1.Bspl(vPara[i]);
			UPVi.push_back(hp);
			//hardy��ֵ
			for (int j = 1; j < uvGrid[i].size() - 1; ++j)
			{
				varray<Vec4> nearpts;
				GetNearPts(m_allInputPts, uvGrid[i][j], a, b, nearpts);
				int L = nearpts.size();
				MatrixXd C(L,1);
				CalHardysC(nearpts, 1, C);
				hp = HardysInterpolate(uvGrid[i][j], nearpts, C, 1);
				UPVi.push_back(hp);
			}
			//ĩλe3��Ӧvֵ��
			hp = fl3.Bspl(vPara[i]);
			UPVi.push_back(hp);
			m_uPtsVec.push_back(UPVi);
		}
	}

	//���u����
	void FitBSplineSurface::FittingUdir()
	{
		FitBSpline fbsl;
		varray<varray<double>> knots;//u�����������߽ڵ�ʸ����
		int uPtVsz = m_uPtsVec.size();
		m_uCtrlPts.clear();
		m_uPara.clear();
		for (int i = 0; i < uPtVsz; ++i)
		{
			fbsl.FittingBspl(m_uPtsVec[i], m_udegree, m_uCtrlPtsNum);
			knots.push_back(fbsl.m_Knots);
			m_uPara.push_back(fbsl.m_oriPtsPara);
			m_uCtrlPts.push_back(fbsl.m_CtrlPts);
		}

		//u����ڵ���ƽ��ֵ
		m_uKnots.clear();
		int knots_sz = knots[0].size();
		for (int i = 0; i < knots_sz; ++i)
		{
			double sumKnot = 0.0;//iλ�����нڵ�ĺ�
			for (int j = 0; j < knots.size(); ++j)
			{
				sumKnot += knots[j][i];
			}
			double knoti = sumKnot / knots.size();
			m_uKnots.push_back(knoti);
		}
	}

	//���v����
	void FitBSplineSurface::FittingVdir()
	{
		varray<varray<Vec4>> uCtrlPts_T = m_uCtrlPts;//u����Ƶ�ת��
		TransPose(uCtrlPts_T);

		varray<varray<double>> knots;
		m_uvCtrlPts.clear();
		int uPtsVec_Tsz = uCtrlPts_T.size();
		FitBSpline fbsl;
		for (int i = 0; i < uPtsVec_Tsz; ++i)
		{
			fbsl.FittingBspl(uCtrlPts_T[i], m_vdegree, m_vCtrlPtsNum);
			knots.push_back(fbsl.m_Knots);
			m_uvCtrlPts.push_back(fbsl.m_CtrlPts);
		}
		//v����ڵ���ƽ��ֵ
		m_vKnots.clear();
		int knot_sz = knots[0].size();
		for (int i = 0; i < knot_sz; ++i)
		{
			double sumKnot = 0.0;//iλ�����нڵ�ĺ�
			for (int j = 0; j < knots.size(); ++j)
			{
				sumKnot += knots[j][i];
			}
			double knoti = sumKnot / knots.size();
			m_vKnots.push_back(knoti);
		}
	}


	//�ж�pnt�Ƿ��ڵ㼯ptsArr��,����������-1Ϊ���ڵ㼯��
	int FitBSplineSurface::InArr(Vec4 pnt, const varray<Vec4>& ptsArr)
	{
		int res = -1;
		for (int i = 0; i < ptsArr.size(); ++i)
		{
			Vec4 dp = pnt - ptsArr[i];
			if (dp.Magnitude() < 0.0001)
				res = i;
		}
		return res;
	}

	//v������8����߽���ж�
	bool FitBSplineSurface::IsVEdgePts(Vec4 pts, const varray<Vec4>& ptsArr)
	{
		//�����
		double dv = 0;
		for (int i = 1; i < m_all4Edge[1].size(); ++i)
		{
			double ds = fabs(m_all4Edge[1][i][m_uTx[1]] - m_all4Edge[1][i - 1][m_uTx[1]]);
			if (ds > dv)
				dv = ds;
		}
		for (int i = 1; i < m_all4Edge[3].size(); ++i)
		{
			double ds = fabs(m_all4Edge[3][i][m_uTx[1]] - m_all4Edge[3][i - 1][m_uTx[1]]);
			if (ds > dv)
				dv = ds;
		}
		double R = 1.5 * dv;//����뾶
		varray<Vec4> Rpts;//�����ڵ㼯
		double r = 1000000;//��С����
		Vec4 rp;//���������
		Vec4 rdp;//���������ʸ����
		for (int i = 0; i < ptsArr.size(); ++i)
		{
			int idx = InArr(pts, ptsArr);
			if (idx != i)
			{
				Vec4 dp = ptsArr[i] - pts;
				double dis = sqrt(dp[m_uTx[0]] * dp[m_uTx[0]] + dp[m_uTx[1]] * dp[m_uTx[1]]);
				if (dis <= R)
				{
					Rpts.push_back(ptsArr[i]);
					if (dis < r)
					{
						r = dis;
						rp = ptsArr[i];
						rdp = dp;
					}
				}
			}
		}
		//��������
		double umin, umax, vmin, vmax;
		umin = pts[m_uTx[0]] - r;
		umax = pts[m_uTx[0]] + r;
		vmin = pts[m_uTx[1]] - r;
		vmax = pts[m_uTx[1]] + r;
		//�ж�
		int qy[8] = { 0,0,0,0,0,0,0,0 };
		for (int i = 0; i < Rpts.size(); ++i)
		{
			Vec4 rptsi = Rpts[i];
			double urptsi = rptsi[m_uTx[0]];
			double vrptsi = rptsi[m_uTx[1]];

			if (urptsi > umin && urptsi < umax && vrptsi > vmax)//����0
				qy[0] += 1;
			else if (urptsi <= umin && vrptsi >= vmax)//����1
				qy[1] += 1;
			else if (vrptsi > vmin && vrptsi < vmax && urptsi < umin)//����2
				qy[2] += 1;
			else if (urptsi <= umin && vrptsi <= vmin)//����3
				qy[3] += 1;
			else if (urptsi > umin && urptsi < umax && vrptsi < vmin)//����4
				qy[4] += 1;
			else if (urptsi >= umax && vrptsi <= vmin)//����5
				qy[5] += 1;
			else if (vrptsi > vmin && vrptsi < vmax && urptsi > umax)//����6
				qy[6] += 1;
			else if (urptsi >= umax && vrptsi >= vmax)//����7
				qy[7] += 1;
			else
			{
				Vec4 dp_rptsi = rptsi - pts;
				Vec4 prdp = Vec4(fabs(dp_rptsi.x), fabs(dp_rptsi.y), fabs(dp_rptsi.z));
				if (prdp[m_uTx[1]] >= prdp[m_uTx[0]])
				{
					if (dp_rptsi[m_uTx[1]] > 0)
						qy[0] += 1;
					else
						qy[4] += 1;
				}
				else
				{
					if (dp_rptsi[m_uTx[0]] > 0)
						qy[6] += 1;
					else
						qy[2] += 1;
				}
			}
		}
		//v��������Ϊ0��Ϊ�߽��
		bool res = false;
		int vdir = m_vDir[m_uTx[1]];
		if (vdir > 0)
		{
			if (qy[4] + qy[3] == 0 || qy[5] + qy[4] == 0)
				res = true;
		}
		else
		{
			if (qy[1] + qy[0] == 0 || qy[0] + qy[7] == 0)
				res = true;
		}
		/*if (qy[0] + qy[1] == 0 || qy[1] + qy[2] == 0 || qy[2] + qy[3] == 0 || qy[3] + qy[4] == 0 ||
		qy[4] + qy[5] == 0 || qy[5] + qy[6] == 0 || qy[6] + qy[7] == 0 || qy[7] + qy[0] == 0)
		res = true;*/
		return res;
	}

	//����̶�v��u����ĵ㼯������m_uPtsVec
	void FitBSplineSurface::CalUPtsVec()
	{
		m_uPtsVec.clear();
		varray<Vec4> uPtsVi;//�����������У�
							//����e0
		uPtsVi = m_all4Edge[0];
		Sorted(uPtsVi, m_uTx[0]);
		m_uPtsVec.push_back(uPtsVi);

		//���������㼯
		varray<Vec4> oriptsForU;
		for (int i = 0; i < m_allInputPts.size(); ++i)//ȥ���߽��
		{
			Vec4 oripts = m_allInputPts[i];
			int isedge = 0;
			for (int j = 0; j < m_all4Edge.size(); ++j)
			{
				int idx = InArr(oripts, m_all4Edge[j]);
				if (idx != -1)
					isedge = 1;
			}
			if (isedge == 0)
				oriptsForU.push_back(oripts);
		}

		int uPV_sz = m_all4Edge[0].size() - 2;//ÿ�д�СΪe0
											  //ɨ�跽��
		Sorted(oriptsForU, m_uTx[1]);//v��С����
		if (m_vDir[m_uTx[1]] < 0)
			Reverse(oriptsForU);//v�ɴ�С

								//ɨ��
		uPtsVi.clear();
		for (int j = 0; j < oriptsForU.size(); ++j)
		{
			Vec4 pnt = oriptsForU[j];
			bool isvp = IsVEdgePts(pnt, oriptsForU);
			if (isvp && InArr(pnt, uPtsVi) == -1)//�ֽ�����
			{
				uPtsVi.push_back(pnt);
				Sorted(uPtsVi, m_uTx[0]);
				if (uPtsVi.size() == uPV_sz)//��һ��
				{
					j = -1;
					//��Ӷ�Ӧ�б߽��
					int upv_sz = m_uPtsVec.size();
					uPtsVi.insert(uPtsVi.begin(), m_all4Edge[1][upv_sz]);//e1�ĵ�upv_sz����
					uPtsVi.push_back(m_all4Edge[3][upv_sz]);//e3�ĵ�upv_sz����
					Sorted(uPtsVi, m_uTx[0]);
					m_uPtsVec.push_back(uPtsVi);
					//ɾ��
					for (int k = 0; k < uPtsVi.size(); ++k)
					{
						Vec4 pi = uPtsVi[k];
						int pos = InArr(pi, oriptsForU);
						if (pos != -1)
							oriptsForU.erase(oriptsForU.begin() + pos);
					}
					uPtsVi.clear();
					//�����ʣһ��
					if (oriptsForU.size() == uPV_sz)
					{
						Sorted(oriptsForU, m_uTx[0]);
						//��Ӷ�Ӧ�б߽��
						int upv_sz = m_uPtsVec.size();
						oriptsForU.insert(oriptsForU.begin(), m_all4Edge[1][upv_sz]);//e1�ĵ�upv_sz����
						oriptsForU.push_back(m_all4Edge[3][upv_sz]);//e3�ĵ�upv_sz����
						Sorted(oriptsForU, m_uTx[0]);
						m_uPtsVec.push_back(oriptsForU);
						break;
					}
				}//��һ��
			}//�ֽ�����
			 //���е�ɨ����ϣ�����û����ɷ���
			if (j == oriptsForU.size() - 1 && m_uPtsVec.size() != m_all4Edge[1].size() - 1)
			{
				using namespace std;
				cout << "CalUPtsVec�������е�ɨ����ϣ�����û����ɷ��飡\n��0�˳�\n��1����" << endl;
				int c_in;
				while (cin >> c_in)
				{
					if (c_in == 0)
						exit(0);
					else if (c_in == 1)
						break;
				}
			}
		}//ɨ�����
		uPtsVi.clear();
		//β��e2
		uPtsVi = m_all4Edge[2];
		Sorted(uPtsVi, m_uTx[0]);
		m_uPtsVec.push_back(uPtsVi);

		//�ָ�˳��
		if (m_uDir[m_uTx[0]] < 0)
		{
			for (int i = 0; i < m_uPtsVec.size(); ++i)
			{
				Reverse(m_uPtsVec[i]);
			}
		}
		//m_uPtsVec�������
	}

	//v���������
	void FitBSplineSurface::VParaPts()
	{
		varray<varray<Vec4>> vPtsVec;
		vPtsVec = m_uPtsVec;
		TransPose(vPtsVec);
		//m_vPtsVec�������
		FitBSpline fbsl;
		m_vPara.clear();
		int uPtsVec_Tsz = vPtsVec.size();
		for (int i = 0; i < uPtsVec_Tsz; ++i)
		{
			fbsl.FittingBspl(vPtsVec[i], m_vdegree, m_vCtrlPtsNum);
			m_vPara.push_back(fbsl.m_oriPtsPara);
		}
		TransPose(m_vPara);
	}

	//����������
	void FitBSplineSurface::CalSurFaceFittingErr()
	{
		VParaPts();
		varray<Vec4> pnt;
		varray<double> pntu;
		varray<double> pntv;
		//ASSERT(m_uPtsVec.size() <= m_uPara.size() && m_uPtsVec.size() <= vparat.size());
		for (int i = 0; i < m_uPtsVec.size(); i++)
		{
			//ASSERT(m_uPtsVec.at(i).size() <= m_uPara.at(i).size() && m_uPtsVec.at(i).size() <= vparat.at(i).size());
			for (int j = 0; j < m_uPtsVec[i].size(); j++)
			{
				pnt.push_back(m_uPtsVec[i][j]);
				pntu.push_back(m_uPara[i][j]);
				pntv.push_back(m_vPara[i][j]);
			}
		}
		double sumeps = 0.0;
		for (int i = 0; i < pnt.size(); ++i)
		{
			sumeps += (pnt[i] - BSF(m_baseSurface,pntu[i], pntv[i])).SquareMagnitude();
		}
		m_FittingErr = sqrt(sumeps / pnt.size());
	}

	/*�������
	allInputdata����ɢ�㼯
	all4edge���߽�㼯
	degree�����ߴ���
	uCtrlNum:u������Ƶ�����
	vCtrlNum:v������Ƶ�����
	*/
	void FitBSplineSurface::FittingBsurface(const varray<Vec4>& allInputdata, const varray<varray<Vec4>>& all4edge,
		int udegree/* = 3 */, int vdegree/* = 3 */, int uCtrlNum/* = 4*/, int vCtrlNum /* = 4*/)
	{
		using namespace std;
		int nsize = allInputdata.size();//��ɢ������

		if (allInputdata.empty() || all4edge.empty())
		{
			cout << "������ɢ�㼯��߽�Ϊ�գ���0�˳�" << endl;
			while (cin)
				return;
		}
		else
		{
			m_allInputPts.clear();
			m_all4Edge.clear();
			m_allInputPts = allInputdata;
			m_all4Edge = all4edge;
			m_udegree = udegree;
			m_vdegree = vdegree;
			m_uCtrlPtsNum = uCtrlNum;
			m_vCtrlPtsNum = vCtrlNum;

			m_baseSurface.udegree = m_udegree;
			m_baseSurface.vdegree = m_vdegree;

			FindOriginPts();
			AligningEdge();
			CalUVDir();
			UVWtoXYZ();
			//�ж����ݵ�����
			bool isUinterpolate = (m_all4Edge[0].size() < m_uCtrlPtsNum);
			bool isVinterpolate = (m_all4Edge[1].size() < m_vCtrlPtsNum);
			FitBSpline fl;
			if (isUinterpolate || isVinterpolate)//��ֵ
			{
				//u�����ֵ
				fl.Interpolate(m_all4Edge[0], m_uCtrlPtsNum);
				fl.Interpolate(m_all4Edge[2], m_uCtrlPtsNum);
				//v�����ֵ
				fl.Interpolate(m_all4Edge[1], m_vCtrlPtsNum);
				fl.Interpolate(m_all4Edge[3], m_vCtrlPtsNum);
			}
			if (nsize != m_all4Edge[0].size()*m_all4Edge[1].size())//Coons��ֵ��������
			{
				//�߽����
				varray<varray<Vec4>> EndgCtrlPts;//�߽���Ƶ㼯
				varray<Vec4> ecpi;
				//e0
				fl.FittingBspl(m_all4Edge[0], m_udegree, m_uCtrlPtsNum);
				ecpi = fl.m_CtrlPts;
				EndgCtrlPts.push_back(ecpi);
				m_uKnots.clear();
				m_uKnots = fl.m_Knots;
				//e1
				fl.FittingBspl(m_all4Edge[1], m_vdegree, m_vCtrlPtsNum);
				ecpi = fl.m_CtrlPts;
				EndgCtrlPts.push_back(ecpi);
				m_vKnots.clear();
				m_vKnots = fl.m_Knots;
				//e2
				fl.FittingBspl(m_all4Edge[2], m_udegree, m_uCtrlPtsNum);
				ecpi = fl.m_CtrlPts;
				EndgCtrlPts.push_back(ecpi);
				//e3
				fl.FittingBspl(m_all4Edge[3], m_vdegree, m_vCtrlPtsNum);
				ecpi = fl.m_CtrlPts;
				EndgCtrlPts.push_back(ecpi);
				//Coons��ֵ
				CoonsInterpolate(EndgCtrlPts, m_uvCtrlPts);
			}
			else//�ǲ�ֵ
			{
				CalUPtsVec();
				FittingUdir();
				FittingVdir();
				m_baseSurface.uKnots = m_uKnots;
				m_baseSurface.vKnots = m_vKnots;
				m_baseSurface.uvCtrlPts = m_uvCtrlPts;
				CalSurFaceFittingErr();
			}
		}
	}

}// !namespace bsl
