#pragma once
#include <../stdafx.h>
#include <../globalFunc.h>

using std::vector;
using std::string;
using std::shared_ptr;
using std::make_shared;

namespace YN
{

	class NurbsLine
	{
	public:
		NurbsLine();
		NurbsLine(const NurbsLine &);
		NurbsLine& operator=(const NurbsLine&);
		//���Ƶ�������������ڵ�ʸ�������Ƶ�
		NurbsLine(int, int, vector<double>, vector<point4d>);
		~NurbsLine();

		int GetUDegree()const;
		shared_ptr<vector<double>> GetUKonts()const;
		int GetUNum()const;
		shared_ptr<vector<point4d>> GetControlPointer()const;

		bool SetControlPoint(const vector<point4d>);
		bool SetControlPoint(shared_ptr<vector<point4d>>);
		bool SetUDegree(const int);
		bool SetUKonts(const vector<double>);
		bool SetUKonts(shared_ptr<vector<double>>);
		bool SetUNum(const int);
		//�������Ƿ���ȷ
		//���Ƶ����+1+����=�ڵ�ʸ������
		bool isParameterCorrect(int, int, int);

		/*����ڵ��±�
		x���ڵ�ֵ
		degree������
		CtrlPtsNum�����Ƶ�����
		knots���ڵ�ʸ��*/
		int FindSpan(const double x, const int degree, const int CtrlPtsNum, const vector<double>& knots)const;

		/*���ݲ���ֵ���ڵ��±꣬���������
		u������ֵ
		k���ڵ��±�
		degree������
		knots���ڵ�ʸ��
		N�����صĻ�����*/
		void BasisFuns(const double u, const int k, const int degree, const vector<double>& knots,
			vector<double>& N)const;

		/*u�����л�����
		u������ֵ
		k���ڵ��±�
		degree������
		knots���ڵ�ʸ��
		ndu�����ص����л�����*/
		void AllBasisFuns(const double u, const int k, const int degree,
			const vector<double>& knots, vector<vector<double>>& ndu);

		/*������n�׵�
		u������ֵ
		k���ڵ��±�
		degree������
		n����ʸ����
		knots���ڵ�ʸ��
		basisDu��������n�׵�*/
		void DerBasisFuns(const double u, const int k, const int degree, const int n, const vector<double>& knots,
			vector<vector<double>>& basisDu);

		/*����u�ڵ��Ӧ�����������
		u���ڵ����*/
		point3d GetLinePoint(const double u);

		/*��ȡ���ϵĵ�
		Unum�����ߵ�����
		linePts�����ߵ�	*/
		void CalLinePoint(const int Unum, vector<point3d>& linePts);

		/*����ʸֵ����A(u)����p�׵���ʸֵ����A(u)��NURBS���߷���ʽ
		u������
		n����ʸ����
		Der��A(u)����p�׵�*/
		void PtsDerivsAu(const double u, const int n, vector<point4d>& Der);

		/*������u������n�׵�
		u������
		n����ʸ����
		Der������n�׵� */
		void PtsDerivs(const double u, const int n, vector<point4d>& Der);

		/*��������n�׵�
		step��������ȷֲ�����
		Der�����ߵ��Ӧ��n�׵�*/
		void CurveDerivs(const int n, const double step, vector<vector<point3d>>& Der);

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
		virtual void KnotsRefine(const vector<double>& u);

		/* ��������
		degree�����׺����*/
		virtual void DegreeElevate(const int degree);

		//�ڽڵ�u�������߽��зָ�
		void Segmentation(const double u, vector<NurbsLine>& lines);

	public:
		//u�������
		int _u_Degree;
		//u����ڵ�ʸ��
		shared_ptr<vector<double>> _u_Knots;
		//u������Ƶ����
		int _u_Num;
		//���Ƶ�
		shared_ptr<vector<point4d>> _ControlPts;
	};

	class NurbsSurface :public NurbsLine
	{
	public:
		NurbsSurface();
		NurbsSurface(int, int, int, int, vector<double>, vector<double>, vector<point4d>);
		NurbsSurface(const NurbsSurface&);
		NurbsSurface& operator=(const NurbsSurface&);
		~NurbsSurface();

		int GetVDegree()const;
		shared_ptr<vector<double>> GetVKonts()const;
		int GetVNum()const;
		point4d getControlPoint(int u, int v);

		bool SetVDegree(const int);
		bool SetVKonts(const vector<double>);
		bool SetVKonts(shared_ptr<vector<double>>);
		bool SetVNum(const int);

		void SetSurface(const int uDegree, const int vDegree, const int uNum, const int vNum,
			const vector<double>& uKnots, const vector<double>& vKnots);

		//�ж��������Ƿ���ͬ
		bool isTwoSurfaceSame(NurbsSurface);

		//����������֮��Ĺ�˹������
		float GetDistanceBetweenTwoSurface(NurbsSurface);

		//���㣨u��v����Ӧ�������ϵĵ�
		point3d GetSurFacePoint(const double u, const double v)const;

		//�����ı�����Ƭ��ʾ����
		//num:�÷������������
		//quads:��Ƭ����
		//lines:�Ȳ�������
		threadParam CalQuads(const int Unum, const int Vnum)const;

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

		//���ݿ��Ƶ��ά��ż���һά���
		int CtrlPtsIdx(const int uIdx, const int vIdx);

		//���Ƶ�����ת��ΪU-V
		void OrderCtrlPts();

		//���Ƶ�����ת��ΪU-V
		void OrderCtrlPts(NurbsSurface& sf);

	public:
		//v�������
		int _v_Degree;
		//v����ڵ�ʸ��
		shared_ptr<vector<double>> _v_Knots;
		//u������Ƶ����
		int _v_Num;
	};


	class NurbsVol :public NurbsSurface
	{
	public:
		NurbsVol();
		~NurbsVol();
		NurbsVol(const NurbsVol &);
		NurbsVol& operator=(const NurbsVol&);

		int GetWDegree()const;
		shared_ptr<vector<double>> GetWKonts()const;
		int GetWNum()const;
		point4d getControlPoint(int u, int v, int w);

		bool SetWDegree(const int);
		bool SetWKonts(const vector<double>);
		bool SetWKonts(shared_ptr<vector<double>>);
		bool SetWNum(const int);
		void SetVol(const int uDegree, const int vDegree, const int wDegree, const int uNum, const int vNum, const int wNum,
			const varray<double>& uKnots, const varray<double>& vKnots, const varray<double>& wKnots);

		//(u,v,w)��Ӧ�����ϵĵ�,
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
		void KnotsRefine(varray<double> Uknot, varray<double> Vknot, varray<double> Wknot);

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

		//����
		//path:·��
		//surfaces0:·������������
		//surfaces1:·���յ��������
		void LoftingNurbsVol(const NurbsLine& path, const NurbsSurface& surfaces0, const NurbsSurface& surfaces1);

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


	public:
		//v�������
		int _w_Degree;
		//v����ڵ�ʸ��
		shared_ptr<vector<double>> _w_Knots;
		//u������Ƶ����
		int _w_Num;
	};

};