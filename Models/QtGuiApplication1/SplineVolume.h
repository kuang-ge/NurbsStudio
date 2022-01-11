#if !defined(_SPLINEVOLUME_H_INCLUDED)
#define _SPLINEVOLUME_H_INCLUDED


#include "SplineSurface.h"
#include"varray.h"
#include"XVec.h"
using namespace base;
namespace base
{
	//ÿ���ı���Ƭ��
	enum surfacePosition{
		U0sfPos = 0,
		U1sfPos = 1,
		V0sfPos = 2,
		V1sfPos = 3,
		W0sfPos = 4,
		W1sfPos = 5
	};
	enum volumeRenderMode{
		iRenderBoundary_ControlPoints = 0,
		iRenderBoundary_Patch = 1,
		iRenderBoundary_WithISOElements = 2,
		iRenderBoundary_All = 3,

		iRenderVolume_ControlPoints = 4,
		iRenderVolume_TinyPatch = 5,
		iRenderVolume_WithISOElements = 6,
		iRenderVolume_All = 7,
        iRenderVolume_ISOWireFrame = 8
	};
#define SEGMENTNUM  60;
class VolumeVertex//��ģ���е�һ�����Ƶ�������Ϣ
{
public:
	Vec4  m_pt;      //�����
	Vec4  m_oript;   //ԭʼ���Ƶ㣬���ڱ����������ǰ�Ŀ��Ƶ㡣
	Vec4  m_matval;  //���ϲ���
	Vec4  m_displacement;  //��ʩ���غɺ��γɵ�λ�ơ�  �ֱ��ʾux��uy,uz
	Vec4  m_stress;    //�ĵ������غɡ�     �ֱ��ʾsigma(x),sigma(y),sigma(z).
	float m_fvonMissStressMat, m_fDisplaceMat;    //��ʾ���������Сֵ�õ���Ӧ����ɫֵ��λ����ɫֵ����Χ��[0,1].��Ҫ��ʾλ�ƻ���Ӧ��ʱ��ֻ�轫��Ӧ��ֵ������ɫ�������ɡ�
public:
	VolumeVertex();
	VolumeVertex(Vec4& vt);
	VolumeVertex(Vec4& vt,Vec4& matt);
	float getVonmissStress(){return sqrt((m_stress.x-m_stress.y)*(m_stress.x-m_stress.y)/2 + (m_stress.y-m_stress.z)*(m_stress.y-m_stress.z)/2 + (m_stress.z-m_stress.x)*(m_stress.z-m_stress.x)/2);};
	~VolumeVertex();    
};
class SplineVolume: public SplineBase
{
public:
	SplineVolume(void);
	~SplineVolume(void);
public:
	int   m_uNum, m_uDegree, m_uRenderNum;
	int   m_vNum, m_vDegree, m_vRenderNum;
	int   m_wNum, m_wDegree, m_wRenderNum;
	bool  m_bHeterogeneous;//������
	varray<Vec4>			 m_CtrlPts;
	varray<VolumeVertex> m_vAllCtrlPts;  //���Ƶ�
	varray<VolumeVertex> m_vVolumePts;   //������ʾ�Ľڵ㡣
	varray<int>  m_CtrlPtsIDinKKKMatrix;   //ÿ�����Ƶ����ܸپ����е�λ�á� ����ǵ�Ƭ��Ĭ�Ͼ���0~m_uNum*m_vNum*m_wNum-1;����Ƕ�Ƭ����Ҫ���±�š�
	bool  m_bHasInterpolated;//������
	varray<double> m_uKnots;
	varray<double> m_vKnots;
	varray<double> m_wKnots;
	varray<SplineSurface> m_6boundarySurface;  //����ǰ���Ƕ�Ӧm_6PatchIdx����Ӧ����ϵ�surface,�洢������ת�������
	bool  m_bHasIGASolution;
	float m_maxJocbian, m_minJocbian;
	int   m_minorJocbianNum;
	double m_isovalAdd;
public:
	void  Clear();
	//��ʾ����
	//void  Render(bool vertshow,bool boundarysurface=false,bool volumedisplay = false);
	/*void  RenderbyMode(int iVolumeMode,int iSurfaceMode);
	void  RenderPatch(bool bBoundary,int iSurfaceMode);
	void  RenderCtrlPts(bool bOnlyBoudary);
	void  RenderIsopoints(bool bOnlyBoudary);
	void  RenderIsoCurves(bool bOnlyBoudary);
	void  RenderIsoSurfaces(bool bOnlyBoudary,float alpha);
	void  RenderTinyCube();*/

	void  CreateAllVolumeRenderPts(int uN, int vN, int wN);//������ʵ��
	void  ReCreateAllVolumeRenderPts();//���¼�����ʵ��
	void  InterpolateInnerControlPtsByBoundarysueface();//ͨ���߽����ֵ�ڲ����Ƶ�
	//evaluate operation
	void  GetAllControlPointsVec4(varray<Vec4>& allVecs);
	Vec4  GetPoint(float u, float v, float w, Vec4& matPt);
	Vec4  GetControlPoint(int ui, int vi, int wi);
	Vec4  GetControlPointMatvalue(int ui, int vi, int wi);
	int   GetControlPointIndex(int ui, int vi, int wi);
	bool  GetIJKIndex(int i, int& ui, int&vi, int& wi, int unum, int vnum, int wnum);
	bool  IsBoundaryPoint(int ui, int vi, int wi, int unum, int vnum, int wnum);
	int   GetRenderPointIndex(int ui, int vi, int wi);
	void  SetControlPtNum(int un, int vn, int wn);//���ÿ��Ƶ����
	void  SetDegree(int ud, int vd, int wd);//���ý���
	void  SetControlPoints(varray<Vec4>& vts);//���ÿ��Ƶ�
	void  SetControlPoints(varray<Vec4>& vts, varray<Vec4>& mts);
	void  SetControlPoints(varray<Vec4>& vts, int udegree, int vdegree, int wdegree, int un, int vn, int wn);//�൱�ڹ��캯��
	bool  JudgeNeedInterpolate();

	//��ȡ�߽�������ĵ�
	void  Get6SurfaceCtrlPtsFromPts(varray<varray<Vec4> >& surfacePts, int uvwmode);//0,1 for u direction surface, 2,3 for v direction surface,4,5 for w direction surface.
	//��ȡ
	void  GetBoundBox(Vec4& leftbot, Vec4& rightTop, int istatus = 0); //0 for render points, 1 for control points.

	//ƽ��ֵ��
	//void  InterpolateAllControlPtsByMeanValue(bool isinitial=true);
	//void  LoopInterpolateAllControlPtsByMeanValue(bool isinitial=true);  //mean value coordinate.	                              
	int   GetIndexForVolumeMatrix(int iu, int iv, int iw, int su, int sv, int sw, int localidx, bool& boundaryflag);
	void  GetLocalScaleForPoint(double scale[], int ciu, int civ, int ciw, bool laplacian = false, bool isinitial = true);
	void  ModifyLocalIndexForCube(int hindex, int& mu, int& mv, int& mw);

	//��ɢ������˹����
	/*void  InterpolateallControlPtsByLaplacian();*/

	//��˹��
	//the following method is referred to ��Farin G.."Discrete coons patches."��
	void  InterpolateAllControlPtsByCoons();  //by Gang Xu's method named discrete coons method.

	//����ÿ������Ƶ����������
	float  GetJacobianValue(float u, float v, float w);
	float  GetFirstPartialDirivatives(float u, float v, float w, Vec4& pu, Vec4& pv, Vec4& pw);
	void   TestMeshQuanity(int segmentNum = 60/*SEGMENTNUM*/);

	void   SaveVolume();

	//�ȼ��η������
	void   SetDisplaceAndLoadvector(int ctrlPtid, Vec4& disVec, Vec4& loadVec);

	void SetVol(const int uDegree, const int vDegree, const int wDegree, const int uNum, const int vNum, const int wNum,
		const varray<double>& uKnots, const varray<double>& vKnots, const varray<double>& wKnots);

	//(u,v,w)��Ӧ�����ϵĵ�
	Vec3 GetVolPoint(const double u, const double v, const double w)const;


	//�����ı�����Ƭ��ʾ����
	//num:�÷������������
	//quads:��Ƭ����
	//lines:�Ȳ�������
	threadParamSplineVOL CalQuads(const int Unum, const int Vnum, const int Wnum,
		varray<varray<varray<Vec3>>>& quads, varray<varray<varray<Vec3>>>& lines)const;

	//���ݿ��Ƶ���ά��ż���һά���
	int CtrlPtsIdx(const int uIdx, const int vIdx, const int wIdx);

	//���Ƶ�����ת��ΪU-V-W
	void OrderCtrlPts();

	//���Ƶ�����ת��ΪU-V-W
	void OrderCtrlPts(SplineVolume& vol);

	//���Ƶ�����ת��ΪU-V-W
	void OrderCtrlPts1(SplineVolume& vol,int mode);

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
	void CreateSweepSplineVolume(const Spline& pathT, const SplineSurface& nurbsSF, const int K);

	/*ƽ��ɨ������Nurbs��ģ�ͣ����治��ֱ��·��
	pathT��ɨ��·��
	nurbsSF����ʼ����*/
	void CreateTransSweepSplineVolume(const Spline& pathT, const SplineSurface& nurbsSF);

	//����
	//path:·��
	//surfaces:��������
	void LoftingSplineVolume(const Spline& path, const varray<SplineSurface>& surfaces);
	//��ȡ�����߽���
	//dir:1->6��ʾ����
	SplineSurface GetSingleSurfaces(int dir) const;
	//����
	//path:·��
	//surfaces0:·������������
	//surfaces1:·���յ��������
	void LoftingSplineVolume(const Spline& path, const SplineSurface& surfaces0, const SplineSurface& surfaces1);
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
		varray<varray<Vec3>>& quads, varray<varray<Vec3>>& lines)const;

	/*ɨ��ʱ�Ľ���ʵ��λ��
	pathT��ɨ��·��
	K������ʵ������
	pos������ʵ��λ��
	NewKnots���½ڵ�ʸ��*/
	void InsLocation(const Spline & pathT, int K, varray<double> & pos);

	/*����ɨ��·���ϵľֲ�����ϵ
	pathT��ɨ��·��
	pos������ʵ��λ��
	TranMat��pos���ľֲ�����ϵ*/
	void LocalCoordinates(const Spline & pathT, const varray<double> & pos,
		varray<varray<Vec4>> & TranMat);



	//��COONS��ֵ
	void VolCoonsInterpolate(const varray<Spline>& EdgeLines);

	/*��ɨ��·���Խ�����б任
	nurbsSF��������Ƶ�
	TranMat���任����
	OriCoordinates:����ľֲ�����ϵ
	allNurbsSF���õ������н���ʵ�����Ƶ�
	*/
	void MatrixTran(const varray<varray<Vec4>> & nurbsSF, const varray<varray<Vec4>>& TranMat,
		varray<varray<varray<Vec4>>>& allNurbsSF);

	//ɨ������Nurbs��
	void SweepSurface(const varray<varray<varray<Vec4>>>& allNurbsSF,
		const varray<varray<varray<double>>>& SFw, const varray<double>& pos);

	//ȡ��ߴ�
	void MaxDegree(const varray<SplineSurface>& surfaces, int& uDegree, int& vDegree);

	//�ڵ�ʸ������
	void KnotsUnify(const varray<SplineSurface>& surfaces, varray<double>& NewUKnots, varray<double>& NewVKnots);
};
}
#endif