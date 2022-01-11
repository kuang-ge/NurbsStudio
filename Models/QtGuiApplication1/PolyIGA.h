#pragma once
#include <algorithm>

#include "SplineSurface.h"
#include "SplineVolume.h"
#include "FittingBSpline.h"
#include "varray.h"
#include "definition.h"
#include "XIOStream.h"
#include "CNurbs.h"

using namespace std;
using namespace base;
using namespace bsl;
namespace base{

	
	class CurveEdge
	{
	public:
		int startidx;
		int endidx;
		int icriticalStart;
		int icriticalEnd;
		varray<int>  PassbyidxinOrder;   //·�������ĵ����š�
		int iStatus;
		int persistence;
		bool bVisited;
		bool delflag;
	public:
		CurveEdge();
		~CurveEdge();
		void Clear() { startidx = -1; endidx = -1; PassbyidxinOrder.clear(); };
		bool HasTooPoints() {
			if (PassbyidxinOrder.size() > 100)
				return true;
			return false;
		};
		bool HasLoop()
		{
			for (int i = 0; i < PassbyidxinOrder.size() - 1; i++)
			{
				for (int j = i + 1; j < PassbyidxinOrder.size(); j++)
				{
					if (PassbyidxinOrder.at(i) == PassbyidxinOrder.at(j))
					{
						return true;
					}
				}
			}
			return false;
		};
	};

	

class Sided4Patch  
{
public:
	Sided4Patch();
	~Sided4Patch();
public:
	varray<CurveEdge> m_curveedges;
	varray<int>       m_includeFaceidx;
	XBaseMesh*        m_pnewMs;
	XBaseMesh*        m_pWholeoriMesh;
	SplineSurface     m_4sidedsplsuf;    //�˴���spline��������ת֮ǰ�ġ�
	bool              m_bHasFitted;
	Vec4              m_geometricCenterPt;
	int               m_twinPatchId;
	bool              m_bEdgeC0ContinuityAdjusted;
	bool              m_bInnerC1ContinuityAdjusted;
private:
	Vec4              m_averageNorm;
	bool              m_bHasComputeAverageNorm;
public:
    void  ConstructBaseMeshFromPatch();
	void  funFitBSplineSurface();
    bool  IstwoPatchHasnorelation(Sided4Patch& sp);
    Vec4  GetAverageNormals();
	void  SetPatchGeometricCenterPt();
private:
	bool  CheckAndMakeValidEdge(varray<varray<Vec4> >& oriedgedata,int& uctrlnum,int& vctrlnum, int degree);
	void  AppendPoints(varray<Vec4>& arr, int ctrlNum);
	//void  OutSurfaceEdgeDataforTest(varray<varray<Vec4> >& edgedata);
};

//ÿ��������Ƭ
class IGAHexCell    
{
public:
	IGAHexCell();
	~IGAHexCell();
public:
	XBaseMesh*        m_WholeoriMesh;
	varray<int>       m_6PatchIdx;  //��άƬ��6��Ƭ��Ƭ��idx���������С�����õ�6��Ƭid��01 for u��23 for v, 45 for w;
	SplineVolume      m_splVol;     //ת��Ϊ�������Ƭ��
    varray<int>       m_faceAdjacentHexIdx;  //�������ڽӵ�����hex�����û�о���Ϊ-1;
    varray<int>       m_controlPtGlobalIDs;  //�洢���Ƶ���ȫ���еı�š���ֵ��0~m_uNum*m_vNum*m_wNum-1
	bool              m_bGenerateGlobalID;   //�Ƿ��Ѿ�����ȫ��ID��
	bool			  m_isorder = false;			//�Ƿ�����
public:
	void  Clear();
	bool  CreateParametricVolume6BoundarySurface(bool efalg=false); //�����Ҫ�����ÿ��patch�Ŀ��Ƶ���˳��
	//Ϊ����splinevolume��׼��
	int   GetIdxIn6PatchofOneHexforNormAccordingtoAxis(int xyzmode, varray<Sided4Patch>& allPatchs);
	bool  OrderAndInput6Patchs(varray<Sided4Patch>& allPatchs,bool uReverse=false,bool vReverse=false,bool wReverse=false);
	void  Out6SurfaceControlPtData(string& filename);
    bool  GetLocalIDByGlobalID(int& uid,int& vid, int& wid,int globalID);
    int   GetGlobalIDByLocalID(int uid, int vid, int wid);
	Vec4  GetControlPtByGlobalID(int globalID);
};

class CPolyParaVolume
{
public:
	CPolyParaVolume(void);
	~CPolyParaVolume(void);
public:
	//void SetMesh(XBaseMesh* pms);
public:
	XBaseMesh*           m_oriWholeMesh;         //��ԭʼ��obj����ģ�͡�
	varray<IGAHexCell>   m_HexVolumes;   //OBJ�����ڰ��������������塣
	varray<Sided4Patch>  m_all4sidePatch;  //OBJ������������ѻ��ֵ��ı���Ƭ��
    bool                 m_showParaVolume;
	int                  m_testInputNum[2];  //���洫�ݲ�����
	float                m_gmaxJoc, m_gminJoc;
	int                  m_gminorJocNum;
	double               m_gIsovalall;

public:
	
	void  ClearData();
	CPolyParaVolume& operator=(const varray<NurbsVol> nurbsVols);
	void  ReadallPachts(string& filename);  //�������ݣ������Ψһ��ڣ������ļ����� 
	void  ReadVolumeTxt(const string filename);
	int   parse_patch_file(const char *file_name, varray<varray<varray<int> > > &patch_edges, varray<varray<int> > &patch_tris); 
	void  OutPatchsData();//���������ƪ��Ƭ������
	void  OutputParaVolumeDataAXL(string filename);
	void  OutputParaVolumeDataTxt(string filename,string Cp);//���volume�������ģ�͵Ĳ�������ʦ�����ļ����������
	void  OutputParaVolumeDataVTK(string path);
	void  OutputConstraintsAndLoadByCtrlID(ofstream& idofs);
	void  OutFittedPatchsData(string filename);
	void  Find6ExtremEndsurfaceControlpointIDs(int surfaceID, varray<int>& CtrlptGlobalIDs);
	void  NearIGAHexCellSameUVW(IGAHexCell *StandarHexCell, IGAHexCell *WaitHexCell,int PatchIndx);//�����е���ģ�ͱ������ٽ���һ��ģ��
	void  Order();//�ڲ�����ָ��һ����ʼ����Ϊ��׼
	
	//��������������.
	void  CreateParametericPatch();
	void  Create2DParametricPatch();
    void  SetTwinPatch();   //��������m_all4sidePatch�е�ֵ��Ϊ��ǰÿ��patchѰ��twinƬ��
	void  SetAdjacentHexIdx();
    void  MakeHexC1Continuity();  //Ŀǰ����hex��8���ǵ㴦��һ������û�п��ǣ����⣬12���ߵ�һ�������Ҳû�п��ǡ�
	void  MakeHexC0Continuity();  //Ŀǰ����hex��8���ǵ㴦û�п��ǡ�
	void  MakeHexContinuity();    //������������
    int   GetHexIdAccordingtoPatchID(int patchid);
public:
	//void  DisplayPolyVolume(int imode);      //��ʾ�����ڲ����������塣
	//void  DisplayBoundary4sidedPatch(int imode);  //��ʾ��������ѻ����ı��档
	bool  GetMirrorIJKOnShareSurface(int volume1, int volume2,int imode,int uid,int vid,int wid,int& muid,int& mvid,int& mwid);
    void  ResetInnerCtrlPoints(int volume1, int volume2,int uid,int vid,int wid,int muid,int mvid,int mwid,Vec4 v0,Vec4 v1,Vec4 v2);
public:
	void  ConstructControlptsIDsforAllBlocks(); //Ϊÿ��volume�е�m_CtrlPtsIDinKKKMatrix���¸�ֵ������ÿ�����Ƶ���װ��������ȥ��
    int   GenerateGlobalIDforOneVolume(int volumeid, int startID);
	void  TestAllblockQuanity();
};
class PolyIGASolution
{
public:
    PolyIGASolution();
	~PolyIGASolution();
public:
	CPolyParaVolume* m_pMultiBlockDomain;
	bool             m_allVolHasSolution;
public:
	void SetMultiVolume(CPolyParaVolume* pBlocks);
    void ReadPatchSolutionData(const char* filename);
    void GenerateSolutionRenderData();
	void GenerateControlpointsDisplacementOrStressColor();  //��һ��Ӧ����λ��ֵ��
	void ChangeRenderDisplaceOrStress(bool disOrStress); //�л�Ӧ����λ�Ƶ���ʾ��
	void ClearSolutionData();

	/*void RenderDisplacementOrStress(int surfaceMode);*/
};
}