
#include "PolyIGA.h"
#include "XFunction.h"
#include "assert.h"
#include <atlconv.h>

namespace base {
CurveEdge::CurveEdge()
	:startidx(-1)
	, endidx(-1)
	, icriticalStart(-1)
	, icriticalEnd(-1)
	, iStatus(-1)
	, persistence(-1)
	, bVisited(false)
	, delflag(false)
{

}
CurveEdge::~CurveEdge()
{
	PassbyidxinOrder.clear();
}
Sided4Patch::Sided4Patch()
{
	m_pWholeoriMesh = NULL;
	m_pnewMs = NULL;
	m_bHasFitted = false;
	m_bHasComputeAverageNorm = false;
	m_twinPatchId = -1;
	m_bEdgeC0ContinuityAdjusted = false;
	m_bInnerC1ContinuityAdjusted = false;
}
Sided4Patch::~Sided4Patch()
{
	if(m_pnewMs!=NULL)
	{
		delete m_pnewMs;
		m_pnewMs = NULL;
	}
	m_curveedges.clear();
	m_includeFaceidx.clear();
	m_pWholeoriMesh = NULL;
}
void  Sided4Patch::SetPatchGeometricCenterPt()
{
	m_geometricCenterPt.x = m_geometricCenterPt.y = m_geometricCenterPt.z = 0;
	if(m_pnewMs == NULL)
		return;
	int vsz = m_pnewMs->GetVSize();
	if(vsz == 0)
		return;
	for(int i=0; i<vsz; i++)
	{
		m_geometricCenterPt += m_pnewMs->GetV(i).Pos();
	}
	m_geometricCenterPt /= vsz;
}
void  Sided4Patch::ConstructBaseMeshFromPatch()
{
	m_pnewMs = new XBaseMesh();
	//ConstructMeshFromFaceidx(m_pWholeoriMesh,m_includeFaceidx,m_pnewMs);
	m_pnewMs->MakeXEdges();
	m_pnewMs->ComputeNormals();
}
void  Sided4Patch::funFitBSplineSurface()
{
	varray<Vec4> allinputdata;
	varray<varray<Vec4> > edgedata;
	for(int j=0; j<m_pnewMs->GetVSize(); j++)
	{
		Vec4& vt = m_pnewMs->GetV(j).Pos();
		allinputdata.push_back(vt);
	}

	int UCtrlPtNum,VCtrlPtNum;
	int UVdegree;
	UCtrlPtNum = VCtrlPtNum = 4;
	UVdegree = 3;
	//�����֮ǰ���ȼ��һ�±߽��ߵĵ㹻����
	edgedata.resize(4);
    CurveEdge& e1 = m_curveedges.at(0);
	for(int j=0; j<e1.PassbyidxinOrder.size(); j++)
	{
		Vec4 vt = m_pWholeoriMesh->GetV(e1.PassbyidxinOrder.at(j)).Pos();
		edgedata.at(0).push_back(vt);
	}
	CurveEdge& e2 = m_curveedges.at(1);
	for(int j=0; j<e2.PassbyidxinOrder.size(); j++)
	{
		Vec4 vt = m_pWholeoriMesh->GetV(e2.PassbyidxinOrder.at(j)).Pos();
		edgedata.at(1).push_back(vt);
	}
	CurveEdge& e3 = m_curveedges.at(2);
	for(int j=0; j<e3.PassbyidxinOrder.size(); j++)
	{
		Vec4 vt = m_pWholeoriMesh->GetV(e3.PassbyidxinOrder.at(j)).Pos();
		edgedata.at(2).push_back(vt);
	}
	CurveEdge& e4 = m_curveedges.at(3);
	for(int j=0; j<e4.PassbyidxinOrder.size(); j++)
	{
		Vec4 vt = m_pWholeoriMesh->GetV(e4.PassbyidxinOrder.at(j)).Pos();
		edgedata.at(3).push_back(vt);
	}
	if(!CheckAndMakeValidEdge(edgedata,UCtrlPtNum,VCtrlPtNum,UVdegree))
		return;
	//////////////////////////////////////////////////////////////////////////
	//OutSurfaceEdgeDataforTest(edgedata);
	//////////////////////////////////////////////////////////////////////////
	FitBSplineSurface  fitInputPatch;
	//fitInputPatch.FittingBsurfaceNew(allinputdata,edgedata,UVdegree,UVdegree,UCtrlPtNum,VCtrlPtNum,2); //���ڱ���ģ�ͣ��������
	fitInputPatch.FittingBsurface(allinputdata,edgedata,UVdegree,UVdegree,UCtrlPtNum,VCtrlPtNum);
	m_4sidedsplsuf.ChangeSurfaceMode(SplineSurface::SURFACEMODE_NONUNI_BSPLINE);
	m_4sidedsplsuf.SetSurfaceDegree(UVdegree,UVdegree);
	m_4sidedsplsuf.AddCtrlPts(fitInputPatch.m_uvCtrlPts);
	m_4sidedsplsuf.AddKnotsVector(fitInputPatch.m_uKnots,true);
	m_4sidedsplsuf.AddKnotsVector(fitInputPatch.m_vKnots,false);
	m_4sidedsplsuf.CreateShowMesh();
	m_bHasFitted = true;
}
//void  Sided4Patch::OutSurfaceEdgeDataforTest(varray<varray<Vec4> >& edgedata)
//{
//	String filename = _T("C:\\Users\\Administrator\\Desktop\\new\\");
//
//	String str1,tempfilename;
//	static int pcount = 0;
//	str1.Format(_T("%d"), pcount++);
//	tempfilename = filename + str1;
//	tempfilename +=  (_T(".txt"));
//	ofstream ofs(tempfilename);
//	ofs << m_pnewMs->GetVSize() << "\n";
//	for(int j=0; j<m_pnewMs->GetVSize(); j++)
//	{
//		Vec4& vt = m_pnewMs->GetV(j).Pos();
//		ofs << "v " << vt.x << " " << vt.y << " "<< vt.z << " " << "\n";
//	}
//	//���´���洢ÿ����Ƭ�������߽硣
//	for(int j=0; j<edgedata.at(0).size();j++)
//	{
//		Vec4& vt = edgedata.at(0).at(j);
//		ofs<<"e1 " << " "<< vt.x << " " << vt.y << " " << vt.z << " " << "\n";
//	}
//	for(int j=0; j<edgedata.at(1).size();j++)
//	{
//		Vec4& vt = edgedata.at(1).at(j);
//		ofs<<"e2 " << " "<< vt.x << " " << vt.y << " " << vt.z << " " << "\n";
//	}
//	for(int j=0; j<edgedata.at(2).size();j++)
//	{
//		Vec4& vt = edgedata.at(2).at(j);
//		ofs<<"e3 " << " "<< vt.x << " " << vt.y << " " << vt.z << " " << "\n";
//	}
//	for(int j=0; j<edgedata.at(3).size();j++)
//	{
//		Vec4& vt = edgedata.at(3).at(j);
//		ofs<<"e4 " << " "<< vt.x << " " << vt.y << " " << vt.z << " " << "\n";
//	}
//}
bool  Sided4Patch::CheckAndMakeValidEdge(varray<varray<Vec4> >& oriedgedata,int& uctrlnum,int& vctrlnum, int degree)
{
	if(oriedgedata.size() != 4)
		return false;
	if(degree < 2)
		return false;
	if(uctrlnum < degree + 1 || vctrlnum < degree +1 )
		return false;
	//���ȸ������ߵĳ����ʵ�������ϵĿ��Ƶ���Ŀ��
	/*int Num1 = oriedgedata.at(0).size();
	int Num2 = oriedgedata.at(1).size();
	int Num3 = oriedgedata.at(2).size();
	int Num4 = oriedgedata.at(3).size();
	float len1 = (oriedgedata.at(0).at(0) - oriedgedata.at(0).at(Num1-1)).Magnitude();
	float len2 = (oriedgedata.at(1).at(0) - oriedgedata.at(1).at(Num2-1)).Magnitude();
	float len3 = (oriedgedata.at(2).at(0) - oriedgedata.at(2).at(Num3-1)).Magnitude();
	float len4 = (oriedgedata.at(3).at(0) - oriedgedata.at(3).at(Num4-1)).Magnitude();
	float len5 = len1 > len2 ? len1 : len2;
	float len6 = len3 > len4 ? len3 : len4;
	float len7 = len5 > len6 ? len5 : len6;
	float len8 = len7 / (uctrlnum < vctrlnum ? uctrlnum : vctrlnum);
	if(len8 > 30)
	{
		uctrlnum *= (int)(len8 / 30);
		vctrlnum *= (int)(len8 / 30);
	}*/
	//������ش��ڸ����Ŀ��Ƶ���
	for(int i=0; i<4; i++)
	{
		int edgeptnum = oriedgedata.at(i).size();
		if(edgeptnum < 2)
			return false;
		if((i==0||i==2)&&edgeptnum < uctrlnum)
		{
			AppendPoints(oriedgedata.at(i),uctrlnum);
		}
		else if((i == 1||i==3) && edgeptnum < vctrlnum)
		{
			AppendPoints(oriedgedata.at(i),vctrlnum);
		}
	}
	//���Ҫ���ڸ�������С������
	int minGivenPointNum = 5;
	for(int i=0; i<4; i++)
	{
		if(oriedgedata.at(i).size() < minGivenPointNum)
		{
             AppendPoints(oriedgedata.at(i),minGivenPointNum);
		}
	}
	return true;
}
void  Sided4Patch::AppendPoints(varray<Vec4>& arr, int ctrlNum)
{
	int ptnumdifference = ctrlNum - arr.size();
	//����㣬ÿ�ζ������
	do   
	{   
		//����Ѱ�ұ߳����ߡ�
		int edgeid = -1;
		float length = -1;
		for(int i=0; i<arr.size()-1 ;i++)
		{
			float templength = (arr.at(i+1) - arr.at(i)).GetLength();
		    if(templength > length)
			{
				length = templength;
				edgeid = i;
			}
		}
		if(edgeid != -1)
		{
			Vec4 vt = (arr.at(edgeid) + arr.at(edgeid+1))/2;
			arr.insert(arr.begin()+edgeid+1,vt);
		}
		ptnumdifference = ctrlNum - arr.size();
	} while (ptnumdifference > -1); //��֤�������ڿ��Ƶ�����1.
}
Vec4  Sided4Patch::GetAverageNormals()
{
	Vec4 vt;
	if(m_bHasComputeAverageNorm)
	{
		vt = m_averageNorm;
	}
	else
	{
		vt.x = vt.y = vt.z = 0.f;
		//m_pnewMs->ComputeNormals();
		for(int i=0; i<m_pnewMs->GetFSize(); i++)
		{
			vt += m_pnewMs->GetF(i).Norm();
		}
		vt = vt / m_pnewMs->GetFSize();
		vt.Normalize();
		m_averageNorm = vt;
		m_bHasComputeAverageNorm = true;
	}
	return vt;
}
bool  Sided4Patch::IstwoPatchHasnorelation(Sided4Patch& sp)  //����Ƭ��ȫ�޵��غϡ�
{
	//ֻ���ж�����patch��Ӧ��mesh����ͬ�㼴�ɡ�����spline��ͬ���жϣ�����Ƚ��鷳������������Ŀ��Ƶ����жϵ�ʱ�򣬿��ܻ������С�
	//int nsize1 = m_pnewMs->GetVSize();
	//int nsize2 = sp.m_pnewMs->GetVSize();
	//for(int i=0; i<nsize1; i++)
	//{
	//	for(int j=0; j<nsize2; j++)
	//	{
	//		 if((m_pnewMs->GetV(i).Pos() -  sp.m_pnewMs->GetV(j).Pos()).Magnitude() < ERRF)  //�����غϱ�׼��
	//			 return false;
	//	}
	//}
	/////////////////////////////////////////////////////
	//����ֱ���ж�����spline�����Ƿ����غϵ㡣
	if(m_4sidedsplsuf.GetUCtrlPtNum()*m_4sidedsplsuf.GetVCtrlPtNum() != sp.m_4sidedsplsuf.GetUCtrlPtNum()*sp.m_4sidedsplsuf.GetVCtrlPtNum())
		return false;
	for(int i=0; i<m_4sidedsplsuf.GetVCtrlPtNum(); i++)
	{
		for(int j=0; j<m_4sidedsplsuf.GetUCtrlPtNum();j++)
		{
			for(int k=0; k<sp.m_4sidedsplsuf.GetVCtrlPtNum(); k++)
			{
				for(int l=0; l<sp.m_4sidedsplsuf.GetUCtrlPtNum();l++)
				{
					if((m_4sidedsplsuf.GetCtrlPt(i,j) - sp.m_4sidedsplsuf.GetCtrlPt(k,l)).Magnitude() < ERRF)
						return false;
				}
			}
		}

	}
	return true;
}
IGAHexCell::IGAHexCell()
{
	m_controlPtGlobalIDs.clear();
	m_bGenerateGlobalID = false;
}
IGAHexCell::~IGAHexCell()
{

}
bool IGAHexCell::GetLocalIDByGlobalID(int& uid,int& vid, int& wid,int globalID)
{
	if(!m_bGenerateGlobalID)
		return false;
	int i;
	uid = vid = wid = -1;
	for(i=0; i < m_controlPtGlobalIDs.size(); i++)
	{
		 if(m_controlPtGlobalIDs.at(i) == globalID)
		 {
			 break;
		 }
	}
	if(i == m_controlPtGlobalIDs.size())
		return false;

	if(!m_splVol.GetIJKIndex(i,uid,vid,wid,m_splVol.m_uNum,m_splVol.m_vNum,m_splVol.m_wNum))
	{
		return false;
	}
	return true;
}
Vec4 IGAHexCell::GetControlPtByGlobalID(int globalID)
{
	int idu,idv,idw;
	GetLocalIDByGlobalID(idu,idv,idw,globalID);
	Vec4 vt = m_splVol.GetControlPoint(idu,idv,idw);
	return vt;
}
int IGAHexCell::GetGlobalIDByLocalID(int uid, int vid, int wid)
{
	if(!m_bGenerateGlobalID)
		return -1;
	int idx = m_splVol.GetControlPointIndex(uid,vid,wid);
	int globalidx = m_controlPtGlobalIDs.at(idx);
	return globalidx;
}
//xyzmode  = -1��1 : Ѱ��x���������Ƭ��
//xyzmode = -2,2�� y��������
//xyzmode = -3,3: z�������
//���º������㷨����ʱ�����ܻ������⡣�޷���ȷ�ж�������
int IGAHexCell::GetIdxIn6PatchofOneHexforNormAccordingtoAxis(int xyzmode,varray<Sided4Patch>& allPatchs)
{
	assert(m_6PatchIdx.size() == 6);  //�ڶ�������ʱ���Ѿ���Ƭ�ı�ŷ��ڴ������С�
	assert(allPatchs.size() >= 6);
    //Ѱ����ӽ�x�Ḻ�������Ϊwidx1�档 Ŀǰֻ���������������
	Vec4 baseNorm;
	if(xyzmode == -1)
	{
		baseNorm = Vec4(-1,0,0);
	}
	else if(xyzmode == 1)
	{
		baseNorm = Vec4(1, 0, 0);
	}
	else if(xyzmode == -2)
	{
		baseNorm = Vec4(0,-1,0);
	}
	else if(xyzmode == 2)
	{
		baseNorm = Vec4(0,1,0);
	}
	else if(xyzmode == -3)
	{
		baseNorm = Vec4(0,0,-1);
	}
	else if(xyzmode == 3)
	{
		baseNorm = Vec4(0,0,1);
	}

	Vec4 pPatchNorm[6];
	for(int i=0; i<6; i++)    
	{
		pPatchNorm[i] = allPatchs.at(m_6PatchIdx.at(i)).GetAverageNormals();
	}
	int maxidx = 0;
	float normangval = pPatchNorm[0].dotmultiple(baseNorm);
	for(int i=1; i<6; i++)
	{
		float val = pPatchNorm[i].dotmultiple(baseNorm);
		if(val >  normangval)
		{
			maxidx = i;
			normangval = val;
		}
	}
	return maxidx;
}
bool  IGAHexCell::OrderAndInput6Patchs(varray<Sided4Patch>& allPatchs,bool uReverse,bool vReverse,bool wReverse)
{
	assert(m_6PatchIdx.size() == 6);  //�ڶ�������ʱ���Ѿ���Ƭ�ı�ŷ��ڴ������С�
	assert(allPatchs.size() >= 6);
	int uidx1,uidx2,vidx1,vidx2,widx1,widx2;
	uidx1 = uidx2 = vidx1 = vidx2 = widx1 = widx2 = -1;

	//���ٷ���������id��ʾ������ǰ������m_6PatchIdx�е�id��
	int idx[6];
	idx[0] = GetIdxIn6PatchofOneHexforNormAccordingtoAxis(-1,allPatchs);
	idx[1] = GetIdxIn6PatchofOneHexforNormAccordingtoAxis(1,allPatchs);
	idx[2] = GetIdxIn6PatchofOneHexforNormAccordingtoAxis(-2,allPatchs);
	idx[3] = GetIdxIn6PatchofOneHexforNormAccordingtoAxis(2,allPatchs);
	idx[4] = GetIdxIn6PatchofOneHexforNormAccordingtoAxis(-3,allPatchs);
	idx[5] = GetIdxIn6PatchofOneHexforNormAccordingtoAxis(3,allPatchs);
	
	bool flag = true;
	for(int i=0; i<6; i++)  //���Ψһ��
	{
		if(flag)
		{
			for(int j=0; j<6; j++)
			{
				if(i == j) continue;
				if(idx[i] == idx[j])
				{
					flag = false;
					break;
				}
			}
		}
	}
	//�����ȷ�ԡ�
	if(flag)
	{
		bool flag1 = allPatchs.at(m_6PatchIdx.at(idx[0])).IstwoPatchHasnorelation(allPatchs.at(m_6PatchIdx.at(idx[1])));
		bool flag2 = allPatchs.at(m_6PatchIdx.at(idx[2])).IstwoPatchHasnorelation(allPatchs.at(m_6PatchIdx.at(idx[3])));
		bool flag3 = allPatchs.at(m_6PatchIdx.at(idx[4])).IstwoPatchHasnorelation(allPatchs.at(m_6PatchIdx.at(idx[5])));
		assert(flag1);
		assert(flag2);
		assert(flag3);
		if(!flag1 || !flag2 || !flag3)
			flag = false;
	}
	
	if(flag)  //�����ڳ�ͻ��id
	{
		uidx1 = m_6PatchIdx.at(idx[0]);
		uidx2 = m_6PatchIdx.at(idx[1]);
		vidx1 = m_6PatchIdx.at(idx[2]);
		vidx2 = m_6PatchIdx.at(idx[3]);
		widx1 = m_6PatchIdx.at(idx[4]);
		widx2 = m_6PatchIdx.at(idx[5]);
	}
	else   //���ڳ�ͻ��id�������Ҫ������һ�ַ��������⡣
	{
		int idx = GetIdxIn6PatchofOneHexforNormAccordingtoAxis(-3,allPatchs);
		widx1 = m_6PatchIdx.at(idx);
		widx2 = vidx1 = vidx2 = uidx1 = uidx2 = -1;
		for(int i=0; i<m_6PatchIdx.size(); i++)
		{
			if(i == widx1) 
				continue;
			if(allPatchs.at(widx1).IstwoPatchHasnorelation(allPatchs.at(m_6PatchIdx.at(i))))
			{
				widx2 = m_6PatchIdx.at(i);
				break;
			}
		}
		assert(widx2 != -1);

		for(int i=0; i<m_6PatchIdx.size(); i++)
		{
			if(m_6PatchIdx.at(i) == widx1 || m_6PatchIdx.at(i) == widx2)
			{
				continue;
			}
			vidx1 = m_6PatchIdx.at(i);
			break;
		}
		assert(vidx1 != -1);

		for(int i=0; i<m_6PatchIdx.size(); i++)
		{
			if(m_6PatchIdx.at(i) == widx1 || m_6PatchIdx.at(i) == widx2 || m_6PatchIdx.at(i) == vidx1)
			{
				continue;
			}
			if(allPatchs.at(vidx1).IstwoPatchHasnorelation(allPatchs.at(m_6PatchIdx.at(i))))
			{
				vidx2 = m_6PatchIdx.at(i);
				break;
			}
		}
		assert(vidx2 != -1);

		for(int i=0; i<m_6PatchIdx.size(); i++)
		{
			if(m_6PatchIdx.at(i) == widx1 || m_6PatchIdx.at(i) == widx2 || m_6PatchIdx.at(i) == vidx1 || m_6PatchIdx.at(i) == vidx2)
			{
				continue;
			}
			uidx1 = m_6PatchIdx.at(i);
			break;
		}
		assert(uidx1 != -1);

		for(int i=0; i<m_6PatchIdx.size(); i++)
		{
			if(m_6PatchIdx.at(i) == widx1 || m_6PatchIdx.at(i) == widx2 || m_6PatchIdx.at(i) == vidx1 || m_6PatchIdx.at(i) == vidx2 || m_6PatchIdx.at(i) == uidx1)
			{
				continue;
			}
			uidx2 = m_6PatchIdx.at(i);
			break;
		}
		assert(uidx2 != -1);

		if(!allPatchs.at(uidx1).IstwoPatchHasnorelation(allPatchs.at(uidx2)))
		{
			assert(true);
			//AfxMessageBox(_T("��������"));

		}
		//����ֻ���г�������id������ÿ��id�������Ǻ��������Ӧ����ö�Ӧһ�¡�widx1��widx2��Z���Ӧ��vidx1��vidx2��Y���Ӧ��uidx1��uidx2��X���Ӧ��
		//����W���Ѿ���Z���Ӧ�����ֻ���ж�v���u��
		float len1 = abs(allPatchs.at(vidx1).GetAverageNormals().dotmultiple(Vec4(0,-1,0)));
		float len2 = abs(allPatchs.at(vidx2).GetAverageNormals().dotmultiple(Vec4(0,1,0)));
		float len3 = abs(allPatchs.at(uidx1).GetAverageNormals().dotmultiple(Vec4(0,-1,0)));
		float len4 = abs(allPatchs.at(uidx2).GetAverageNormals().dotmultiple(Vec4(0,1,0)));
		float len5 = len1 > len2 ? len1 : len2;
		float len6 = len3 > len4 ? len3 : len4;
		if(len5 < len6)
		{
			swap(vidx1,uidx1);
			swap(vidx2,uidx2);
		}
	}

	if(uidx1 == -1 || uidx2 == -1 || vidx1 == -1 || vidx2 == -1 || widx1 == -1 || widx2 == -1)
		return false;

	//////////////////////////////////////////////////////////////////////////
	//���ڷ�������һ��׼ȷ��������ݷ��������㷨��һ�������ԡ�
	//�����һ��ȷ��uidx1��uidx2,�Լ�vidx1��vidx2,widx1��widx2��������
	/*float tempval = allPatchs.at(uidx1).GetAverageNormals().dotmultiple(Vec4(-1,0,0));
	if(tempval)
	{
	swap(uidx1,uidx2);
	}
	tempval = allPatchs.at(uidx2).GetAverageNormals().dotmultiple(Vec4(1,0,0));
	assert(tempval > 0);
	tempval = allPatchs.at(vidx1).GetAverageNormals().dotmultiple(Vec4(0,-1,0));
	if(tempval < 0)
	{
	swap(vidx1,vidx2);
	}
	tempval = allPatchs.at(vidx2).GetAverageNormals().dotmultiple(Vec4(0,1,0));
	assert(tempval > 0);*/
	//���ݰ�Χ�����жϡ�
	//������������ȽϹ�������˿������������������жϴ���λ�á�
	//---------------------------------------------------------------------------------------------------\
	//-----------------------------------------------------------------------------------------------
	//����deckelģ�ͣ���73,80,105,106�ŵ�Ԫ��u���жϴ�������һ��Ҫswap
	if(uReverse || allPatchs.at(uidx1).m_geometricCenterPt.x > allPatchs.at(uidx2).m_geometricCenterPt.x)
	{
		swap(uidx1,uidx2);
	}
	if(vReverse ||allPatchs.at(vidx1).m_geometricCenterPt.y > allPatchs.at(vidx2).m_geometricCenterPt.y)
	{
		swap(vidx1,vidx2);
	}
	if(wReverse || allPatchs.at(widx1).m_geometricCenterPt.z > allPatchs.at(widx2).m_geometricCenterPt.z)
	{
		swap(widx1,widx2);
	}
	//////////////////////////////////////////////////////////////////////////

	varray<int> temp6idx = m_6PatchIdx;
	m_6PatchIdx.clear();
	m_6PatchIdx.resize(6);
	m_6PatchIdx.at(U0sfPos) = uidx1;
	m_6PatchIdx.at(U1sfPos) = uidx2;
	m_6PatchIdx.at(V0sfPos) = vidx1;
	m_6PatchIdx.at(V1sfPos) = vidx2;
	m_6PatchIdx.at(W0sfPos) = widx1;
	m_6PatchIdx.at(W1sfPos) = widx2;
	m_splVol.m_6boundarySurface.clear();
	m_splVol.m_6boundarySurface.resize(6);
	m_splVol.m_6boundarySurface.at(U0sfPos) = allPatchs.at(uidx1).m_4sidedsplsuf;
	m_splVol.m_6boundarySurface.at(U1sfPos) = allPatchs.at(uidx2).m_4sidedsplsuf;
	m_splVol.m_6boundarySurface.at(V0sfPos) = allPatchs.at(vidx1).m_4sidedsplsuf;
	m_splVol.m_6boundarySurface.at(V1sfPos) = allPatchs.at(vidx2).m_4sidedsplsuf;
	m_splVol.m_6boundarySurface.at(W0sfPos) = allPatchs.at(widx1).m_4sidedsplsuf;
	m_splVol.m_6boundarySurface.at(W1sfPos) = allPatchs.at(widx2).m_4sidedsplsuf;
	/*CString filename = _T("C:\\Users\\Administrator\\Desktop\\new\\test_beforerotate.txt");
	Out6SurfaceControlPtData(filename);*/
	return true;
}
bool  IGAHexCell::CreateParametricVolume6BoundarySurface(bool efalg)
{
	//����surface�Ѿ���������
	int istate[5];
	int posidx[5];
	m_splVol.m_6boundarySurface.at(W0sfPos).RotateBaseSurface();
	if(efalg)
	    m_splVol.m_6boundarySurface.at(W0sfPos).RotateControlPoints(3);//��������p33ģ�ͣ�ת�÷�������
	istate[0] = m_splVol.m_6boundarySurface.at(W0sfPos).GetTwoSurfaceEdgeshareStateNew(m_splVol.m_6boundarySurface.at(U0sfPos));
	istate[1] = m_splVol.m_6boundarySurface.at(W0sfPos).GetTwoSurfaceEdgeshareStateNew(m_splVol.m_6boundarySurface.at(U1sfPos));
	istate[2] = m_splVol.m_6boundarySurface.at(W0sfPos).GetTwoSurfaceEdgeshareStateNew(m_splVol.m_6boundarySurface.at(V0sfPos));
	istate[3] = m_splVol.m_6boundarySurface.at(W0sfPos).GetTwoSurfaceEdgeshareStateNew(m_splVol.m_6boundarySurface.at(V1sfPos));

	posidx[0] =  m_splVol.m_6boundarySurface.at(U0sfPos).RotateTheWholeSurfaceAccordingtoBasesurface(istate[0]);
	posidx[1] =  m_splVol.m_6boundarySurface.at(U1sfPos).RotateTheWholeSurfaceAccordingtoBasesurface(istate[1]);
	posidx[2] =  m_splVol.m_6boundarySurface.at(V0sfPos).RotateTheWholeSurfaceAccordingtoBasesurface(istate[2]);	
	posidx[3] =  m_splVol.m_6boundarySurface.at(V1sfPos).RotateTheWholeSurfaceAccordingtoBasesurface(istate[3]); 

	istate[4] = m_splVol.m_6boundarySurface.at(U0sfPos).GetTwoSurfaceEdgeshareStateNew(m_splVol.m_6boundarySurface.at(W1sfPos));
	posidx[4] =  m_splVol.m_6boundarySurface.at(W1sfPos).RotateLastSuface(istate[4]); 

	bool flag = true;
	for(int i=0; i<5; i++)  //���һ�²��ܴ��ڲ�����ıߡ�
	{
		assert(istate[i] != 0);
	    if(istate[i] == 0)
		{
			 flag = false;
		}
	}

	for(int i=0; i<6; i++)
	{
		m_splVol.m_6boundarySurface.at(i).CreateShowMesh(true);
	}

	/*CString filename = _T("C:\\Users\\Administrator\\Desktop\\new\\test_afterrotate.txt");
	Out6SurfaceControlPtData(filename);*/
	for(int j=0; j<6; j++)
	{
		m_splVol.m_6boundarySurface.at(j).Generate4SidedBoudaryCurve();
	}
	return flag;
}
CPolyParaVolume::CPolyParaVolume(void)
{
	m_oriWholeMesh = NULL;
	m_showParaVolume = false;
	m_testInputNum[0] = 0;
	m_testInputNum[1] = 6;
	m_gminorJocNum = 0;
	m_gminJoc = 0;
	m_gmaxJoc = 0;
	m_gIsovalall = 0;
}


CPolyParaVolume::~CPolyParaVolume(void)
{

}
void  CPolyParaVolume::ClearData()
{
	m_oriWholeMesh = NULL;
	m_HexVolumes.clear(); 
	m_all4sidePatch.clear(); 
	m_showParaVolume = false;
}
void  CPolyParaVolume::ReadallPachts(string& filename)
{
	m_HexVolumes.clear(); 
	m_all4sidePatch.clear(); 
	m_showParaVolume = false;

	USES_CONVERSION;//����ת����
	const char* filenm = filename.c_str();
	varray<varray<int> >      facelist;
	varray<varray<varray<int> > >    edgelist;
	parse_patch_file(filenm,edgelist,facelist);
	int patchnum = facelist.size();
	int startid = 0;
	//////////////////////////////////////////////////////////////////////////
	//���ֻ��������Ƭ����ȡ���������仰��ע�͡�
	startid = m_testInputNum[0];
	patchnum = m_testInputNum[1];
	//////////////////////////////////////////////////////////////////////////
	for(int i=startid; i<patchnum; i++)
	{
		if(i < startid)
			continue;
		Sided4Patch pat;
		pat.m_pWholeoriMesh = m_oriWholeMesh;
		pat.m_includeFaceidx.clear();
		//����ÿƬ�ĵ㡣
		for(int j=0; j<facelist.at(i).size(); j++)
		{
			pat.m_includeFaceidx.push_back(facelist.at(i).at(j));
		}
		
		//����ÿƬ�ı�
		pat.m_curveedges.resize(edgelist.at(i).size());
		for(int j=0;j<edgelist.at(i).size(); j++)
		{
			for(int k=0; k<edgelist.at(i).at(j).size();k++)
			{
				pat.m_curveedges.at(j).PassbyidxinOrder.push_back(edgelist.at(i).at(j).at(k));
			}
			pat.m_curveedges.at(j).startidx = pat.m_curveedges.at(j).PassbyidxinOrder.at(0);
			pat.m_curveedges.at(j).endidx = pat.m_curveedges.at(j).PassbyidxinOrder.back();
		}
		m_all4sidePatch.push_back(pat);
	}

	//��ʼ��ÿ��patch
	m_HexVolumes.clear();
	int volumeNum = m_all4sidePatch.size()/6;
	m_HexVolumes.resize(volumeNum);
	for(int i=0; i<volumeNum;  i++)
	{
		m_HexVolumes.at(i).m_WholeoriMesh = m_oriWholeMesh;
		m_HexVolumes.at(i).m_6PatchIdx.clear();
		for(int j=6*i; j<6*(i+1); j++)
		{
			m_HexVolumes.at(i).m_6PatchIdx.push_back(j);
		}
	}
	//����ÿƬ��mesh
	int patchsize = m_all4sidePatch.size();
	for(int i=0; i<patchsize; i++)
	{
		m_all4sidePatch.at(i).ConstructBaseMeshFromPatch();
	}
	//OutPatchsData();
}
void  CPolyParaVolume::ReadVolumeTxt(const string filename)
{
	//��Ŀ���ļ�
	ifstream ifs;
	ifs.open(filename);
	assert(ifs.is_open());

	int i;

	while (!ifs.eof())
	{
		string str;
		ifs >> str;
	contin:
		if (str == "PN")
		{
			int count;
			ifs >> count;
			m_HexVolumes.resize(count);
		}
		if (str == "PI")
		{
			ifs >> i;
		}
		if (str == "OD")
		{
			ifs >> m_HexVolumes.at(i).m_splVol.m_uDegree;
			ifs >> m_HexVolumes.at(i).m_splVol.m_vDegree;
			ifs >> m_HexVolumes.at(i).m_splVol.m_wDegree;
			//�����������Ĵ�������Ϊ��ȡ�ļ���ʱ����㲻��ȷ������ֶ�������ȷ�ˣ����԰������ע��
			m_HexVolumes.at(i).m_splVol.m_uDegree--;
			m_HexVolumes.at(i).m_splVol.m_vDegree--;
			m_HexVolumes.at(i).m_splVol.m_wDegree--;
		}
		if (str == "UK")
		{
			int count;
			ifs >> count;
			m_HexVolumes.at(i).m_splVol.m_uKnots.resize(count);
			for (int j = 0; j < count; j++)
			{
				double temp;
				ifs >> temp;
				m_HexVolumes.at(i).m_splVol.m_uKnots[j] = temp;
			}
		}
		if (str == "VK")
		{
			int count;
			ifs >> count;
			m_HexVolumes.at(i).m_splVol.m_vKnots.resize(count);
			for (int j = 0; j < count; j++)
			{
				double temp;
				ifs >> temp;
				m_HexVolumes.at(i).m_splVol.m_vKnots[j] = temp;
			}
		}
		if (str == "WK")
		{
			int count;
			ifs >> count;
			m_HexVolumes.at(i).m_splVol.m_wKnots.resize(count);
			for (int j = 0; j < count; j++)
			{
				double temp;
				ifs >> temp;
				m_HexVolumes.at(i).m_splVol.m_wKnots[j] = temp;
			}
		}
		if (str == "CP")
		{
			ifs >> m_HexVolumes.at(i).m_splVol.m_uNum;
			ifs >> m_HexVolumes.at(i).m_splVol.m_vNum;
			ifs >> m_HexVolumes.at(i).m_splVol.m_wNum;
			ifs >> str;
			while(str!="PI")
			{
				VolumeVertex v;
				v.m_pt.x = stof(str);
				ifs >> v.m_pt.y;
				ifs >> v.m_pt.z;
				m_HexVolumes.at(i).m_splVol.m_vAllCtrlPts.push_back(v);
				ifs >> str;
				if (ifs.eof())
				{
					return ;
				}
			}
			goto contin;
		}
	}
}

CPolyParaVolume& CPolyParaVolume::operator=(const varray<NurbsVol> nurbsVols)
{
	m_HexVolumes.resize(nurbsVols.size());
	for (int i = 0; i != nurbsVols.size(); i++)
	{
		//u�ڵ�����
		for (int j = 0; j != nurbsVols.at(i).m_uKnots.size(); ++j)
		{
			//m_HexVolumes.at(i).m_splVol.m_uKnots.push_back(nurbsVols.at(i).m_uKnots[j]);
			m_HexVolumes.at(i).m_splVol.m_vKnots.push_back(nurbsVols.at(i).m_uKnots[j]);
		}
		//v�ڵ�����
		for (int j = 0; j != nurbsVols.at(i).m_vKnots.size(); ++j)
		{
			//m_HexVolumes.at(i).m_splVol.m_vKnots.push_back(nurbsVols.at(i).m_vKnots[j]);
			m_HexVolumes.at(i).m_splVol.m_uKnots.push_back(nurbsVols.at(i).m_vKnots[j]);
		}
		//w�ڵ�����
		for (int j = 0; j != nurbsVols.at(i).m_wKnots.size(); ++j)
		{
			m_HexVolumes.at(i).m_splVol.m_wKnots.push_back(nurbsVols.at(i).m_wKnots[j]);
		}
		//uvw����
		/*m_HexVolumes.at(i).m_splVol.m_uDegree = nurbsVols.at(i).m_uDegree;
		m_HexVolumes.at(i).m_splVol.m_vDegree = nurbsVols.at(i).m_vDegree;*/
		m_HexVolumes.at(i).m_splVol.m_vDegree = nurbsVols.at(i).m_uDegree;
		m_HexVolumes.at(i).m_splVol.m_uDegree = nurbsVols.at(i).m_vDegree;
		m_HexVolumes.at(i).m_splVol.m_wDegree = nurbsVols.at(i).m_wDegree;

		/*m_HexVolumes.at(i).m_splVol.m_uNum = nurbsVols.at(i).m_uNum;
		m_HexVolumes.at(i).m_splVol.m_vNum = nurbsVols.at(i).m_vNum;*/
		m_HexVolumes.at(i).m_splVol.m_vNum = nurbsVols.at(i).m_uNum;
		m_HexVolumes.at(i).m_splVol.m_uNum = nurbsVols.at(i).m_vNum;
		m_HexVolumes.at(i).m_splVol.m_wNum = nurbsVols.at(i).m_wNum;

		//���Ƶ�
		for (int j = 0; j != nurbsVols.at(i).m_CtrlPts.size(); ++j)
		{
			VolumeVertex v;
			v.m_pt.x = nurbsVols.at(i).m_CtrlPts[j].x;
			v.m_pt.y = nurbsVols.at(i).m_CtrlPts[j].y;
			v.m_pt.z = nurbsVols.at(i).m_CtrlPts[j].z;
			v.m_pt.w = nurbsVols.at(i).m_CtrlPts[j].w;
			m_HexVolumes.at(i).m_splVol.m_vAllCtrlPts.push_back(v);
		}
	}

	return *this;
}

void  IGAHexCell::Clear()
{
	m_6PatchIdx.clear();
	m_splVol.Clear();
	m_WholeoriMesh = NULL;
	m_faceAdjacentHexIdx.clear();
}
void  IGAHexCell::Out6SurfaceControlPtData(string& filename)
{
	ofstream ofs;
	ofs.open(filename.data(), ios::out, 0);
	for(int i=0; i<m_splVol.m_6boundarySurface.size(); i++)
	{
		for(int j=0; j<m_splVol.m_6boundarySurface.at(i).GetVCtrlPtNum();j++)
		{
			for(int k=0; k<m_splVol.m_6boundarySurface.at(i).GetUCtrlPtNum();k++)
			{
				Vec4 vt = m_splVol.m_6boundarySurface.at(i).GetCtrlPt(j,k);
				ofs << "v " << vt.x << " " << vt.y << " "<< vt.z << " " << "\n";
			}
			ofs <<  "\n";
		}	
		ofs << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" <<  "\n";
	}
}
void  CPolyParaVolume::Create2DParametricPatch()
{
	//�ڶ���Ƭ����ʱ���Ѿ���Ƭ���뵽ÿ��hex�У���˴˴�������ݡ�
	m_HexVolumes.clear();

	//����ÿƬ��mesh
	int patchsize = m_all4sidePatch.size();
	for(int i=0; i<patchsize; i++)
	{
		m_all4sidePatch.at(i).funFitBSplineSurface();  //���
		m_all4sidePatch.at(i).m_4sidedsplsuf.Generate4SidedBoudaryCurve();
		m_all4sidePatch.at(i).m_4sidedsplsuf.CreateShowMesh(true);
	}

	//���Ƚ���twinƬ���ڲ�Ƭɾ�����ڲ�Ƭһ������twinƬ����������û��Ƭ��
	//�������֮���ټ���twinƬ��ʹ����ϼ��������࣬�Ժ��ٸĽ���
	SetTwinPatch();  //�˱�־���ж�����Ƭ����ɾ�����ڲ���Ƭ���ܻ������⡣

	///*���´���ɾ��ʱ�������⡣�����ǵȺ�û�����ء�
	//patchsize = m_all4sidePatch.size(); 
	//for(int i=0; i<patchsize; i++)
	//{
	//if(m_all4sidePatch.at(i).m_twinPatchId != -1)
	//{
	//m_all4sidePatch.erase(m_all4sidePatch.begin()+i);
	//patchsize--;
	//i--;
	//}
	//}*/
	m_showParaVolume = false;
}
void  CPolyParaVolume::CreateParametericPatch()
{
	//����ÿƬ��mesh
	int patchsize = m_all4sidePatch.size();
	for(int i=0; i<patchsize; i++)
	{
	     m_all4sidePatch.at(i).funFitBSplineSurface();  //���
	}
	//OutFittedPatchsData();
	int volumeNum = patchsize/6;   //Ŀǰÿ��patch����6���档
	m_HexVolumes.resize(volumeNum);
	for(int i=0; i<volumeNum; i++)
	{
		 m_HexVolumes.at(i).m_WholeoriMesh = m_oriWholeMesh;
		 bool uflag = false;
		 bool vflag = false;
		 bool wflag = false;
		 bool eflag = false;

		 //////////////////////////////////////////////////////////////////////////
		 //�˴�ȫΪdeckelģ��
		 if(i==73 || i==80 || i==105 || i ==106)   //for deckel ��438-444�������γɵ�block�����û���⣬�ڲ�ֵ�����⡣��Ҫu��ǿ�е�����
			 uflag = true;
		 //////////////////////////////////////////////////////////////////////////
		 //�˴�ȫΪp33����ģ��
		 //if(i==90 || i==91 || i==114|| i==116 || i==132 || i==143||i==144 || i==148) //for ����p33ģ��,��Щ��ŵ�block��v��ȫ������
			// vflag = true;
		 //if(i==115)  //for ����p33ģ�ͣ�w�򷭵���
			// wflag = true; 
		 //if(i==142)  //�˴�����ת������
			// eflag=true;
		 //////////////////////////////////////////////////////////////////////////
		 if(!m_HexVolumes.at(i).OrderAndInput6Patchs(m_all4sidePatch,uflag,vflag,wflag))  //����
			 continue;
         if(!m_HexVolumes.at(i).CreateParametricVolume6BoundarySurface(eflag))  //���ɱ�����������档
		 {
			 //string str;
			 //str.Format(_T("��%d����"),i);
			 //str += _T("���ڲ�����ı�");
			 //AfxMessageBox(str);
		 }
	}
	SetAdjacentHexIdx();
	MakeHexC0Continuity();
	for(int i=0; i<volumeNum; i++)
	{
	     m_HexVolumes.at(i).m_splVol.InterpolateInnerControlPtsByBoundarysueface();   //�����ڲ����Ƶ�
	}
	m_showParaVolume = true;
}
void  CPolyParaVolume::OutFittedPatchsData(string filename)
{
	//CString filename = _T("C:\\Users\\Administrator\\Desktop\\new\\allfittedControlpts.txt");
	ofstream ofs(filename);
	for(int i=0; i<m_all4sidePatch.size(); i++)
	{
		ofs << i << "\n";
		for(int j=0; j<m_all4sidePatch.at(i).m_4sidedsplsuf.GetVCtrlPtNum();j++)
		{
			for(int k=0; k<m_all4sidePatch.at(i).m_4sidedsplsuf.GetUCtrlPtNum();k++)
			{
				Vec4 vt = m_all4sidePatch.at(i).m_4sidedsplsuf.GetCtrlPt(j,k);
				ofs << "v " << vt.x << " " << vt.y << " "<< vt.z << " " << "\n";
			}
		}	
	}
}
void  CPolyParaVolume::OutPatchsData()
{
	//����ļ�·��������ʹ��
	string filename = "C:\\Users\\Administrator\\Desktop\\new\\";
	filename += "allPatchData.txt";

	ofstream ofs(filename);
	ofs << "��Ƭ����" << m_all4sidePatch.size() << "\n";
	for(int i=0; i<m_all4sidePatch.size(); i++)
	{
		//���´��뽫ÿ������洢��һ���ļ��С�
		string str1,tempfilename;
		str1 += to_string(i);
		//str1.Format(_T("%d"), i);
		tempfilename = filename + str1;
		tempfilename +=  ".txt";
		//SaveObj(filename,m_all4sidePatch.at(i).m_pnewMs);
		ofstream ofs(tempfilename);
		ofs << m_all4sidePatch.at(i).m_pnewMs->GetVSize() << "\n";
		for(int j=0; j<m_all4sidePatch.at(i).m_pnewMs->GetVSize(); j++)
		{
			Vec4& vt = m_all4sidePatch.at(i).m_pnewMs->GetV(j).Pos();
			ofs << "v " << vt.x << " " << vt.y << " "<< vt.z << " " << "\n";
		}
		//���´���洢ÿ����Ƭ������obj��ʽ��
		/*for(int j=0; j<m_all4sidePatch.at(i).m_pnewMs->GetFSize(); j++)
		{
			XFace& vf = m_all4sidePatch.at(i).m_pnewMs->GetF(j);
			ofs<<"f " << vf.GetIndex(0) << " " << vf.GetIndex(1) << " " << vf.GetIndex(2) << " " << "\n";
		}*/
		//���´��뽫���е���Ƭ�洢��һ���ļ��С�
		//ofs <</* "��" << i << "Ƭ����Ŀ��" <<*/ m_all4sidePatch.at(i).m_pnewMs->GetVSize() << "\n";
		//for(int j=0; j<m_all4sidePatch.at(i).m_pnewMs->GetVSize(); j++)
		//{
		//	Vec4& vt = m_all4sidePatch.at(i).m_pnewMs->GetV(j).Pos();
		//	ofs << "v " << vt.x << " " << vt.y << " "<< vt.z << " " << "\n";
		//}
		//���´���洢ÿ����Ƭ�������߽硣
        CurveEdge& e1 = m_all4sidePatch.at(i).m_curveedges.at(0);
		for(int j=0; j<e1.PassbyidxinOrder.size(); j++)
		{
			Vec4 vt = m_oriWholeMesh->GetV(e1.PassbyidxinOrder.at(j)).Pos();
			ofs<<"e1 " << " "<< vt.x << " " << vt.y << " " << vt.z << " " << "\n";
		}
		CurveEdge& e2 = m_all4sidePatch.at(i).m_curveedges.at(1);
		for(int j=0; j<e2.PassbyidxinOrder.size(); j++)
		{
			Vec4 vt = m_oriWholeMesh->GetV(e2.PassbyidxinOrder.at(j)).Pos();
			ofs<<"e2 " << " "<< vt.x << " " << vt.y << " " << vt.z << " " << "\n";
		}
		CurveEdge& e3 = m_all4sidePatch.at(i).m_curveedges.at(2);
		for(int j=0; j<e3.PassbyidxinOrder.size(); j++)
		{
			Vec4 vt = m_oriWholeMesh->GetV(e3.PassbyidxinOrder.at(j)).Pos();
			ofs<<"e3 " << " "<< vt.x << " " << vt.y << " " << vt.z << " " << "\n";
		}
		CurveEdge& e4 = m_all4sidePatch.at(i).m_curveedges.at(3);
		for(int j=0; j<e4.PassbyidxinOrder.size(); j++)
		{
			Vec4 vt = m_oriWholeMesh->GetV(e4.PassbyidxinOrder.at(j)).Pos();
			ofs<<"e4 " << " "<< vt.x << " " << vt.y << " " << vt.z << " " << "\n";
		}
	}
}
//surfaceID = 0 for xmin;
//surfaceID = 1 for xmax;
//surfaceID = 2 for ymin;
//surfaceID = 3 for ymax;
//surfaceID = 4 for zmin;
//surfaceID = 5 for zmax;
//CtrlptGlobalIDs �����ظ��Ŀ��Ƶ��š�
void  CPolyParaVolume::Find6ExtremEndsurfaceControlpointIDs(int surfaceID, varray<int>& CtrlptGlobalIDs)
{
	CtrlptGlobalIDs.clear();

	float extremVal = 0;
	int volumeNum = m_HexVolumes.size();
	int volCtrlnum;
	Vec4 vt;
	float valtorlerence = 5.0;
	if(surfaceID == 0 || surfaceID == 1)
	{
		if(surfaceID == 0)
		{
			extremVal = 1.0e6;
			for(int i=0; i<volumeNum; i++)
			{
				IGAHexCell& hexcell = m_HexVolumes.at(i);
				volCtrlnum = hexcell.m_controlPtGlobalIDs.size();
				for(int j=0; j<volCtrlnum; j++)
				{
					vt = hexcell.GetControlPtByGlobalID(hexcell.m_controlPtGlobalIDs.at(j));
					if(vt.x < extremVal)
					{
						extremVal = vt.x;
					}
				}
			}
		}
		else 
		{
			extremVal = -1.0e6;
			for(int i=0; i<volumeNum; i++)
			{
				IGAHexCell& hexcell = m_HexVolumes.at(i);
				volCtrlnum = hexcell.m_controlPtGlobalIDs.size();
				for(int j=0; j<volCtrlnum; j++)
				{
					vt = hexcell.GetControlPtByGlobalID(hexcell.m_controlPtGlobalIDs.at(j));
					if(vt.x > extremVal)
					{
						extremVal = vt.x;
					}
				}
			}
		}
		for(int volid=0; volid<volumeNum; volid++)
		{
			IGAHexCell& hexcell = m_HexVolumes.at(volid);
			for(int k=0; k<hexcell.m_splVol.m_wNum; k++)
			{
				for(int j=0; j<hexcell.m_splVol.m_vNum; j++)
				{
					for(int i=0; i<hexcell.m_splVol.m_uNum; i++)
					{
						if(abs(hexcell.m_splVol.GetControlPoint(i,j,k).x - extremVal) < valtorlerence)
						{
							CtrlptGlobalIDs.push_back(hexcell.GetGlobalIDByLocalID(i,j,k));
						}
					}
				}
			}
		}
	}
	else if(surfaceID == 2 || surfaceID == 3)
	{
		if(surfaceID == 2)
		{
			extremVal = 1.0e6;
			for(int i=0; i<volumeNum; i++)
			{
				IGAHexCell& hexcell = m_HexVolumes.at(i);
				volCtrlnum = hexcell.m_controlPtGlobalIDs.size();
				for(int j=0; j<volCtrlnum; j++)
				{
					vt = hexcell.GetControlPtByGlobalID(hexcell.m_controlPtGlobalIDs.at(j));
					if(vt.y < extremVal)
					{
						extremVal = vt.y;
					}
				}
			}
		}
		else 
		{
			extremVal = -1.0e6;
			for(int i=0; i<volumeNum; i++)
			{
				IGAHexCell& hexcell = m_HexVolumes.at(i);
				volCtrlnum = hexcell.m_controlPtGlobalIDs.size();
				for(int j=0; j<volCtrlnum; j++)
				{
					vt = hexcell.GetControlPtByGlobalID(hexcell.m_controlPtGlobalIDs.at(j));
					if(vt.y > extremVal)
					{
						extremVal = vt.y;
					}
				}
			}
		}
		for(int volid=0; volid<volumeNum; volid++)
		{
			IGAHexCell& hexcell = m_HexVolumes.at(volid);
			for(int k=0; k<hexcell.m_splVol.m_wNum; k++)
			{
				for(int j=0; j<hexcell.m_splVol.m_vNum; j++)
				{
					for(int i=0; i<hexcell.m_splVol.m_uNum; i++)
					{
						if(abs(hexcell.m_splVol.GetControlPoint(i,j,k).y - extremVal) < valtorlerence)
						{
							CtrlptGlobalIDs.push_back(hexcell.GetGlobalIDByLocalID(i,j,k));
						}
					}
				}
			}
		}
	}
	else if(surfaceID == 4 || surfaceID == 5)
	{
		if(surfaceID == 4)
		{
			extremVal = 1.0e6;
			for(int i=0; i<volumeNum; i++)
			{
				IGAHexCell& hexcell = m_HexVolumes.at(i);
				volCtrlnum = hexcell.m_controlPtGlobalIDs.size();
				for(int j=0; j<volCtrlnum; j++)
				{
					vt = hexcell.GetControlPtByGlobalID(hexcell.m_controlPtGlobalIDs.at(j));
					if(vt.z < extremVal)
					{
						extremVal = vt.z;
					}
				}
			}
		}
		else 
		{
			extremVal = -1.0e6;
			for(int i=0; i<volumeNum; i++)
			{
				IGAHexCell& hexcell = m_HexVolumes.at(i);
				volCtrlnum = hexcell.m_controlPtGlobalIDs.size();
				for(int j=0; j<volCtrlnum; j++)
				{
					vt = hexcell.GetControlPtByGlobalID(hexcell.m_controlPtGlobalIDs.at(j));
					if(vt.z > extremVal)
					{
						extremVal = vt.z;
					}
				}
			}
		}
		for(int volid=0; volid<volumeNum; volid++)
		{
			IGAHexCell& hexcell = m_HexVolumes.at(volid);
			for(int k=0; k<hexcell.m_splVol.m_wNum; k++)
			{
				for(int j=0; j<hexcell.m_splVol.m_vNum; j++)
				{
					for(int i=0; i<hexcell.m_splVol.m_uNum; i++)
					{
						if(abs(hexcell.m_splVol.GetControlPoint(i,j,k).z - extremVal) < valtorlerence)
						{
							CtrlptGlobalIDs.push_back(hexcell.GetGlobalIDByLocalID(i,j,k));
						}
					}
				}
			}
		}
	}
	//ɾ���������ظ������ݡ�
	int idCount;
	idCount = CtrlptGlobalIDs.size();
	for(int i=0; i<idCount; i++)
	{
		for(int j=i+1; j<idCount;j++)
		{
			if(CtrlptGlobalIDs.at(i) == CtrlptGlobalIDs.at(j))
			{
				CtrlptGlobalIDs.erase(CtrlptGlobalIDs.begin()+j);
				idCount--;
			}
		}
	}
}
void  CPolyParaVolume::NearIGAHexCellSameUVW(IGAHexCell *StandarHexCell, IGAHexCell *WaitHexCell, int PatchIndx)
{
	Sided4Patch temp;
	vector<vector<int>> continor, continor2;
	vector<int> IndxOrder, IndxOrder_old, IndxOrder_orient, FinalOrder;
	vector<Vec4> tt, tt1, tt2, tt3;
	int PatchIndx_r;
	switch (PatchIndx)
	{
	case 0:
		PatchIndx_r = 1;
		break;
	case 1:
		PatchIndx_r = 0;
		break;
	case 2:
		PatchIndx_r = 3;
		break;
	case 3:
		PatchIndx_r = 2;
		break;
	case 4:
		PatchIndx_r = 5;
		break;
	case 5:
		PatchIndx_r = 4;
		break;
	default:
		break;
	}

	//��άתһά
	auto TwoToOne = [](varray<varray<Vec4>> a)
	{
		vector<Vec4> b;
		for (int i = 0; i != a.size(); ++i)
		{
			for (int j = 0; j != a[i].size(); ++j)
			{
				b.push_back(a[i][j]);
			}
		}
		return b;
	};
	//�ж��������������Ƿ�һ��
	auto isSame = [](vector<Vec4> v1, vector<Vec4> v2)
	{
		for (int i = 0; i != v1.size(); ++i)
		{
			bool a = false;
			for (int j = 0; j != v2.size(); ++j)
			{
				if (abs(v1[i].x - v2[j].x) < 0.00001&&abs(v1[i].y - v2[j].y) < 0.00001&& abs(v1[i].z - v2[j].z) < 0.00001)
				{
					a = true;
					break;
				}
				else {
					continue;
				}
			}
			if (a == false)
			{
				return false;
			}
		}
		return true;
	};

	//����һ�����ҵ���ͬ������ͬ���棬���ر��
	auto findSurface = [=]() {
		varray<varray<Vec4>>a = m_all4sidePatch[StandarHexCell->m_6PatchIdx[PatchIndx]].m_4sidedsplsuf.m_UVCtrlPts;
		vector<Vec4> aa = TwoToOne(a);
		for (int i = 0; i != 6; ++i)
		{
			varray<varray<Vec4>> b = m_all4sidePatch[WaitHexCell->m_6PatchIdx[i]].m_4sidedsplsuf.m_UVCtrlPts;
			vector<Vec4> bb = TwoToOne(b);
			if (bb.size() == aa.size() && isSame(aa, bb))
			{
				return i;
			}
			else
			{
				continue;
			}
		}
		return -1;
	};
	//���Ƶ��Ӧ����ģ���ڱ��
	auto findNumber = [](vector<Vec4> tt, IGAHexCell* cell)
	{
		vector<int> order;
		for (int i = 0; i != tt.size(); ++i)
		{
			for (int j = 0; j != cell->m_splVol.m_vAllCtrlPts.size(); ++j)
			{
				if (abs(tt[i].x - cell->m_splVol.m_vAllCtrlPts[j].m_pt.x) < 0.00001&&abs(tt[i].y - cell->m_splVol.m_vAllCtrlPts[j].m_pt.y) < 0.00001&&abs(tt[i].z - cell->m_splVol.m_vAllCtrlPts[j].m_pt.z) < 0.00001)
				{
					order.push_back(j);
				}
			}
		}
		return order;
	};

	//����һ��
	int waitIndx = findSurface();

	auto GetFinalOrder = [&]()
	{
		//�������
		temp = m_all4sidePatch[StandarHexCell->m_6PatchIdx[PatchIndx]];
		//���ϵĿ��Ƶ�ת��Ϊһά����
		tt1 = TwoToOne(temp.m_4sidedsplsuf.m_UVCtrlPts);
		//�ڴ����������Ѱ�ҵ���Щ���Ƶ�,�õ�һ���ɿ��Ƶ��±�
		IndxOrder_old = findNumber(tt1, WaitHexCell);

		temp = m_all4sidePatch[WaitHexCell->m_6PatchIdx[waitIndx]];
		tt2 = TwoToOne(temp.m_4sidedsplsuf.m_UVCtrlPts);
		IndxOrder_orient = findNumber(tt2, WaitHexCell);

		//finalorder�ļ���
		for (int i = 0; i != IndxOrder_orient.size(); ++i)
		{
			for (int j = 0; j != IndxOrder_old.size(); ++j)
			{
				if (IndxOrder_old[j] == IndxOrder_orient[i])
				{
					FinalOrder.push_back(j);
					break;
				}
			}
		}
	};



	auto MySwitch = [&](int x, IGAHexCell *b)
	{
		vector<vector<int>> c;
		switch (x)
		{

		case 0://uv0����
			for (int i = b->m_splVol.m_wNum - 1; i != -1; i--)
			{
				for (int j = 0; j != b->m_splVol.m_vNum; ++j)
				{
					for (int k = 0; k != b->m_splVol.m_uNum; ++k)
					{
						tt3.push_back(b->m_splVol.GetControlPoint(k, j, i));
					}
				}
				c.push_back(findNumber(tt3, b));
				tt3.clear();
			};
			break;
		case 1://uv1����
			for (int i = 0; i != b->m_splVol.m_wNum; i++)
			{
				for (int j = 0; j != b->m_splVol.m_vNum; ++j)
				{
					for (int k = 0; k != b->m_splVol.m_uNum; ++k)
					{
						tt3.push_back(b->m_splVol.GetControlPoint(k, j, i));
					}
				}
				c.push_back(findNumber(tt3, b));
				tt3.clear();
			};
			break;
		case 2://vw0����
			for (int i = b->m_splVol.m_uNum - 1; i != -1; i--)
			{
				for (int j = 0; j != b->m_splVol.m_wNum; ++j)
				{
					for (int k = 0; k != b->m_splVol.m_vNum; ++k)
					{
						tt3.push_back(b->m_splVol.GetControlPoint(i, k, j));
					}
				}
				c.push_back(findNumber(tt3, b));
				tt3.clear();
			};
			break;
		case 3://vw1����
			for (int i = 0; i != b->m_splVol.m_uNum; i++)
			{
				for (int j = 0; j != b->m_splVol.m_wNum; ++j)
				{
					for (int k = 0; k != b->m_splVol.m_vNum; ++k)
					{
						tt3.push_back(b->m_splVol.GetControlPoint(i, k, j));
					}
				}
				c.push_back(findNumber(tt3, b));
				tt3.clear();
			};
			break;
		case 4://uw0����
			for (int i = b->m_splVol.m_vNum - 1; i != -1; i--)
			{
				for (int j = 0; j != b->m_splVol.m_wNum; ++j)
				{
					for (int k = 0; k != b->m_splVol.m_uNum; ++k)
					{
						tt3.push_back(b->m_splVol.GetControlPoint(k, i, j));
					}
				}
				c.push_back(findNumber(tt3, b));
				tt3.clear();
			};
			break;
		case 5://uw1����
			for (int i = 0; i != b->m_splVol.m_vNum; i++)
			{
				for (int j = 0; j != b->m_splVol.m_wNum; ++j)
				{
					for (int k = 0; k != b->m_splVol.m_uNum; ++k)
					{
						tt3.push_back(b->m_splVol.GetControlPoint(k, i, j));
					}
				}
				c.push_back(findNumber(tt3, b));
				tt3.clear();
			};
			break;
		default:
			break;
		}
		return c;
	};



	//����һ��
	GetFinalOrder();

	continor = MySwitch(PatchIndx_r, StandarHexCell);
	continor2 = MySwitch(waitIndx, WaitHexCell);

	//Ȼ�����finalOrder����һ������
	for (int i = 0; i != continor2.size(); ++i)
	{
		vector<int> c;
		c.resize(continor2[i].size());
		for (int j = 0; j != continor2[i].size(); ++j)
		{
			c[FinalOrder[j]] = continor2[i][j];
			//c.push_back(continor2[i][FinalOrder[j]]);
		}
		continor2[i] = c;
	}

	//���������ģ�ͽ�������
	varray<VolumeVertex> replace;
	replace.resize(WaitHexCell->m_splVol.m_vAllCtrlPts.size());
	for (int i = 0; i != continor.size(); ++i)
	{
		for (int j = 0; j != continor[i].size(); ++j)
		{
			replace[continor[i][j]] = WaitHexCell->m_splVol.m_vAllCtrlPts[continor2[i][j]];
		}
	}
	//����滻
	WaitHexCell->m_splVol.m_vAllCtrlPts = replace;

	//����滻֮����Ҫ��ԭ�ȵ����ݽ��и���

	//����ȡ����������
	vector<Sided4Patch> temp_Sided4Patch;

	//��һ����
	temp.m_4sidedsplsuf.m_UVCtrlPts.clear();
	temp.m_4sidedsplsuf.m_vKnots = WaitHexCell->m_splVol.m_vKnots;//�ڵ�ʸ��
	temp.m_4sidedsplsuf.m_uKnots = WaitHexCell->m_splVol.m_uKnots;
	temp.m_4sidedsplsuf.m_uNum = WaitHexCell->m_splVol.m_uNum;
	temp.m_4sidedsplsuf.m_vNum = WaitHexCell->m_splVol.m_vNum;
	temp.m_4sidedsplsuf.SetSurfaceDegree(WaitHexCell->m_splVol.m_uDegree, WaitHexCell->m_splVol.m_vDegree);//����
	//��һ����Ŀ��Ƶ�
	for (int j = 0; j != WaitHexCell->m_splVol.m_vNum; ++j)
	{
		varray<Vec4> t;
		for (int k = 0; k != WaitHexCell->m_splVol.m_uNum; ++k)
		{
			Vec3 aa = WaitHexCell->m_splVol.GetControlPoint(k, j, 0);
			t.push_back(WaitHexCell->m_splVol.GetControlPoint(k, j, 0));
		}
		temp.m_4sidedsplsuf.m_UVCtrlPts.push_back(t);
	}
	temp_Sided4Patch.push_back(temp);
	temp.m_4sidedsplsuf.m_UVCtrlPts.clear();

	for (int j = 0; j != WaitHexCell->m_splVol.m_vNum; ++j)
	{
		varray<Vec4> t;
		for (int k = 0; k != WaitHexCell->m_splVol.m_uNum; ++k)
		{
			Vec3 tt = WaitHexCell->m_splVol.GetControlPoint(k, j, WaitHexCell->m_splVol.m_wNum - 1);
			t.push_back(WaitHexCell->m_splVol.GetControlPoint(k, j, WaitHexCell->m_splVol.m_wNum - 1));
		}
		temp.m_4sidedsplsuf.m_UVCtrlPts.push_back(t);
	}
	temp_Sided4Patch.push_back(temp);
	temp.m_4sidedsplsuf.m_UVCtrlPts.clear();
	//��������
	temp.m_4sidedsplsuf.m_uKnots = WaitHexCell->m_splVol.m_vKnots;
	temp.m_4sidedsplsuf.m_vKnots = WaitHexCell->m_splVol.m_wKnots;
	temp.m_4sidedsplsuf.m_uNum = WaitHexCell->m_splVol.m_vNum;
	temp.m_4sidedsplsuf.m_vNum = WaitHexCell->m_splVol.m_wNum;
	temp.m_4sidedsplsuf.SetSurfaceDegree(WaitHexCell->m_splVol.m_vDegree, WaitHexCell->m_splVol.m_wDegree);
	//��������Ŀ��Ƶ�
	for (int j = 0; j != WaitHexCell->m_splVol.m_wNum; ++j)
	{
		varray<Vec4> t;
		for (int k = 0; k != WaitHexCell->m_splVol.m_vNum; ++k)
		{
			t.push_back(WaitHexCell->m_splVol.GetControlPoint(0, k, j));
		}
		temp.m_4sidedsplsuf.m_UVCtrlPts.push_back(t);
	}
	temp_Sided4Patch.push_back(temp);
	temp.m_4sidedsplsuf.m_UVCtrlPts.clear();
	//���ĸ���
	for (int j = 0; j != WaitHexCell->m_splVol.m_wNum; ++j)
	{
		varray<Vec4> t;
		for (int k = 0; k != WaitHexCell->m_splVol.m_vNum; ++k)
		{
			t.push_back(WaitHexCell->m_splVol.GetControlPoint(WaitHexCell->m_splVol.m_uNum - 1, k, j));
		}
		temp.m_4sidedsplsuf.m_UVCtrlPts.push_back(t);
	}
	temp_Sided4Patch.push_back(temp);
	temp.m_4sidedsplsuf.m_UVCtrlPts.clear();

	//�������
	temp.m_4sidedsplsuf.m_uKnots = WaitHexCell->m_splVol.m_uKnots;
	temp.m_4sidedsplsuf.m_vKnots = WaitHexCell->m_splVol.m_wKnots;
	temp.m_4sidedsplsuf.m_uNum = WaitHexCell->m_splVol.m_uNum;
	temp.m_4sidedsplsuf.m_vNum = WaitHexCell->m_splVol.m_wNum;
	temp.m_4sidedsplsuf.SetSurfaceDegree(WaitHexCell->m_splVol.m_uDegree, WaitHexCell->m_splVol.m_wDegree);
	for (int j = 0; j != WaitHexCell->m_splVol.m_wNum; ++j)
	{
		varray<Vec4> t;
		for (int k = 0; k != WaitHexCell->m_splVol.m_uNum; ++k)
		{
			t.push_back(WaitHexCell->m_splVol.GetControlPoint(k, 0, j));
		}
		temp.m_4sidedsplsuf.m_UVCtrlPts.push_back(t);
	}
	temp_Sided4Patch.push_back(temp);
	temp.m_4sidedsplsuf.m_UVCtrlPts.clear();
	//��������
	for (int j = 0; j != WaitHexCell->m_splVol.m_wNum; ++j)
	{
		varray<Vec4> t;
		for (int k = 0; k != WaitHexCell->m_splVol.m_uNum; ++k)
		{
			t.push_back(WaitHexCell->m_splVol.GetControlPoint(k, WaitHexCell->m_splVol.m_vNum - 1, j));
		}
		temp.m_4sidedsplsuf.m_UVCtrlPts.push_back(t);
	}
	temp_Sided4Patch.push_back(temp);

	//�ҵ�ƥ���������
	vector<int> tt5;
	for (int i = 0; i != WaitHexCell->m_6PatchIdx.size(); ++i)
	{
		for (int j = 0; j != temp_Sided4Patch.size(); ++j)
		{
			vector<Vec4> a = TwoToOne(m_all4sidePatch[WaitHexCell->m_6PatchIdx[i]].m_4sidedsplsuf.m_UVCtrlPts);
			vector<Vec4> b = TwoToOne(temp_Sided4Patch[j].m_4sidedsplsuf.m_UVCtrlPts);
			if (isSame(a, b))
			{
				tt5.push_back(j);
				break;
			}
		}
	}

	//�������ҵ�֮���ԭ��������������滻��
	for (int i = 0; i != 6; ++i)
	{
		m_all4sidePatch[WaitHexCell->m_6PatchIdx[i]] = temp_Sided4Patch[i];
	}

	//���ݷ���������������
	varray<int> replace_t;
	replace_t.resize(6);
	for (int i = 0; i != WaitHexCell->m_faceAdjacentHexIdx.size(); ++i)
	{
		replace_t[tt5[i]] = WaitHexCell->m_faceAdjacentHexIdx[i];
	}
	WaitHexCell->m_faceAdjacentHexIdx = replace_t;
}

void CPolyParaVolume::Order()
{
	//ǰ�����������
	std::list<IGAHexCell*> List;

	//����
	/*for (int i = 0; i != m_HexVolumes.size(); ++i)
	{
		List.push_back(&m_HexVolumes[i]);
	}*/

	m_HexVolumes[0].m_isorder = true;
	List.push_back(&m_HexVolumes[0]);
	while (!List.empty())
	{
		IGAHexCell* a = List.front();
		List.pop_front();
		for (int i = 0; i != a->m_faceAdjacentHexIdx.size(); ++i)
		{
			int m = a->m_faceAdjacentHexIdx[i];
			if (a->m_faceAdjacentHexIdx[i] != -1 && m_HexVolumes[a->m_faceAdjacentHexIdx[i]].m_isorder != true)
			{
				NearIGAHexCellSameUVW(a, &m_HexVolumes[a->m_faceAdjacentHexIdx[i]], i);
				m_HexVolumes[a->m_faceAdjacentHexIdx[i]].m_isorder = true;
				List.push_back(&m_HexVolumes[a->m_faceAdjacentHexIdx[i]]);
			}
		}
	}
}

void  CPolyParaVolume::OutputParaVolumeDataVTK(string path)
{
	ofstream file(path); //�������ʽ���ļ�
	//file.open("C:\\r\\hexs.txt");
	//
	int uSegnum,vSegnum,wSegnum,patchPtNumCount,patchNumCount;
	uSegnum = vSegnum = wSegnum = 10;
	patchPtNumCount = (uSegnum+1)*(vSegnum+1)*(wSegnum+1);  //�����Ŀ
	int patchNum = m_HexVolumes.size();
	Vec4 pt,mpt;
	int started = 0;

	file<<"# vtk DataFile Version 2.0"<<"\n"<<"TET"<<"\n"<<"ASCII"<<"\n";
	file<<"\n"<<"DATASET UNSTRUCTURED_GRID"<<"\n";
	file<<"POINTS "<< patchNum*patchPtNumCount <<" float"<<"\n";

	for(int patID=0; patID < patchNum; patID++)
	{
		SplineVolume& sv = m_HexVolumes.at(patID).m_splVol;
		for(int k=0; k<=wSegnum;k++)
		{
			for(int j=0; j<=vSegnum; j++)
			{
				for(int i=0; i<=uSegnum; i++)
				{
					pt = sv.GetPoint(i*1.f/uSegnum,j*1.f/vSegnum,k*1.f/wSegnum,mpt);
					file << pt.x << " " << pt.y << " " << pt.z << " " << "\n";
				}
			}
		}
	}

	patchNumCount = uSegnum*vSegnum*wSegnum;   //��Ԫ��Ŀ��
	file << "CELLS "<<patchNum * patchNumCount <<" "<< 9* patchNum * patchNumCount <<"\n";
	for(int patID=0; patID < patchNum; patID++)
	{
		SplineVolume& sv = m_HexVolumes.at(patID).m_splVol;
		started = patID * patchPtNumCount;
		for(int k=0; k<wSegnum;k++)
		{
			for(int j=0; j<vSegnum; j++)
			{
				for(int i=0; i<uSegnum; i++)
				{
					file<<"8 "<<started + k*(uSegnum+1)*(vSegnum+1) + (uSegnum+1)*j+i<< " "
						      <<started + k*(uSegnum+1)*(vSegnum+1) + (uSegnum+1)*j+i+1<<" "
						      <<started + k*(uSegnum+1)*(vSegnum+1) + (uSegnum+1)*(j+1)+i+1<<" "
						      <<started + k*(uSegnum+1)*(vSegnum+1) + (uSegnum+1)*(j+1)+i<<" "
						      <<started + (k+1)*(uSegnum+1)*(vSegnum+1) + (uSegnum+1)*j+i<< " "
						      <<started + (k+1)*(uSegnum+1)*(vSegnum+1) + (uSegnum+1)*j+i+1<<" "
						      <<started + (k+1)*(uSegnum+1)*(vSegnum+1) + (uSegnum+1)*(j+1)+i+1<<" "
						      <<started + (k+1)*(uSegnum+1)*(vSegnum+1) + (uSegnum+1)*(j+1)+i<<"\n";
				}
			}
		}
	}
	
	file << "CELL_TYPES "<<patchNum * patchNumCount<<"\n";
	for (int l=0;l<patchNum * patchNumCount;l++)
	{
		file<<"12"<<"\n";
	}
	file.close();
}
void  CPolyParaVolume::OutputParaVolumeDataTxt
(string filename,string Cp)
{
	//CString filename = _T("C:\\Users\\Administrator\\Desktop\\new\\Controlpoints.txt");
	ofstream ofs(filename);
	ofs << "PN " <<" " << m_HexVolumes.size() << "\n";  //patch number
	for(int i=0; i<m_HexVolumes.size(); i++)
	{
		ofs << "PI" << " " << i << "\n";   //��ǰƬ��id�š�
		
		SplineVolume& sv = m_HexVolumes.at(i).m_splVol;
		ofs << "OD" <<" " << sv.m_uDegree+1 << " "<< sv.m_vDegree+1 << " "<< sv.m_wDegree+1 <<"\n"; //Order

		ofs << "UK" << " " << sv.m_uKnots.size() << "\n";
		for(int i=0; i<sv.m_uKnots.size(); i++)
			ofs << sv.m_uKnots.at(i) << " ";   //�ڵ�����
		ofs <<  "\n";
		ofs << "VK" << " " << sv.m_vKnots.size() << "\n";
		for(int i=0; i<sv.m_vKnots.size(); i++)
			ofs << sv.m_vKnots.at(i) << " ";   //�ڵ�����
		ofs <<  "\n";
		ofs << "WK" << " " << sv.m_wKnots.size() << "\n";
		for(int i=0; i<sv.m_wKnots.size(); i++)
			ofs << sv.m_wKnots.at(i) << " ";   //�ڵ�����
		ofs <<  "\n";

		varray<Vec4> ptvarry;
		sv.GetAllControlPointsVec4(ptvarry);
		ofs << "CP"<<" " << sv.m_uNum << " "<< sv.m_vNum << " "<< sv.m_wNum  << "\n"; //control point number
		for(int i=0; i<ptvarry.size(); i++)
			ofs << ptvarry.at(i).x << " "<< ptvarry.at(i).y << " "<< ptvarry.at(i).z <<"\n";
	}
	ofs.close();
	//���Ƶ������
	//CString filename1 = _T("C:\\Users\\Administrator\\Desktop\\new\\ControlpointsGlobalIDs.txt");
	ofstream idofs(Cp);
	ofs << "PN " <<" " << m_HexVolumes.size() << "\n";  //patch number
	for(int i=0; i<m_HexVolumes.size(); i++)
	{
		int idxnum = m_HexVolumes.at(i).m_controlPtGlobalIDs.size();
		idofs << "PN" <<" "<< i << "\n";
		idofs << "PI" <<" "<< idxnum << "\n";
		for(int j=0; j<idxnum; j++)
		{
			idofs << m_HexVolumes.at(i).m_controlPtGlobalIDs.at(j) << " ";
		}
		idofs << "\n";
	}
	//�����Լ�����������Ŀ��Ƶ��š�
	OutputConstraintsAndLoadByCtrlID(idofs);
	idofs.close();
}
void  CPolyParaVolume::OutputConstraintsAndLoadByCtrlID(ofstream& idofs)
{
	//for sphere
	if(m_HexVolumes.size() == 1)
	{
		idofs << "WC" << " "; //��Լ���ġ��µ���
		for(int j=0; j<m_HexVolumes.at(0).m_splVol.m_vNum; j++)
		{
			for(int i=0; i<m_HexVolumes.at(0).m_splVol.m_uNum;i++)
			{
				int idx = m_HexVolumes.at(0).GetGlobalIDByLocalID(i,j,0);
				idofs << idx << " ";
			}
		}
		idofs << "\n";
		idofs << "WF" << " "; //���������ġ��ϵ���
		for(int j=0; j<m_HexVolumes.at(0).m_splVol.m_vNum; j++)
		{
			for(int i=0; i<m_HexVolumes.at(0).m_splVol.m_uNum;i++)
			{
				int idx = m_HexVolumes.at(0).GetGlobalIDByLocalID(i,j,m_HexVolumes.at(0).m_splVol.m_wNum-1);
				idofs << idx << " ";
			}
		}
		idofs << "\n";
	}


	//for L-shape
	if(m_HexVolumes.size() == 3)
	{
		idofs << "WC" << " "; //��Լ���ġ� �����
		for(int k=0; k<m_HexVolumes.at(2).m_splVol.m_wNum; k++)
		{
			for(int j=0; j<m_HexVolumes.at(2).m_splVol.m_vNum;j++)
			{
				int idx = m_HexVolumes.at(2).GetGlobalIDByLocalID(0,j,k);
				idofs << idx << " ";
			}
		}
		idofs << "\n";
		idofs << "WF" << " "; //���������ġ�  y���϶���
		for(int k=0; k<m_HexVolumes.at(0).m_splVol.m_wNum; k++)
		{
			for(int i=0; i<m_HexVolumes.at(0).m_splVol.m_uNum;i++)
			{
				int idx = m_HexVolumes.at(0).GetGlobalIDByLocalID(i,m_HexVolumes.at(0).m_splVol.m_vNum-1,k);
				idofs << idx << " ";
			}
		}
		idofs << "\n";
	}


	//for torus
	if(m_HexVolumes.size() == 8)
	{
		idofs << "WC" << " "; //��Լ���ġ�
		//m_HexVolumes.at(1) �µ��档
		for(int j=0; j<m_HexVolumes.at(1).m_splVol.m_vNum; j++)
		{
			for(int i=0; i<m_HexVolumes.at(1).m_splVol.m_uNum;i++)
			{
				int idx = m_HexVolumes.at(1).GetGlobalIDByLocalID(i,j,0);
				idofs << idx << " ";
			}
		}
		idofs << "\n";
		idofs << "WF" << " "; //���������ġ�
		//m_HexVolumes.at(6) �ϵ��档
		for(int j=0; j<m_HexVolumes.at(6).m_splVol.m_vNum; j++)
		{
			for(int i=0; i<m_HexVolumes.at(6).m_splVol.m_uNum;i++)
			{
				int idx = m_HexVolumes.at(6).GetGlobalIDByLocalID(i,j,m_HexVolumes.at(6).m_splVol.m_wNum-1);
				idofs << idx << " ";
			}
		}
		idofs << "\n";
	}

	//for deckel
	if(m_HexVolumes.size() == 220)
	{
		idofs << "WC" << " "; //��Լ���ġ�
		//m_HexVolumes.at(6)  //�����    //�����ң���Ҫд�������ҡ�
		varray<int> ConstraintsIDs;
		Find6ExtremEndsurfaceControlpointIDs(0,ConstraintsIDs);
		for(int ii=0; ii<ConstraintsIDs.size();ii++)
		{
			idofs << ConstraintsIDs.at(ii) << " ";
		}
		idofs << "\n";

		idofs << "WF" << " "; //���������ġ�
		//m_HexVolumes.at(126) //�������档
		/*for(int k=0; k<m_HexVolumes.at(126).m_splVol.m_wNum; k++)
		{
		for(int j=0; j<m_HexVolumes.at(126).m_splVol.m_vNum;j++)
		{
		int idx = m_HexVolumes.at(126).GetGlobalIDByLocalID(m_HexVolumes.at(126).m_splVol.m_uNum-1,j,k);
		idofs << idx << " ";
		}
		}*/
		varray<int> loadIDs;
		Find6ExtremEndsurfaceControlpointIDs(1,loadIDs);
		for(int jj=0; jj<loadIDs.size(); jj++)
		{
			idofs << loadIDs.at(jj) << " ";
		}
		idofs << "\n";
	}

	//for block
	if(m_HexVolumes.size() == 24)
	{
		//m_HexVolumes.at(8)  //�µ���
		idofs << "WC" << " "; //��Լ���ġ�
		for(int j=0; j<m_HexVolumes.at(8).m_splVol.m_vNum; j++)
		{
			for(int i=0; i<m_HexVolumes.at(8).m_splVol.m_uNum;i++)
			{
				int idx = m_HexVolumes.at(8).GetGlobalIDByLocalID(i,j,0);
				idofs << idx << " ";
			}
		}
		/*varray<int> ConstraintsIDs;
		Find6ExtremEndsurfaceControlpointIDs(4,ConstraintsIDs);
		for(int ii=0; ii<ConstraintsIDs.size();ii++)
		{
		idofs << ConstraintsIDs.at(ii) << " ";
		}*/

		idofs << "\n";

		idofs << "WF" << " "; //���������ġ�
		//m_HexVolumes.at(15)  //�ϵ���
		for(int j=0; j<m_HexVolumes.at(15).m_splVol.m_vNum; j++)
		{
			for(int i=0; i<m_HexVolumes.at(15).m_splVol.m_uNum;i++)
			{
				int idx = m_HexVolumes.at(15).GetGlobalIDByLocalID(i,j,m_HexVolumes.at(15).m_splVol.m_wNum-1);
				idofs << idx << " ";
			}
		}
		/*varray<int> loadIDs;
		Find6ExtremEndsurfaceControlpointIDs(5,loadIDs);
		for(int jj=0; jj<loadIDs.size(); jj++)
		{
			idofs << loadIDs.at(jj) << " ";
		}*/
		idofs << "\n";
	}

}
void  CPolyParaVolume::OutputParaVolumeDataAXL(string filename)
{
	//CString filename = _T("C:\\Users\\Administrator\\Desktop\\new\\controlpoints.axl");
	ofstream ofs(filename);
	ofs << "patch number= " << m_HexVolumes.size() << "\n";
	ofs << "<axl>" <<  "\n";
	for(int i=0; i<m_HexVolumes.size(); i++)
	{
		string str1,tempfilename;
		str1 += to_string(i);
		//str1.Format(_T("%d"), i);
		tempfilename = "bv_" + str1;
		ofs << "<volume type=" <<  "bspline" << " " << "name=" <<  /*(LPCTSTR)*/tempfilename << ">"<<"\n";
		SplineVolume& sv = m_HexVolumes.at(i).m_splVol;
		ofs << "<number>" << sv.m_uNum << " "<< sv.m_vNum << " "<< sv.m_wNum << "</number>" << "\n";
		ofs << "<order>" << sv.m_uDegree+1 << " "<< sv.m_vDegree+1 << " "<< sv.m_wDegree+1 << "</order>" << "\n";
		ofs << "<knots>";
		for(int i=0; i<sv.m_uKnots.size(); i++)
		    ofs << sv.m_uKnots.at(i)<<" ";
		ofs << "</knots>" << "\n";
		ofs << "<knots>";
		for(int i=0; i<sv.m_vKnots.size(); i++)
			ofs << sv.m_vKnots.at(i)<<" ";
		ofs << "</knots>" << "\n";
		ofs << "<knots>";
		for(int i=0; i<sv.m_wKnots.size(); i++)
			ofs << sv.m_wKnots.at(i)<<" ";
		ofs << "</knots>" << "\n";
		varray<Vec4> ptvarry;
		sv.GetAllControlPointsVec4(ptvarry);
		ofs << "<points>";
		for(int i=0; i<ptvarry.size(); i++)
		    ofs << ptvarry.at(i).x << " "<< ptvarry.at(i).y << " "<< ptvarry.at(i).z << "\n";
		ofs << "</points>" << "\n";
		ofs << "</volume>" << "\n";
	}
	ofs << "<script run="<< "true"<< "></script>" << "\n";
	ofs << "<camera>" << "\n";
	ofs << "<Parameters fieldOfView=" << "0.785398" << "Type="<< "PERSPECTIVE" << "zNearCoefficient="<< "0.005" << "orthoCoef="<< "0.414214" << "zClippingCoefficient="<< "1.73205"<< "/>" << "\n";
	ofs << "<Stereo distToScreen="<< "0.5" << "focusDistance="<< "48.2843" << "physScreenWidth="<< "0.4" << "IODist="<< "0.062"<< "/>" << "\n";
	ofs << "<ManipulatedCameraFrame>" << "\n";
	ofs << "<position x="<< "-0.74886" << "y="<< "0.548172" << "z="<< "1.17488"<< "/>" << "\n";
	ofs << "<orientation q0="<< "0.00386867" << "q1="<< "-0.314504" << "q2="<< "-0.275119" << "q3="<< "0.908505"<< "/>" << "\n";
	ofs << "<ManipulatedParameters transSens="<< "1" << "rotSens="<< "1" << "wheelSens="<< "1" << "spinSens="<< "<< "<< "0.3"<< "/>" << "\n";
	ofs << "<ManipulatedCameraParameters flySpeed="<< "0.2"<< ">" << "\n";
	ofs << "<flyUpVector x="<< "0" << "y="<< "1" << "z="<< "0"<< "/>" << "\n";
	ofs << "</ManipulatedCameraParameters>" << "\n";
	ofs << "</ManipulatedCameraFrame>" << "\n";
	ofs << "</camera>" << "\n";
	ofs << "</axl>" << "\n";
}
int CPolyParaVolume::parse_patch_file(const char *file_name,
	varray<varray<varray<int> > >& patch_edges,
	varray<varray<int> > & patch_tris) 
{
	ifstream ifs;
	ifs.open(file_name);
	assert(ifs.is_open());
	if (ifs.fail())
		return -1;
	int patch_num;
	ifs >> patch_num;
	patch_edges.resize(patch_num);
	patch_tris.resize(patch_num);
	for (int i = 0; i < patch_num; ++i)
	{
		patch_edges[i].resize(4);
		for (int j = 0; j < 4; ++j) 
		{
			int edge_pt_num;
			ifs >> edge_pt_num;
			patch_edges[i][j].resize(edge_pt_num);
			for (int k = 0; k < edge_pt_num; ++k)
				ifs >> patch_edges[i][j][k];
		}
		int patch_tri_num;
		ifs >> patch_tri_num;
		patch_tris[i].resize(patch_tri_num);
		for (int j = 0; j < patch_tri_num; ++j)
			ifs >> patch_tris[i][j];
	}
	if (ifs.fail())
		return -1;
	return 0;
}
//void  CPolyParaVolume::DisplayPolyVolume(int imode)
//{
//	//��Ҫ��ʾÿƬ�������߽���ڲ���Ƭ��
//	if(m_showParaVolume)
//	{
//		int volumeNum = m_HexVolumes.size();
//		for(int i=0; i<volumeNum; i++)
//		{
//			if(m_HexVolumes.at(i).m_splVol.m_bHasIGASolution)  //����еȼ��η���ֵ���Ͳ��ڴ˴���ʾ������������ʾ��
//				continue;
//			if(m_HexVolumes.at(i).m_splVol.m_bHasInterpolated)
//			{
//				 m_HexVolumes.at(i).m_splVol.RenderbyMode(iRenderBoundary_Patch,imode);
//				 m_HexVolumes.at(i).m_splVol.RenderbyMode(iRenderBoundary_ControlPoints,imode);
//				 //m_HexVolumes.at(i).m_splVol.RenderbyMode(iRenderVolume_WithISOElements,imode);
//			}
//			else
//			{
//                  for(int j=0; j<6; j++)
//			      m_HexVolumes.at(i).m_splVol.m_6boundarySurface.at(j).DisplaySplineSurface(imode,true/*WIREFRAME_MODE*/);
//			}
//		}
//	}
//	else
//	{
//        if(m_all4sidePatch.size() > 0)
//	         DisplayBoundary4sidedPatch(imode);
//	}
//};
//void  CPolyParaVolume::DisplayBoundary4sidedPatch(int imode)
//{
//	//��ʾ4���ߡ�
//	if(m_oriWholeMesh == NULL)
//		return;
//	float oldColor[4] = {0,0,0,0};	
//
//	glGetFloatv(GL_CURRENT_COLOR, oldColor);
//
//	int npatch = m_all4sidePatch.size();
//	if(npatch == 0)
//		return;
//
//	glColor4fv(oldColor);
//
//	for(int i=0; i<npatch;i++)
//	{
//		if(m_all4sidePatch.at(i).m_bHasFitted)
//		{	
//			if(m_all4sidePatch.at(i).m_twinPatchId != -1)  //���ڴ�ʱֻ�Ǳ�ʾ��������棬����ڲ����治��Ҫ��ʾ��
//				continue;
//			m_all4sidePatch.at(i).m_4sidedsplsuf.DisplaySplineSurface(imode,true); //��ʾ���������ڵ��档
//			//m_all4sidePatch.at(i).m_4sidedsplsuf.DrawSurfaceCtrlPts(RGB(255,0,0));
//		}
//		else
//		{
//			if(i % 6 == 0)
//			    m_all4sidePatch.at(i).m_pnewMs->RenderMeshWire(COLORREF(RGB(255,0,0)),1);
//			else if(i % 6 == 1)
//				m_all4sidePatch.at(i).m_pnewMs->RenderMeshWire(COLORREF(RGB(0,255,0)),1);
//			else if(i % 6 == 2)
//				m_all4sidePatch.at(i).m_pnewMs->RenderMeshWire(COLORREF(RGB(0,0,255)),1);
//			else if(i % 6 == 3)
//				m_all4sidePatch.at(i).m_pnewMs->RenderMeshWire(COLORREF(RGB(255,255,0)),1);
//			else if(i % 6 == 4)
//				m_all4sidePatch.at(i).m_pnewMs->RenderMeshWire(COLORREF(RGB(0,255,255)),1);
//			else if(i % 6 == 5)
//				m_all4sidePatch.at(i).m_pnewMs->RenderMeshWire(COLORREF(RGB(255,0,255)),1);
//			//��ʾ�����߽硣
//			for(int j =0; j<4; j++)
//			{
//				const CurveEdge& e = m_all4sidePatch.at(i).m_curveedges.at(j);
//				int psize = e.PassbyidxinOrder.size();
//				glColor3f(1.f,0.f,0.f);
//				glBegin(GL_LINES);	
//				for(int k=0; k<psize-1; k++)
//				{
//					const Vec4& pt1 = m_oriWholeMesh->GetV(e.PassbyidxinOrder.at(k)).Pos();
//					const Vec4& pt2 = m_oriWholeMesh->GetV(e.PassbyidxinOrder.at(k+1)).Pos();
//					glVertex3f(pt1.x, pt1.y, pt1.z);
//					glVertex3f(pt2.x, pt2.y, pt2.z);	
//				}
//				glEnd();		
//			}
//		}
//	}
//};
//��������m_all4sidePatch�е�ֵ��Ϊ��ǰÿ��patchѰ��twinƬ��
//patch��洢������ת֮ǰ��splinesurface
//��ÿ��hexcell��洢����ת֮���splinesurface��
void  CPolyParaVolume::SetTwinPatch()
{
	int patchnum = m_all4sidePatch.size();
	if(patchnum < 2)
		return;
	for(int i=0; i<patchnum; i++)
	{
		m_all4sidePatch.at(i).m_twinPatchId = -1;
	}
	varray<int> idxarr;
	for(int i=0; i<patchnum-1; i++)
	{
		Sided4Patch& sp1 = m_all4sidePatch.at(i);
		if(sp1.m_twinPatchId != -1 || !sp1.m_bHasFitted) //�������֮���������twinƬ��
			continue;
		idxarr.clear();
		for(int j=i+1; j<patchnum; j++)
		{
			Sided4Patch& sp2 = m_all4sidePatch.at(j);
			if(sp2.m_twinPatchId != -1)
				continue;
			if(sp1.m_4sidedsplsuf.IstwoSurfaceSame(sp2.m_4sidedsplsuf))
			{
				idxarr.push_back(j);
			}
		}
		if(idxarr.size() < 1) 
			continue;
		int minidx = idxarr.at(0);
		if(idxarr.size() > 1)  //�ҵ����治ֹһ����Ŀǰ��������ġ�����ҵ�����Ƭ��ֹһ����������Ծ�����С���档
		{
			float minlen = sp1.m_4sidedsplsuf.GetDistanceBetweenTwoSurface(m_all4sidePatch.at(idxarr.at(0)).m_4sidedsplsuf);
			for(int k=1; k<idxarr.size();k++)
			{
				float templen = sp1.m_4sidedsplsuf.GetDistanceBetweenTwoSurface(m_all4sidePatch.at(idxarr.at(k)).m_4sidedsplsuf);
				if(templen < minlen)
				{
					minidx = idxarr.at(k);
					minlen = templen;
				}
			}
		}
		sp1.m_twinPatchId = minidx;
		m_all4sidePatch.at(minidx).m_twinPatchId = i;
	}
}
//Ϊÿ��hexcell�ж�6����������
//�������m_faceAdjacentHexIdx
void  CPolyParaVolume::SetAdjacentHexIdx()
{
	SetTwinPatch();
	int adjacentVolID;
	for(int i=0; i<m_HexVolumes.size(); i++)
	{
		m_HexVolumes.at(i).m_faceAdjacentHexIdx.clear();
		//����patch��id�ǰ���m_6PatchIdx�е�˳�����е�
		for(int j=0; j<6; j++)
		{
			int idx = m_all4sidePatch.at(m_HexVolumes.at(i).m_6PatchIdx.at(j)).m_twinPatchId;
			if(idx == -1)
				m_HexVolumes.at(i).m_faceAdjacentHexIdx.push_back(-1);
			else
			{
				adjacentVolID = GetHexIdAccordingtoPatchID(idx);
				m_HexVolumes.at(i).m_faceAdjacentHexIdx.push_back(adjacentVolID);
			}
		}
	}
}
int  CPolyParaVolume::GetHexIdAccordingtoPatchID(int patchid)
{
	int hexcellID = -1;
	if(patchid < 0 || patchid >= m_all4sidePatch.size())
		return -1;
	for(int i=0; i<m_HexVolumes.size();i++)
	{
		if(IsInTheVarray(patchid,m_HexVolumes.at(i).m_6PatchIdx))
		{
			hexcellID = i;
		}
	}
	return hexcellID;
}
//������ʹ��֮ǰ������ȷ��������Ƭ�Ѿ��ź���������ת��λ����������ȷ��
void  CPolyParaVolume::MakeHexC0Continuity()
{
	int volumeNum = m_HexVolumes.size();
	for(int i=0; i<volumeNum; i++)
	{
		for(int j=0; j<m_HexVolumes.at(i).m_faceAdjacentHexIdx.size();j++)
		{
			Sided4Patch& sf1 = m_all4sidePatch.at(m_HexVolumes.at(i).m_6PatchIdx.at(j));
			SplineSurface& ssf1 = m_HexVolumes.at(i).m_splVol.m_6boundarySurface.at(j);
			if(m_HexVolumes.at(i).m_faceAdjacentHexIdx.at(j) == -1 || sf1.m_bEdgeC0ContinuityAdjusted )
				continue;
			IGAHexCell& sp = m_HexVolumes.at(m_HexVolumes.at(i).m_faceAdjacentHexIdx.at(j));
			int secondPatchID;
			if(j == U0sfPos)
				secondPatchID = U1sfPos;
			else if(j == U1sfPos)
				secondPatchID = U0sfPos;
			else if(j == V0sfPos)
				secondPatchID = V1sfPos;
			else if(j == V1sfPos)
				secondPatchID = V0sfPos;
			else if(j == W0sfPos)
				secondPatchID = W1sfPos;
			else if(j == W1sfPos)
				secondPatchID = W0sfPos;
			Sided4Patch& sf2 = m_all4sidePatch.at(sp.m_6PatchIdx.at(secondPatchID));
			SplineSurface& ssf2 = sp.m_splVol.m_6boundarySurface.at(secondPatchID);
			if(sf2.m_bEdgeC0ContinuityAdjusted)
				continue;
			//��sf1��sf2�����п��Ƶ��غϡ�
			//���ڽǵ�ĵط��漰����Ƭ�Ƚ϶࣬��ˣ��ô����Ƶ����Ҳ����������ʱȷ�����е����߱����ڸô��غϡ�
			int uctrlNum = ssf1.GetUCtrlPtNum();
			int vctrlNum = ssf1.GetVCtrlPtNum();
			assert(uctrlNum == ssf2.GetUCtrlPtNum());
			assert(vctrlNum == ssf2.GetVCtrlPtNum());
			for(int k=0;k<vctrlNum;k++)
			{
				for(int l=0; l<uctrlNum;l++)
				{
					if((k==0 && l==0) || (k==0 && l == uctrlNum-1) || (k==vctrlNum-1 && l == 0) ||(k==vctrlNum-1 && l == uctrlNum-1))
						continue;
					Vec4& tempVec1 = ssf1.GetCtrlPt(k,l);
					Vec4& tempVec2 = ssf2.GetCtrlPt(k,l);
					Vec4 vt =  (tempVec1+tempVec2)/2;
					ssf1.SetCtrlPt(k,l,vt);
					ssf2.SetCtrlPt(k,l,vt);
				}
			}
			ssf1.CreateShowMesh(true);
			ssf2.CreateShowMesh(true);
			sf1.m_bEdgeC0ContinuityAdjusted = true;
			sf2.m_bEdgeC0ContinuityAdjusted = true;
		}
	}
}
void  CPolyParaVolume::MakeHexC1Continuity()
{
	//��������ÿƬ���ڽ�id�š�
	int v,f;
	int mk,mj,mi;
	Vec4 vt0,vt1,vt2;
	for(v=0; v<m_HexVolumes.size();v++)
	{
		SplineVolume& sp1 = m_HexVolumes.at(v).m_splVol;
		for(f=0; f<m_HexVolumes.at(v).m_faceAdjacentHexIdx.size();f++)
		{
			int idx = m_HexVolumes.at(v).m_faceAdjacentHexIdx.at(f);
			if(idx == -1)
				continue;
			SplineVolume& sp2 = m_HexVolumes.at(idx).m_splVol;

			//for the u����
			if(f == U0sfPos && !m_all4sidePatch.at(m_HexVolumes.at(v).m_6PatchIdx.at(U0sfPos)).m_bInnerC1ContinuityAdjusted)
			{
				for(int i=1; i<sp1.m_wNum-1;i++)
				{
					for(int j=1; j<sp1.m_vNum-1;j++)
					{
						vt1 = sp1.GetControlPoint(1,j,i);
						vt0 = sp1.GetControlPoint(0,j,i);
						if(!GetMirrorIJKOnShareSurface(v,idx,0,1,j,i,mk,mj,mi))
							continue;
						vt2 = sp2.GetControlPoint(mk,mj,mi);
						ResetInnerCtrlPoints(v,idx,1,j,i,mk,mj,mi,vt0,vt1,vt2);
					}
				}
				m_all4sidePatch.at(m_HexVolumes.at(v).m_6PatchIdx.at(U0sfPos)).m_bInnerC1ContinuityAdjusted = true;
				m_all4sidePatch.at(m_HexVolumes.at(idx).m_6PatchIdx.at(U1sfPos)).m_bInnerC1ContinuityAdjusted = true;
			}
			
			//for the u������
			else if(f == U1sfPos && !m_all4sidePatch.at(m_HexVolumes.at(v).m_6PatchIdx.at(U1sfPos)).m_bInnerC1ContinuityAdjusted)
			{
				for(int i=1; i<sp1.m_wNum-1;i++)
				{
					for(int j=1; j<sp1.m_vNum-1;j++)
					{
						vt1 = sp1.GetControlPoint(sp1.m_uNum-2,j,i);
						vt0 = sp1.GetControlPoint(sp1.m_uNum-1,j,i);
						if(!GetMirrorIJKOnShareSurface(v,idx,0,sp1.m_uNum-2,j,i,mk,mj,mi))
							continue;
						vt2 = sp2.GetControlPoint(mk,mj,mi);
						ResetInnerCtrlPoints(v,idx,sp1.m_uNum-2,j,i,mk,mj,mi,vt0,vt1,vt2);
					}
				}
				m_all4sidePatch.at(m_HexVolumes.at(v).m_6PatchIdx.at(U1sfPos)).m_bInnerC1ContinuityAdjusted = true;
				m_all4sidePatch.at(m_HexVolumes.at(idx).m_6PatchIdx.at(U0sfPos)).m_bInnerC1ContinuityAdjusted = true;
			}
			//for the v����
			else if(f == V0sfPos && !m_all4sidePatch.at(m_HexVolumes.at(v).m_6PatchIdx.at(V0sfPos)).m_bInnerC1ContinuityAdjusted)
			{
				for(int i=1; i<sp1.m_wNum-1;i++)
				{
					for(int k=1; k<sp1.m_uNum-1;k++)
					{
						vt1 = sp1.GetControlPoint(k,1,i);
						vt0 = sp1.GetControlPoint(k,0,i);
						if(!GetMirrorIJKOnShareSurface(v,idx,1,k,1,i,mk,mj,mi))
							continue;
						vt2 = sp2.GetControlPoint(mk,mj,mi);
						ResetInnerCtrlPoints(v,idx,k,1,i,mk,mj,mi,vt0,vt1,vt2);
					}
				}
				m_all4sidePatch.at(m_HexVolumes.at(v).m_6PatchIdx.at(V0sfPos)).m_bInnerC1ContinuityAdjusted = true;
				m_all4sidePatch.at(m_HexVolumes.at(idx).m_6PatchIdx.at(V1sfPos)).m_bInnerC1ContinuityAdjusted = true;
			}
			//for the v������
			else if(f == V1sfPos && !m_all4sidePatch.at(m_HexVolumes.at(v).m_6PatchIdx.at(V1sfPos)).m_bInnerC1ContinuityAdjusted)
			{
				for(int i=1; i<sp1.m_wNum-1;i++)
				{
					for(int k=1; k<sp1.m_uNum-1;k++)
					{
						vt1 = sp1.GetControlPoint(k,sp1.m_vNum-2,i);
						vt0 = sp1.GetControlPoint(k,sp1.m_vNum-1,i);
						if(!GetMirrorIJKOnShareSurface(v,idx,1,k,sp1.m_vNum-2,i,mk,mj,mi))
							continue;
						vt2 = sp2.GetControlPoint(mk,mj,mi);
						ResetInnerCtrlPoints(v,idx,k,sp1.m_vNum-2,i,mk,mj,mi,vt0,vt1,vt2);
					}
				}
				m_all4sidePatch.at(m_HexVolumes.at(v).m_6PatchIdx.at(V1sfPos)).m_bInnerC1ContinuityAdjusted = true;
				m_all4sidePatch.at(m_HexVolumes.at(idx).m_6PatchIdx.at(V0sfPos)).m_bInnerC1ContinuityAdjusted = true;
			}
			//for the w����
			else if(f == W0sfPos && !m_all4sidePatch.at(m_HexVolumes.at(v).m_6PatchIdx.at(W0sfPos)).m_bInnerC1ContinuityAdjusted)
			{
				for(int j=1; j<sp1.m_vNum-1;j++)
				{
					for(int k=1; k<sp1.m_uNum-1;k++)
					{
						vt1 = sp1.GetControlPoint(k,j,1);
						vt0 = sp1.GetControlPoint(k,j,0);
						if(!GetMirrorIJKOnShareSurface(v,idx,2,k,j,1,mk,mj,mi))
							continue;
						vt2 = sp2.GetControlPoint(mk,mj,mi);
						ResetInnerCtrlPoints(v,idx,k,j,1,mk,mj,mi,vt0,vt1,vt2);
					}
				}
				m_all4sidePatch.at(m_HexVolumes.at(v).m_6PatchIdx.at(W0sfPos)).m_bInnerC1ContinuityAdjusted = true;
				m_all4sidePatch.at(m_HexVolumes.at(idx).m_6PatchIdx.at(W1sfPos)).m_bInnerC1ContinuityAdjusted = true;
			}
			//for the w������
			else if(f == W1sfPos && !m_all4sidePatch.at(m_HexVolumes.at(v).m_6PatchIdx.at(W1sfPos)).m_bInnerC1ContinuityAdjusted)
			{
				for(int j=1; j<sp1.m_vNum-1;j++)
				{
					for(int k=1; k<sp1.m_uNum-1;k++)
					{
						vt1 = sp1.GetControlPoint(k,j,sp1.m_wNum-2);
						vt0 = sp1.GetControlPoint(k,j,sp1.m_wNum-1);
						if(!GetMirrorIJKOnShareSurface(v,idx,2,k,j,sp1.m_wNum-2,mk,mj,mi))
							continue;
						vt2 = sp2.GetControlPoint(mk,mj,mi);
						ResetInnerCtrlPoints(v,idx,k,j,sp1.m_wNum-2,mk,mj,mi,vt0,vt1,vt2);
					}
				}
				m_all4sidePatch.at(m_HexVolumes.at(v).m_6PatchIdx.at(W1sfPos)).m_bInnerC1ContinuityAdjusted = true;
				m_all4sidePatch.at(m_HexVolumes.at(idx).m_6PatchIdx.at(W0sfPos)).m_bInnerC1ContinuityAdjusted = true;
			}
		}
		
	}	
}
void  CPolyParaVolume::ResetInnerCtrlPoints(int volume1, int volume2,int uid,int vid,int wid,int muid,int mvid,int mwid,Vec4 v0,Vec4 v1,Vec4 v2)
{
	//���ȵõ��µ�v1��v2����ֵ��
	Vec4 vdir = v2 - v1;
	vdir = vdir.Normalize();
	v2 = (v2 - v0).Magnitude()*vdir + v0;
	vdir *= -1.f;
	v1= (v1 - v0).Magnitude()*vdir + v0;

	//Ȼ��ֱ����á�
	int idx = m_HexVolumes.at(volume1).m_splVol.GetControlPointIndex(uid,vid,wid);
	m_HexVolumes.at(volume1).m_splVol.m_vAllCtrlPts.at(idx) = v1;
	idx = m_HexVolumes.at(volume2).m_splVol.GetControlPointIndex(muid,mvid,mwid);
	m_HexVolumes.at(volume2).m_splVol.m_vAllCtrlPts.at(idx) = v2;
}
//imode �Ǿ����ģʽ����imode=0��������u��Գƣ�����imode =1 ��������v��Գƣ�����imode = 2��������w��Գ�
//uid,vid,wid ��volume1�еı�š���Ҫ��ȡ����volume2�ж�Ӧ�ı�š�
bool  CPolyParaVolume::GetMirrorIJKOnShareSurface(int volume1, int volume2, int imode,int uid,int vid,int wid,int& muid,int& mvid,int& mwid)
{
	muid = mvid = mwid = -1;
	if(volume1 < 0 || volume1 >= m_HexVolumes.size() || volume2 < 0 || volume2 >= m_HexVolumes.size())
		return false;
	muid = uid;
	mvid = vid;
	mwid = wid;
	if(imode == 0)  //��������u��Գơ�
	{
		muid = m_HexVolumes.at(volume2).m_splVol.m_uNum - uid - 1;
	}
	else if(imode == 1)//��������v��Գơ�
	{
		mvid = m_HexVolumes.at(volume2).m_splVol.m_vNum - vid - 1;
	}
	else if(imode == 2)//��������w��Գơ�
	{
		mwid = m_HexVolumes.at(volume2).m_splVol.m_wNum - wid - 1;
	}
	
	return (muid >= 0 && mvid >= 0 && mwid >= 0 && muid < m_HexVolumes.at(volume2).m_splVol.m_uNum && mvid < m_HexVolumes.at(volume2).m_splVol.m_vNum && mwid < m_HexVolumes.at(volume2).m_splVol.m_wNum);
}
void  CPolyParaVolume::MakeHexContinuity()
{
	MakeHexC1Continuity();
	int volumeNum = m_HexVolumes.size();
	for(int i=0; i<volumeNum; i++)
	{
		m_HexVolumes.at(i).m_splVol.ReCreateAllVolumeRenderPts();
	}
}
void  CPolyParaVolume::TestAllblockQuanity()
{
	m_gmaxJoc = -1.0e6;
	m_gminJoc = 1.0e6;
	m_gminorJocNum = 0;
	m_gIsovalall = 0;

	int volnum = m_HexVolumes.size();
	for(int i=0; i<volnum; i++)
	{
		SplineVolume& sv = m_HexVolumes.at(i).m_splVol;
		sv.TestMeshQuanity(10);

		if (sv.m_minJocbian < 0) {
			cout << "��" << i << "Ƭ������" << "\n";
			cout << "Jocbian is " << sv.m_minJocbian << "\n";
		}

		if(sv.m_maxJocbian > m_gmaxJoc)
			m_gmaxJoc = sv.m_maxJocbian;

		if (sv.m_minJocbian < m_gminJoc)
		{
			m_gminJoc = sv.m_minJocbian;
		}
		m_gminorJocNum += sv.m_minorJocbianNum;
		m_gIsovalall += sv.m_isovalAdd;
	}

	cout << "����ſɱ�" << m_gmaxJoc << "\n";
	cout << "��С�ſɱ�" << m_gminJoc << "\n";
	cout << "�ſɱ�Ϊ������" << m_gminorJocNum << "\n";

	//�öԻ������ʽ�����
	/*CString ostr,str;
	ostr.Format(_T("�����ſ˱�ֵ��%f;"),m_gmaxJoc);
	str.Format(_T("��С���ſ˱�ֵ��%f;"),m_gminJoc);
	ostr += str;
	str.Format(_T("�ſ˱�ֵΪ�����ĸ�����%d;"),m_gminorJocNum);
	ostr += str;
	str.Format(_T("�Ȳ�ֵ�ĺͣ�%f;"),m_gIsovalall);
	AfxMessageBox(ostr);*/
}
void CPolyParaVolume::ConstructControlptsIDsforAllBlocks()
{
	int volumeNum = m_HexVolumes.size();
	int startIdx = 0;
	for(int i=0; i<volumeNum; i++)
	{
		int idx = GenerateGlobalIDforOneVolume(i,startIdx);
		startIdx = idx;
	}
};
//startID�ǵ�ǰ�������µı�ſ�������ʼ��ŵ�ֵ�� �ڶԵ�ǰvolume��ֵ����Ժ󣬷���һ�������ŵ�ֵ��
int  CPolyParaVolume::GenerateGlobalIDforOneVolume(int volumeid, int startID)
{
	if(startID < 0 || volumeid < 0 || volumeid >= m_HexVolumes.size())
		return 0;
	//���ȳ�ʼ�����Ƶ��Ÿ�ֵ��'
	IGAHexCell& vol = m_HexVolumes.at(volumeid);
	int unum = vol.m_splVol.m_uNum;
	int vnum = vol.m_splVol.m_vNum;
	int wnum = vol.m_splVol.m_wNum;

	//���ȸ���ֵ��
	int maxsize = unum*vnum*wnum;
	vol.m_controlPtGlobalIDs.resize(maxsize);
	for(int i=0; i<maxsize; i++)
	{
		vol.m_controlPtGlobalIDs.at(i) = -1;
	}
	int pcount = 0;
	int idx = 0;

	int twinPatchID[6],adjacenHexID[6];

	twinPatchID[U0sfPos] = m_all4sidePatch.at(vol.m_6PatchIdx.at(U0sfPos)).m_twinPatchId;
	twinPatchID[U1sfPos] = m_all4sidePatch.at(vol.m_6PatchIdx.at(U1sfPos)).m_twinPatchId;
	twinPatchID[V0sfPos] = m_all4sidePatch.at(vol.m_6PatchIdx.at(V0sfPos)).m_twinPatchId;
	twinPatchID[V1sfPos] = m_all4sidePatch.at(vol.m_6PatchIdx.at(V1sfPos)).m_twinPatchId;
	twinPatchID[W0sfPos] = m_all4sidePatch.at(vol.m_6PatchIdx.at(W0sfPos)).m_twinPatchId;
	twinPatchID[W1sfPos] = m_all4sidePatch.at(vol.m_6PatchIdx.at(W1sfPos)).m_twinPatchId;

	adjacenHexID[U0sfPos] = vol.m_faceAdjacentHexIdx.at(U0sfPos);
	adjacenHexID[U1sfPos] = vol.m_faceAdjacentHexIdx.at(U1sfPos);
	adjacenHexID[V0sfPos] = vol.m_faceAdjacentHexIdx.at(V0sfPos);
	adjacenHexID[V1sfPos] = vol.m_faceAdjacentHexIdx.at(V1sfPos);
	adjacenHexID[W0sfPos] = vol.m_faceAdjacentHexIdx.at(W0sfPos);
	adjacenHexID[W1sfPos] = vol.m_faceAdjacentHexIdx.at(W1sfPos);

	for(int k=0; k<wnum; k++)
	{
		for(int j=0; j<vnum; j++)
		{
			for(int i=0; i<unum; i++)
			{
				//���´����Ƚ���������������
				if(i == 0)
				{
					if(adjacenHexID[U0sfPos] != -1 && m_HexVolumes.at(adjacenHexID[U0sfPos]).m_bGenerateGlobalID)
					{
						IGAHexCell& adjcentBlock = m_HexVolumes.at(adjacenHexID[U0sfPos]);
						int muidx,mvidx,mwidx;
						GetMirrorIJKOnShareSurface(volumeid,adjacenHexID[U0sfPos],0,i,j,k,muidx,mvidx,mwidx);
						int adidx = m_HexVolumes.at(adjacenHexID[U0sfPos]).GetGlobalIDByLocalID(muidx,mvidx,mwidx);
						vol.m_controlPtGlobalIDs.at(idx++) = adidx;
						continue;
					}
				}
				else if(i == unum-1)
				{
					if(adjacenHexID[U1sfPos] != -1 && m_HexVolumes.at(adjacenHexID[U1sfPos]).m_bGenerateGlobalID)
					{
						IGAHexCell& adjcentBlock = m_HexVolumes.at(adjacenHexID[U1sfPos]);
						int muidx,mvidx,mwidx;
						GetMirrorIJKOnShareSurface(volumeid,adjacenHexID[U1sfPos],0,i,j,k,muidx,mvidx,mwidx);
						int adidx = m_HexVolumes.at(adjacenHexID[U1sfPos]).GetGlobalIDByLocalID(muidx,mvidx,mwidx);
						vol.m_controlPtGlobalIDs.at(idx++) = adidx;
						continue;
					}
				}
				if(j == 0)
				{
					if(adjacenHexID[V0sfPos] != -1 && m_HexVolumes.at(adjacenHexID[V0sfPos]).m_bGenerateGlobalID)
					{
						IGAHexCell& adjcentBlock = m_HexVolumes.at(adjacenHexID[V0sfPos]);
						int muidx,mvidx,mwidx;
						GetMirrorIJKOnShareSurface(volumeid,adjacenHexID[V0sfPos],1,i,j,k,muidx,mvidx,mwidx);
						int adidx = m_HexVolumes.at(adjacenHexID[V0sfPos]).GetGlobalIDByLocalID(muidx,mvidx,mwidx);
						vol.m_controlPtGlobalIDs.at(idx++) = adidx;
						continue;
					}
				}
				else if(j == vnum-1)
				{
					if(adjacenHexID[V1sfPos] != -1 && m_HexVolumes.at(adjacenHexID[V1sfPos]).m_bGenerateGlobalID)
					{
						IGAHexCell& adjcentBlock = m_HexVolumes.at(adjacenHexID[V1sfPos]);
						int muidx,mvidx,mwidx;
						GetMirrorIJKOnShareSurface(volumeid,adjacenHexID[V1sfPos],1,i,j,k,muidx,mvidx,mwidx);
						int adidx = m_HexVolumes.at(adjacenHexID[V1sfPos]).GetGlobalIDByLocalID(muidx,mvidx,mwidx);
						vol.m_controlPtGlobalIDs.at(idx++) = adidx;
						continue;
					}
				}
				if(k == 0)
				{
					if(adjacenHexID[W0sfPos] != -1 && m_HexVolumes.at(adjacenHexID[W0sfPos]).m_bGenerateGlobalID)
					{
						IGAHexCell& adjcentBlock = m_HexVolumes.at(adjacenHexID[W0sfPos]);
						int muidx,mvidx,mwidx;
						GetMirrorIJKOnShareSurface(volumeid,adjacenHexID[W0sfPos],2,i,j,k,muidx,mvidx,mwidx);
						int adidx = m_HexVolumes.at(adjacenHexID[W0sfPos]).GetGlobalIDByLocalID(muidx,mvidx,mwidx);
						vol.m_controlPtGlobalIDs.at(idx++) = adidx;
						continue;
					}
				}
				else if(k == wnum-1)
				{
					if(adjacenHexID[W1sfPos] != -1 && m_HexVolumes.at(adjacenHexID[W1sfPos]).m_bGenerateGlobalID)
					{
						IGAHexCell& adjcentBlock = m_HexVolumes.at(adjacenHexID[W1sfPos]);
						int muidx,mvidx,mwidx;
						GetMirrorIJKOnShareSurface(volumeid,adjacenHexID[W1sfPos],2,i,j,k,muidx,mvidx,mwidx);
						int adidx = m_HexVolumes.at(adjacenHexID[W1sfPos]).GetGlobalIDByLocalID(muidx,mvidx,mwidx);
						vol.m_controlPtGlobalIDs.at(idx++) = adidx;
						continue;
					}
				}
				//��������Ѿ�������ϣ���ʣ��������ġ�
				vol.m_controlPtGlobalIDs.at(idx++) = pcount + startID;
				pcount++;
			}
		}
	}
	//Ȼ����Ѱ�ܱ��Ѿ���ֵ��hex��ֱ�ӽ����Ѹ�ֵ�ı��ֵȡ������
	startID += pcount;
	vol.m_bGenerateGlobalID = true;

	//���һ�������Ƿ���ȷ��
	bool checkflag = true;
	for(int i=0; i<maxsize; i++)
	{
		if(vol.m_controlPtGlobalIDs.at(i) == -1)
		{
			checkflag = false;
			break;
		}
	}
	assert(checkflag);   //��ǰ��Ų��ܴ���-1.
	checkflag = false;
	for(int i=0; i<maxsize-1; i++)
	{
		for(int j=i+1; j<maxsize; j++)
		{
			if(vol.m_controlPtGlobalIDs.at(i) == vol.m_controlPtGlobalIDs.at(j))
			{
				checkflag = true; 
				break;
			}
		}
	}
	assert(!checkflag);  //��ǰ��Ų������ظ��ĺ��롣
	int maxid = -1;
	maxsize = vol.m_controlPtGlobalIDs.size();
	for(int i=0; i<maxsize; i++)
	{
		if(vol.m_controlPtGlobalIDs.at(i) > maxid)
			maxid = vol.m_controlPtGlobalIDs.at(i);
	}
	assert(maxid+1 == startID);   //��ǰ�������Ϊ�¸������С��š�
	//��󷵻�һ��ֵ����ֵ������һ��volume�ı����ʼID��
	return startID;
};

PolyIGASolution::PolyIGASolution()
{
	m_pMultiBlockDomain = NULL;
	m_allVolHasSolution = false;
}
PolyIGASolution::~PolyIGASolution()
{

}
void PolyIGASolution::SetMultiVolume(CPolyParaVolume* pBlocks)
{
	if(pBlocks == NULL)
		return;
	m_pMultiBlockDomain = pBlocks;
}
void PolyIGASolution::ClearSolutionData()
{

}
void PolyIGASolution::ReadPatchSolutionData(const char* filename)
{
	//�ڶ�֮ǰ�������ǰһ�������Ľ����

	ClearSolutionData();
	//��ʼ�����ݡ�
	FILE *fp = NULL;
	fp = fopen(filename, "r");
	if(fp == NULL)
	{
		return;
	}
	char	str[200];
	char    tempstr[20];
	Vec4	disv,loadv;
	int     h,pcount;
	int     ctrlPtIdx;
	fgets(str, sizeof(str), fp);
	sscanf( str, "%d", &pcount);
	if(pcount != m_pMultiBlockDomain->m_HexVolumes.size())
	{
	   // AfxMessageBox(_T("�����Ƭ��������Ƭ������ȣ���������"));
	}
	
	h = 0;
	while (fgets(str, sizeof(str), fp) != NULL)
	{
		sscanf( str, "%s", tempstr);
		if(tempstr[0] == 'P')
		{
			h++;
			ctrlPtIdx = 0;
			m_pMultiBlockDomain->m_HexVolumes.at(h-1).m_splVol.m_bHasIGASolution = true;
			continue;
		}
		else
		{
            sscanf( str, "%f%f%f%f%f%f",&disv.x, &disv.y,&disv.z,&loadv.x,&loadv.y,&loadv.z);
			m_pMultiBlockDomain->m_HexVolumes.at(h-1).m_splVol.SetDisplaceAndLoadvector(ctrlPtIdx,disv,loadv);
			ctrlPtIdx++;
		}
	}
	fclose(fp);
	GenerateSolutionRenderData();
}
void PolyIGASolution::GenerateControlpointsDisplacementOrStressColor()
{
	float fminDisplace = 1.0e6;
	float fMaxDisplace = -1.0e6;
	float fminStress = 1.0e6;
	float fMaxStress = -1.0e6;
    int patchNum = m_pMultiBlockDomain->m_HexVolumes.size();
	for(int i=0; i<patchNum; i++)
	{
		SplineVolume& sv = m_pMultiBlockDomain->m_HexVolumes.at(i).m_splVol;
		if(!sv.m_bHasIGASolution)
		{
			sv.m_bHeterogeneous = false;
			continue;
		}
		
		//�����ҵ������С��ֵ��Ȼ���ٹ�һ����
		
		int ctrlNum = sv.m_vAllCtrlPts.size();
		for(int j=0; j<ctrlNum; j++)
		{
			VolumeVertex& svvt = sv.m_vAllCtrlPts.at(j);
			float dis = svvt.m_displacement.Magnitude();
			if(dis > fMaxDisplace)
				fMaxDisplace = dis;
			if(dis < fminDisplace)
				fminDisplace = dis;
			float fstress = svvt.getVonmissStress();
			if(fstress > fMaxStress)
				fMaxStress = fstress;
			if(fstress < fminStress)
				fminStress = fstress;
		}
	}
	for(int i=0; i<patchNum; i++)
	{
		SplineVolume& sv = m_pMultiBlockDomain->m_HexVolumes.at(i).m_splVol;
		int ctrlNum = sv.m_vAllCtrlPts.size();
		for(int j=0; j<ctrlNum; j++)
		{
			VolumeVertex& svvt = sv.m_vAllCtrlPts.at(j);
			float dis = svvt.m_displacement.Magnitude();
		    svvt.m_fDisplaceMat = (dis - fminDisplace)/(fMaxDisplace - fminDisplace);
			float fstress = svvt.getVonmissStress();
			svvt.m_fvonMissStressMat = (fstress - fminStress)/(fMaxStress - fminStress);
			sv.m_vAllCtrlPts.at(j).m_oript = sv.m_vAllCtrlPts.at(j).m_pt;
		}
	}
}
void PolyIGASolution::ChangeRenderDisplaceOrStress(bool disOrStress) //�л�Ӧ����λ�Ƶ���ʾ��
{
	//Ȼ�����ɵȲ���ʾ�㴦�����ݡ�
	if(m_pMultiBlockDomain == NULL)
		return;
	int patchNum = m_pMultiBlockDomain->m_HexVolumes.size();
	for(int i=0; i<patchNum; i++)
	{
		SplineVolume& sv = m_pMultiBlockDomain->m_HexVolumes.at(i).m_splVol;
		if(!sv.m_bHasIGASolution)
		{
			sv.m_bHeterogeneous = false;
			continue;
		}
		sv.m_bHeterogeneous = true;
		int ctrlNum = sv.m_vAllCtrlPts.size();
		for(int j=0; j<ctrlNum; j++)
		{
			if(disOrStress)
			{
				sv.m_vAllCtrlPts.at(j).m_matval = Vec4(sv.m_vAllCtrlPts.at(j).m_fDisplaceMat,0,0);
				sv.m_vAllCtrlPts.at(j).m_pt = sv.m_vAllCtrlPts.at(j).m_oript;
				sv.m_vAllCtrlPts.at(j).m_pt += sv.m_vAllCtrlPts.at(j).m_displacement;
			}
			else
			{
				sv.m_vAllCtrlPts.at(j).m_matval = Vec4(sv.m_vAllCtrlPts.at(j).m_fvonMissStressMat,0,0);
				sv.m_vAllCtrlPts.at(j).m_pt = sv.m_vAllCtrlPts.at(j).m_oript;
			}
		}
		sv.ReCreateAllVolumeRenderPts();
	}
}
void PolyIGASolution::GenerateSolutionRenderData()
{
	//���Ƚ����Ƶ㴦��Ӧ����Ӧ���һ����
	GenerateControlpointsDisplacementOrStressColor();  
	m_allVolHasSolution = true;

	//Ȼ�����ɵȲ���ʾ�㴦�����ݡ�
	if(m_pMultiBlockDomain == NULL)
		return;
	int patchNum = m_pMultiBlockDomain->m_HexVolumes.size();
	for(int i=0; i<patchNum; i++)
	{
		SplineVolume& sv = m_pMultiBlockDomain->m_HexVolumes.at(i).m_splVol;
		if(!sv.m_bHasIGASolution)
		{
			sv.m_bHeterogeneous = false;
			continue;
		}
		sv.m_bHeterogeneous = true;
		int ctrlNum = sv.m_vAllCtrlPts.size();
		for(int j=0; j<ctrlNum; j++)  
		{
			if(true)   //Ĭ��������ʾλ�ơ�
			{
				sv.m_vAllCtrlPts.at(j).m_matval = Vec4(sv.m_vAllCtrlPts.at(j).m_fDisplaceMat,0,0);
				sv.m_vAllCtrlPts.at(j).m_pt = sv.m_vAllCtrlPts.at(j).m_oript;
				sv.m_vAllCtrlPts.at(j).m_pt += sv.m_vAllCtrlPts.at(j).m_displacement;
			}
			else
			{
				sv.m_vAllCtrlPts.at(j).m_matval = Vec4(sv.m_vAllCtrlPts.at(j).m_fvonMissStressMat,0,0);
				sv.m_vAllCtrlPts.at(j).m_pt = sv.m_vAllCtrlPts.at(j).m_oript;
			}
		}
		sv.ReCreateAllVolumeRenderPts();
	}
}
//void PolyIGASolution::RenderDisplacementOrStress(int surfaceMode)
//{
//	if(m_pMultiBlockDomain == NULL)
//		return;
//	int patchNum = m_pMultiBlockDomain->m_HexVolumes.size();
//	for(int i=0; i<patchNum; i++)
//	{
//		SplineVolume& sv = m_pMultiBlockDomain->m_HexVolumes.at(i).m_splVol;
//		if(!sv.m_bHasIGASolution)
//		{
//			sv.m_bHeterogeneous = false;
//			continue;
//		}
//		sv.m_bHeterogeneous = true;
//		sv.RenderbyMode(iRenderBoundary_WithISOElements,surfaceMode);
//	}
//}
};
	