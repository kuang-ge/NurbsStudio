#pragma once

#include "CNurbs.h"
#include "globalFunc.h"
#include <iostream>
#include <sstream>

#include "RWGeometric.h"

#include <vector>
#include "sisl.h"
#include "NurbsTrans.h"

#include <list>
#include <queue>
#include <stack>
#include <map>
#include <algorithm>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/intersections.h>
//SISL
#include<ctime>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Point_2                                         Point_2;
typedef Kernel::Segment_2										Segment_22;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Curve_2                                       Segment_2;
typedef Kernel::Intersect_2										Intersect_2;

using namespace std;

const bool largeAngle = true;//�Ƿ���д���ж�
const double lag = (double)(160 / 180.0)*PI;//��ǶȵĶ���>=160��
const bool byW = false;//��true���������߷�ʽ(����������)����ͳһ��һ����wֵ��С��������
const int sharpMode = 1;//��1��ȥ����㣬�����Ǽ�����
const bool preProcessQuadPol = true;//Ĭ��true��false��ʾ��ǰ���ı�����ȡ���ıߵļ���
const bool priorityQuadPol = false;//Ĭ��false��true���ı��ηŵ��������ʷ�
const bool cgalJudge = true;//true���õ�һ���ж���������õڶ����ж�����

//���߱��
struct PartNurbsLine :public Spline
{
public:
	//��Spline����PartNurbsLine
	void CreatPartNurbsLine(Spline inputLine, int ifSeg);

	//��PartNurbsLine���´���NurnsLine
	void CreatNurbsLine(Spline &newLine);

	//��u���ָ�����
	void PartSegmentation(const double u, varray<PartNurbsLine>& lines);

	//�����׶��㷽������
	Vec3 CalBeginDirecVec();
	Vec3 CalBeginDirecVec()const;

	//����β���㷽������
	Vec3 CalEndDirecVec();
	Vec3 CalEndDirecVec()const;

public:
	int ifSegLine;                          //1Ϊ�ɷָ��ߣ�0Ϊ���ɷָ���
};


//����
class SubPoint
{
public:
	SubPoint() {}
	SubPoint(int isoutp, Vec3 p) :ifOutPoint(isoutp), vetx(p) {}

public:
	int ifOutPoint;                     //1Ϊ���������㣬0Ϊ����������,-1��ʾ�����ضϴ�����
	Vec3 vetx;                          //��������ֵ	dads
};


//����ʹ��(�����ɽ���ش���򻯡���ʹ��stl���Դ���sort������������)
class DandI
{
public:
	double funcVar;
	int order;
};

//�߱�ڵ�
struct EdgeNode
{
	int adjvex;						//�ڽӵ��򣬴洢�ö�����±�
	PartNurbsLine EdgeLine;			//��Ӧ��Nurbs������
	EdgeNode* next;					//����ָ����һ���ڽӵ�
};

//�����ڵ�
struct VertexNode
{
	SubPoint vexdata;				//�����򣬴洢������Ϣ
	EdgeNode* firstedge;			//�߱�ͷָ��
};

//����ͼ
struct GraphAdjList
{
	varray<VertexNode> adjlist;		//���㼯��
	int numNodes, numEdges;			//ͼ�е�ǰ����ͱ���
};

//���ӵ�
struct VisiblePoints
{
	int firstPoint;					//��һ��
	int secondPoint;				//�ڶ���
	double w;						//�������߷�ʽ��Ȩֵ
	int ifConvex;					//�Ƿ���ڰ���,1��ʾ���ڰ���
	int ifSharp;					//�Ƿ���ڼ��,0��ʾ�����ڣ� 1��ʾ����
};

//������е�һ��(pt)����Ӧ�Ŀɷָ���
struct VisibleLine
{
	int pt;							//����
	int beginP;						//���ض���ͷ��
	int endP;						//���ض���β��
	double u;						//�ָ�Ľڵ�ʸ��
	double w;						//Ȩֵ
};

//���������
class MyPolygon
{
public:
	//���һ���ǰ�����
	void GetFBpNumber(int pointNumber, int &frontP, int &behindP);

	//���ĳ���ǰ�����
	int GetFrontNumber(int pointNumber);

public:
	varray<SubPoint> p_Points;				//����
	varray<PartNurbsLine> p_NurbsLines;		//������
	varray<int> objectInAllLines;			//��Ӧ������all_NurbsLines�����е�λ��
	varray<int> convexHull;					//0����͹�㣬1�����㣬180�㹲�ߵ�Ҳ�㰼��
	//varray<VisiblePoints> visiblePointPair;
};


//�ʷ���
class QuadPart 
{
public:
	QuadPart() {}

	QuadPart(varray<Spline> outLines, varray<Spline> inLines, varray<Spline> addLines, varray<int> outQuality, const varray<varray<Vec3>> removePol);

	QuadPart(varray<Spline> outLines, varray<Spline> inLines, varray<Spline> addLines, varray<int> outQuality);

public:
	varray<SubPoint> in_Points;					//���������߶���
	varray<SubPoint> out_Points;				//���������߶���
	varray<SubPoint> all_Points;				//���ж���
	varray<PartNurbsLine> in_NurbsLines;		//����������
	varray<PartNurbsLine> out_NurbsLines;		//����������
	varray<PartNurbsLine> all_NurbsLines;		//��������

	varray<int> out_Number;						//���������߶�Ӧ���
	varray<varray<Vec3>> removepoints_Pol;		//���������εĵ�����
	varray<bool> remove_genus;					//���������εĿ�����
	varray<varray<int>> remove_seg;				//���������ε����߷ָ���

	GraphAdjList g_Adjlist;						//����ͼ

	varray<MyPolygon> quad_Polygon;				//�ı��μ���
	varray<MyPolygon> large_Polygon;			//�����ıߵĶ���μ���
	varray<MyPolygon> remove_Polygon;			//ȥ�������м���Ķ���μ���


public:
	//��ʼ��
	void Init();

	//ִ�������ʷ�����
	void ExecuteQuadPart();

	//�Ƴ������ʷֵĶ����
	void RemovePolygon();

	//�ʷֺ���
	void ResetData();

	//��������������(��)������������(��)�������������㼯���������㼯���Լ����������߼������������߼�
	//outQuality:��outLine��Ŀ��ȣ���Ϊ1���������߿ɷָ��Ϊ0�����ɷָ�
	void CreateAllPointAndLine(varray<Spline> outLines, varray<Spline> inLines, varray<Spline> addLines, varray<int> &outQuality);

	//���������غ��ж�������
	/*
	IL:�ڲ����߼���
	OL:�����߼���
	*/
	void OperateOutandInLines(varray<Spline> &IL, varray<Spline> &OL, varray<int>& outQuality);

	//�����������㰴��ʱ���������
	void OrderOutPointsAntioclock();

	//����ͼ���ڽӱ�ṹ
	void CreateALGraph();

	//�������е㼯���߼������бպ϶����������������ʱ�뷽��)
	void CreatePolygon();

	//��������ļ���
	Vec3 CalGravity(varray<SubPoint> polygenPoints);

	//����β�����(����Զ�����������ĵĳ̶�����)��ͬʱ���ı������ֳ���
	void OrderPolygens();
	
	//�����趨��������ݽṹ(����������)
	void ResetPolygens();

	//�ı��κͶ���Σ��ɷָ��벻�ɷָ�����趨(��ʱ����Ҫ��һ��)

	//����ο��ӵ��Ѱ�ң�����ηŵ�����Ĺ������ٴ���
	void FindVisiblePoint(MyPolygon &pendingPolygen,varray<VisiblePoints> &visiblePointPair);

	//Ѱ�ҵڶ�����ӵ��(������γ������Σ����ڷָ����ߵ�Ԥ����)
	void FindAnotherVp(const MyPolygon &inputPolygen, varray<VisiblePoints> &anotherVp);

	//����εİ����ж�
	void JudgeConcavePoint(MyPolygon &pendingPolygen);

	//���ӵ�Եİ����ж�
	void VisiblePointConcaveJudge(MyPolygon &inputPolygen, varray<VisiblePoints> &visiblePointPair);

	//���ӵ�Եļ���ж�
	void VisiblePointSharpJudeg(MyPolygon &inputPolygen, varray<VisiblePoints>&visiblePointPair, const varray<PartNurbsLine> &allLines);

	//�����߶�֮��uֵ����������µļ���ж�
	bool LinkLineSharpJudge(const Vec3 &frontVec, const Vec3 &lastVec, const Vec3 &pt, const Vec3 &pu);

	//���һ���ǰ�����
	void GetFBpNumber(const MyPolygon &inputPolygen, int pointNumber, int &frontP, int &behindP);

	//���һ���ǰ������
	void GetFBcoord(const MyPolygon &inputPolygen, int pointNumber, Vec3 &frontCoord, Vec3 &behindCoord);

	//���һ���ǰ�������
	void GetVecs(const MyPolygon &inputPolygen, int pointNumber, Vec3 &frontVec, Vec3 &behindVec);

	//��һ���Ȩ����ֵ(���ڰ�������)
	double GetFuncSita(const MyPolygon &inputPolygen, int pointNumber, Vec3 vecAB);

	//��һ��Ȩ����ֵ(����)������
	/*mid:�õ�
	begin:ǰһ��
	end:��һ��
	pt:����һ��*/
	double GetFuncSita(const Vec3 &begin, const Vec3 &end, const Vec3 &mid, const Vec3 &pt);

	//����͹���Ȩ����(ֻ����͹��,��ʱ����д��д��OrderVisiblePoints������)
	double GetConvexFunc();

	//�Կ��ӵ�Խ�������
	void OrderVisiblePoints(const MyPolygon &inputPolygen, varray<VisiblePoints> &visiblePointPair);

	//(ûд��)���ӵ�Եķ���(��������)�����ڼ��&�����ڼ�㣬����ڲ����࣬���հ�͹���Լ�Ȩ������������
	/*���visiblePointPairֻ���Ǽ��ļ���
	sharpVisiblePointPairֻ����㼯��
	*/
	void SortVisiblePoints(MyPolygon &inputPolygen, varray<VisiblePoints> &visiblePointPair, varray<VisiblePoints> &sharpVisiblePointPair);

	//��ö����ĳ�����һ����
	int GetNextPoint(const MyPolygon &inputPolygen, const int &pNumber);

	//���ָ���uֵ
	/*begin:ǰһ��;
	end:��һ��
	crossp:���������е�ĳһ��*/
	double GetPartPointU(const MyPolygon &inputPolygen, const varray<PartNurbsLine> &allLines, const Vec3 &begin, const Vec3 &end, const Vec3 &crossp);

	//Ѱ�Ҷ�����еĿɷָ��ߣ����ҵ����ڵĿ��ӵ�(δ���������������)���������������u(�����ȴ�ֱ���ٽ�ƽ�֣����е�)��
	//�������������Ȩֵ������
	void CalPartLines(const MyPolygon &inputPolygen, const varray<PartNurbsLine> &allLines, varray<VisibleLine> &partLine, const int &mode );

	//ͨ�����㴴��2��Nurbs����(������Ҫ�������������)
	void CreatPartNurbsLine(const Vec3 &p1, const Vec3 &p2, PartNurbsLine &line);

	//�������߲���(Ĭ�ϵ�һ����ȵڶ�����λ�ÿ�ǰ)
	bool ConncetTwoPoints(MyPolygon &pendingPol, varray<MyPolygon> &quadPol, varray<MyPolygon> &morePol, varray<PartNurbsLine> &allLines, VisiblePoints &visibleP);
	
	//���ӵ��ֱ��u����
	bool ConnectLineWithPoint(MyPolygon &pendingPol, varray<MyPolygon> &quadPol, varray<MyPolygon> &morePol, varray<PartNurbsLine> &allLines, VisibleLine &visibleL);



	////////////////////////////////////////////////////////////////////////
	//�ʷֲ���
	/*
	mode:����ģʽ��0Ϊ����ģʽ��1Ϊ�ϸ�ģʽ(��Լ��������)*/
	bool QuadOperation(const varray<MyPolygon> &quadPolygen, const varray<MyPolygon> &morePolygen, const varray<PartNurbsLine> &allLines, varray<MyPolygon> &allPolygen, varray<PartNurbsLine> &endLines, int mode);
};

//����������ڵ�
class SfCtainTreeNode
{
public:
	SfCtainTreeNode();

	//outlines:�������Χ������
	SfCtainTreeNode(const varray<Spline> &outlines);

	//conlines:���������������
	SfCtainTreeNode(const varray<Spline> &outlines, const varray<Spline> &conlines);

	SfCtainTreeNode(const varray<Spline> &outlines, const varray<Spline> &conlines, varray<int> seg);

	//genus:�Ƿ�Ϊ��������,false��ʾ������
	SfCtainTreeNode(const varray<Spline> &outlines, const varray<Spline> &conlines, varray<int> seg, bool genus);

	//��������
	~SfCtainTreeNode();


	//��ȡ��ǰ�ڵ�����������
	void GetOutPoints(varray<Vec3>& ps);

	//��ȡ��ǰ�ڵ�����������
	void GetOutLines(varray<Spline>& outline);

	//��ȡ��ǰ�ڵ�����
	void GetSurfs(varray<SplineSurface> &surfs);


public:
	varray<Spline> outLines;			//��Χ������
	varray<Spline> conLines;			//����������
	varray<Vec3> vertex;				//��������������
	varray<int> isSeg;					//���������ߵĿɷָ���
	bool isGenus;						//�Ƿ����

	varray<Spline> allLines;			//��������
	varray<int> outNumber;				//��Χ���������(�Ѿ�����ʷֺ��)
	varray<varray<int>> quadPolNumber;	//�����ı��ζ�Ӧ�������
	//varray<int> conNumber;			//�������������
	//SfCtainTreeNode* child;			//�ӽڵ�
	//SfCtainTreeNode* brother;			//�ֵܽڵ�
	list<SfCtainTreeNode*> childs;		//�������������ӽڵ�
	bool isSort;						//�����Ƿ񾭹�����
	int num;							//������ţ�����ʱ��ʱʹ��
	bool quaded;						//�Ƿ�����ʷ�
};


//����ȥ��
void DecoincideNurbsLine(varray<Spline>& nl);

//����ÿ���������ڵ�(��Ϊ���ڵ�ʱ)
void ReSetTreeNode(SfCtainTreeNode* father);

//�������������
SfCtainTreeNode* CreateSurfContainTree(const varray<varray<Spline>> &surf, const varray<Spline> conlines, varray<varray<int>> seg, varray<bool> genus);
//�������������(���������ʹ�øú���)
//surf:�������ߣ�[0]λ��Ϊ����������
//conlines:��������������
//seg:���߿ɽض���
//genus:������(�Ƿ�Ϊ��)trueΪ�� falseΪʵ��
//sfnum:�������
SfCtainTreeNode* CreateSurfContainTree(const varray<varray<Spline>> &surf, const varray<Spline> conlines, varray<varray<int>> seg, varray<bool> genus, varray<int> sfnum);

/////////////////////////////////////////
//
//
//�ʷֵ�ǰ�ڵ��������Σ�ÿһ�����ʷֺ���
void QuadTreeNodeSurf(SfCtainTreeNode* root);

//������������������ʷֲ���
void QuadWithContainTree(SfCtainTreeNode* root);

//��root�ڵ㿪ʼ��������������Ϣ
void GetQuadSurfaces(SfCtainTreeNode* root, varray<SplineSurface> &surf);

//ִ�����а���������
void StartTreeOperation(const varray<varray<Spline>>& surf, const varray<Spline>& conlines, const varray<varray<int>> seg, varray<bool> genus, varray<Spline> &allLines, varray<SplineSurface> &allSurf);

//ͨ�����㴴��2��Nurbs����(������Ҫ�������������)
Spline CreatNurbsLine(const Vec3 &p1, const Vec3 &p2);


//���һ�����ĳ��ǰ������
void GetLinesByPoint(const MyPolygon &inputPolygen, int pointNumber, const varray<PartNurbsLine> &allLine, Spline &frLine, Spline &backLine);

//��ȡһ�����ĳ��ǰ�����߶�Ӧ������
void GetVecsByPoint(const MyPolygon &inputPolygen, int pointNumber, const varray<PartNurbsLine> &allLine, Vec3 &frVec, Vec3 &backVec);

//���һ�����ĳ���ǰ������
void GetFBcoord(const MyPolygon &inputPolygen, int pointNumber, Vec3 &frontCoord, Vec3 &behindCoord);


//ͨ�����㴴��2��Nurbs����,��������˵����(���������������)
//pol:������ɵĶ����
PartNurbsLine GetPartNurbsLine(const Vec3 &p1, const Vec3 &p2, const varray<PartNurbsLine> &allLine, MyPolygon pol);

//ͨ�����㴴��2��Nurbs����,���һ���˵�,һ�����ض��������
//p1:���߶˵�����
//frontp:���ض���������Ӧ�ڶ������ǰһ������
//u:���ضϵ�Ĳ���
PartNurbsLine GetPartNurbsLine(const Vec3 &p1, const Vec3 &frontp, double u, const varray<PartNurbsLine> &allLine, MyPolygon pol);


//���ָ���uֵ
//begin:ǰһ��;
//end:��һ��
//crossp:���������е�ĳһ��
double GetPartU (const Spline &lines, const Vec3 &begin, const Vec3 &end, const Vec3 &crossp);

//����ǰ������(��ʱ��)�Լ�Բ�����꣬�����NURBSԲ��(��xyƽ����)
void CreatArc(const Vec3 &p1, const Vec3 &p2, const Vec3 &centre, Spline &line);

//����һ�ı����ı߽�������(������Coons��ֵ������)
void QuadOrder(const MyPolygon &quad,varray<PartNurbsLine> &allLine ,varray<Spline> &coonsLine);

//����һ�ı����ı߽�������(������Coons��ֵ������)
void QuadOrder(const varray<int> &quadNums, const varray<Spline> &allLine, varray<Spline> &coonsLine);

//�ж������Ƿ��غ�
bool JudgeTwoPointsCoincide(const Vec3 &p1, const Vec3 &p2);

//�ж���Nurbs�����Ƿ��غ�(ֻ��ͷβ����)
bool JudgeTwoLinesCoincide(const Spline &L1, const Spline &L2);

//����������ĽǶ�
double GetAngle(Vec3 vec1, Vec3 vec2);

//����ֱ�ߵĽ���(��ά������ͬһƽ�棬zֵ���),δ������ֱ���غ����
Vec3 GetTwoLinesCrossPoint(const Vec3 &p1, const Vec3 &p2, const Vec3 &p3, const Vec3 &p4);

//��ά�ռ�㵽ֱ�ߵĴ���(�������߶�����)
Vec3 GetFootPerpendicular(const Vec3 &begin, const Vec3 &end, const Vec3 &pt);

//�ж��Ƿ���ڸý���,�Ҳ�Ϊ������
bool JudegeFootPerpendicular(const Vec3 &begin, const Vec3 &end, const Vec3 &crossPoint);

//����������ƽ��������
Vec3 GetAngularBisectorVec(const Vec3 &leftVec, const Vec3 &rigthtVec);

//����ƽ������������ֱ�ߺ�ԭֱ�ߵĽ���(֮���ٽ������������ж��Ƿ���ڸý���)
Vec3 GetAngularBisectorPoint(const Vec3 &begin, const Vec3 &end,const Vec3 &pt, const Vec3 &vec);

//����(�����߱��)�ռ�任(�����������޸�)
void SurfaceConverse(varray<Spline> &lines, Vec3 targetPts, Vec3 targetVec, Matrix4d &conMx, Matrix4d &reconMx);

//���߿ռ�任
void LineConverse(varray<Spline> &lines, Matrix4d &conMx);

//���һ�����u�Ƿ�Ϊ���߶˵�
int CalU(double u);

//����ɱպ϶���ε��������߰���ʱ������,������ʱ�����������
//polLines:���������α����߼���
varray<Vec3> OrderLinesAntioclock(varray<Spline> &polLines);


//�����ȷu
bool CalRightu(double u, double y, Spline Nline, vector<double> &us);

bool CalRightu(vector<double> u, double y, Spline Nline, vector<double> &us);

bool CalRightu(vector<float> u, float y, Spline Nline, vector<float> &us);

//�߶���NURBS(SISL)�����󽻵�
//���ؽ�����
int LineIntersectNurbs(SISLCurve *& curve, const varray<Vec3>& strLine, const Spline * nurbsCurve = nullptr);

//���һ��������NURBS����Χ�ɵĶ����(��/��/��)
//����ֵΪ-1:error
//����ֵΪ0:�ڶ������(��Ϊ�˵㴦)
//����ֵΪ1:�ڶ�����ڲ�
//����ֵΪ2:�ڶ������
//����ֵΪ3:�ڶ������(�Ҳ�Ϊ�˵㴦)
int PointRelatePolygon(Vec3 p, varray<Spline> &pol);

//��NURBS������Χ����ΰ�����ϵ,pola�Ƿ���polb�ڲ�
bool PolygonsRelationship(varray<Spline> pola, varray<Spline> polb);



