#include "stdafx.h"
#include "Fragment.h"

Fragment::Fragment()
{

}

Fragment::Fragment(float a, float b, float c, float begin, float end, float degree)
	:_A(a)
	,_B(b)
	,_C(c)
	,_begin(begin)
	,_end(end)
	,_degree(degree)
{
	//����_degreeҪ��һ����Ϊ�г�����ƽ��ÿ�һЩ��������������棬�͸ոպ���ͷ��һ�£�β����һ��
	_each = (_end - _begin) / (_degree-1) ;
}



Fragment::~Fragment()
{
}

//vector<float> Fragment::ReturnCurve()
//{
//	vector<float> temp;
//	for (int i = 0; i < Patchnum; i++)
//	{
//		for (int j = 0; j < PatchHex[i].size(); j++)
//		{
//
//			for (int k = 0; k < 8; k++)
//			{
//				temp.push_back(PatchHex[i][j].m_ptidx[k].x);
//				temp.push_back(PatchHex[i][j].m_ptidx[k].y);
//				temp.push_back(PatchHex[i][j].m_ptidx[k].z);
//			}
//		}
//	}
//	return temp;
//}

vector<float> Fragment::ReturnSupportPoints()
{
	vector<float> temp;
	for (int i = 0; i < SupportPointsCube.size(); i++)
	{
		temp.push_back(SupportPointsCube[i].x);
		temp.push_back(SupportPointsCube[i].y);
		temp.push_back(SupportPointsCube[i].z);
	}
	return temp;
}

//vector<float> Fragment::ReturnFourBoundary()
//{
//	vector<float> temp;
//	for (int i = 0; i < AllBoundaryEdges.size(); i++)
//	{
//		for (int j = 0; j < AllBoundaryEdges[i].size(); j++)
//		{
//			temp.push_back(AllBoundaryEdges[i][j].NodePoint[0].x);
//			temp.push_back(AllBoundaryEdges[i][j].NodePoint[0].y);
//			temp.push_back(AllBoundaryEdges[i][j].NodePoint[0].z);
//			temp.push_back(AllBoundaryEdges[i][j].NodePoint[1].x);
//			temp.push_back(AllBoundaryEdges[i][j].NodePoint[1].y);
//			temp.push_back(AllBoundaryEdges[i][j].NodePoint[1].z);
//		}
//	}
//	return temp;
//}

//vector<float> Fragment::ReturnSurface()
//{
//	vector<float> temp;
//	for (int i = 0; i < m_allsurface.size(); i++)
//	{
//		for (int j = 0; j < m_allsurface[i].Intersectpoints.size(); j++)
//		{
//			if (m_allsurface[i].Intersectpoints[j].size() > 0)
//			{
//				for (int k = 0; k < m_allsurface[i].Intersectpoints[j].size(); k++)
//				{
//					temp.push_back(m_allsurface[i].Intersectpoints[j][k].x);
//					temp.push_back(m_allsurface[i].Intersectpoints[j][k].y);
//					temp.push_back(m_allsurface[i].Intersectpoints[j][k].z);
//				}
//			}
//		}
//	}
//	return temp;
//}

point3d Fragment::alf_mult_point(double alf, point3d p)
{
	p.x = alf * p.x;
	p.y = alf * p.y;
	p.z = alf * p.z;
	return p;
}

//int Fragment::Findnum(int left, int right, double x, varray<double> kont)
//{
//	int mid;
//	if (x >= kont[right])
//		return right - 1;
//	mid = (left + right) / 2;
//	while (x < kont[mid] || x >= kont[mid + 1])
//	{
//		if (x < kont[mid])
//			right = mid;
//		else
//			left = mid;
//		mid = (right + left) / 2;
//	}
//	return mid;
//}
//
//void Fragment::QuadsAndHexs()
//{
//	for (int i = 0; i < Patchnum; i++)//����
//	{
//		GetAllCell(u_num - 1, v_num - 1, w_num - 1, u_order - 1, v_order - 1, w_order - 1, u, v, w, PatchControlPoints[i]); //�洢�Ȳ�������Ԫ
//	}
//}
//
//void Fragment::SolvingSection()
//{
//	for (int i = 0; i < 1; i++) //��������
//	{
//		_D = _D - 8;
//		//����ƽ����Ƭ
//		CIntersecPoints();               //�������ϵĵ�
//		if (m_allIntersetPoints.size() > 0)
//		{
//			SurfaceFitting();            //ƽ�����
//			CoordinateTransformation();  //ͶӰ��XOYƽ��
//			FindBoundaryPoints();
//			ClearHex();
//		}
//	}
//}
//
//void Fragment::SolvingSupport()
//{
//	FindBoundaryCurve();
//	GetSupporPoints();
//	ReturnSupport();
//}
//
//void Fragment::ReturnSupport()
//{
//	point3d pt;
//	for (int i = 0; i < m_allSupportPoints.size(); i++)
//	{
//		for (int j = 0; j < m_allSupportPoints[i].size(); j++)
//		{
//			pt.x = Tran[i][0] * m_allSupportPoints[i][j].x + Tran[i][1] * m_allSupportPoints[i][j].y + Tran[i][2] * m_allSupportPoints[i][j].z + Tran[i][3];
//			pt.y = Tran[i][4] * m_allSupportPoints[i][j].x + Tran[i][5] * m_allSupportPoints[i][j].y + Tran[i][6] * m_allSupportPoints[i][j].z + Tran[i][7];
//			pt.z = Tran[i][8] * m_allSupportPoints[i][j].x + Tran[i][9] * m_allSupportPoints[i][j].y + Tran[i][10] * m_allSupportPoints[i][j].z + Tran[i][11];
//			SupportPointsCube.push_back(pt);
//		}
//	}
//}
//
//void Fragment::GetSupporPoints()
//{
//	varray<point3d> Support;
//
//	for (int i = 0; i < m_allsurface.size(); i++)
//	{
//		for (int l = 0; l < m_allsurface[i].Transformpoints.size(); l++)
//		{
//			if (m_allsurface[i].Transformpoints[l].size() > 0)
//			{
//				if (i == 0)
//				{
//					for (int m = 0; m < m_allsurface[i].Transformpoints[l].size() - 1; m++)
//					{
//						Support.push_back(m_allsurface[i].Transformpoints[l][m]);
//					}
//				}
//				else
//				{
//					for (int m = 0; m < m_allsurface[i].Transformpoints[l].size(); m++)
//					{
//						if (isPolygonContainsPoint(AllBoundaryEdges[i - 1], m_allsurface[i].Transformpoints[l][m]) == false &&
//							isPointInPolygonBoundary(AllBoundaryEdges[i - 1], m_allsurface[i].Transformpoints[l][m]) == false)
//						{
//							Support.push_back(m_allsurface[i].Transformpoints[l][m]);
//						}
//					}
//				}
//			}
//		}
//		m_allSupportPoints.push_back(Support);
//		Support.clear();
//	}
//}
//
//bool Fragment::isPointInPolygonBoundary(varray<Quadedge> curvepoints, point3d pending)
//{
//	for (int i = 0; i < curvepoints.size(); i++)
//	{
//		if (pending.y < min(curvepoints[i].NodePoint[0].y, curvepoints[i].NodePoint[1].y))
//		{
//			continue;
//		}
//		if (pending.y > max(curvepoints[i].NodePoint[0].y, curvepoints[i].NodePoint[1].y))
//		{
//			continue;
//		}
//		if (fabs(curvepoints[i].NodePoint[0].y - curvepoints[i].NodePoint[1].y) < 0.001)
//		{
//			double minX = min(curvepoints[i].NodePoint[0].x, curvepoints[i].NodePoint[1].x);
//			double maxX = max(curvepoints[i].NodePoint[0].x, curvepoints[i].NodePoint[1].x);
//			if ((fabs(pending.y - curvepoints[i].NodePoint[0].y) < 0.001) && (fabs(pending.y - curvepoints[i].NodePoint[1].y) < 0.001) && (pending.x > minX && pending.x < maxX))
//			{
//				return true;
//			}
//		}
//		else
//		{
//			double x = (pending.y - curvepoints[i].NodePoint[0].y)*(curvepoints[i].NodePoint[1].x - curvepoints[i].NodePoint[0].x) / (curvepoints[i].NodePoint[1].y - curvepoints[i].NodePoint[0].y) + curvepoints[i].NodePoint[0].x;
//			if (fabs(x - pending.x) < 0.001)
//			{
//				return true;
//			}
//		}
//	}
//	return false;
//}
//
//bool Fragment::isPolygonContainsPoint(varray<Quadedge> curvepoints, point3d pending)
//{
//	int num = 0;
//
//	for (int i = 0; i < curvepoints.size(); i++)
//	{
//		if (fabs(curvepoints[i].NodePoint[0].y - curvepoints[i].NodePoint[1].y) < 0.001)
//		{
//			continue;
//		}
//		if (pending.y < min(curvepoints[i].NodePoint[0].y, curvepoints[i].NodePoint[1].y))
//		{
//			continue;
//		}
//		if (pending.y > max(curvepoints[i].NodePoint[0].y, curvepoints[i].NodePoint[1].y))
//		{
//			continue;
//		}
//		double x = (pending.y - curvepoints[i].NodePoint[0].y)*
//			(curvepoints[i].NodePoint[1].x - curvepoints[i].NodePoint[0].x) 
//			/ (curvepoints[i].NodePoint[1].y - curvepoints[i].NodePoint[0].y) + curvepoints[i].NodePoint[0].x;
//		if (x > pending.x)
//		{
//			num++;
//		}
//	}
//	return  (num % 2 == 1);
//}
//
//bool Fragment::CompareTwoEggs(Quadedge One, Quadedge Two)
//{
//	if ((fabs(One.NodePoint[0].x - Two.NodePoint[0].x) < 0.001 && fabs(One.NodePoint[0].y - Two.NodePoint[0].y) < 0.001 && fabs(One.NodePoint[1].x - Two.NodePoint[1].x) < 0.001 && fabs(One.NodePoint[1].y - Two.NodePoint[1].y) < 0.001) ||
//		(fabs(One.NodePoint[0].x - Two.NodePoint[1].x) < 0.001 && fabs(One.NodePoint[0].y - Two.NodePoint[1].y) < 0.001 && fabs(One.NodePoint[1].x - Two.NodePoint[0].x) < 0.001 && fabs(One.NodePoint[1].y - Two.NodePoint[0].y) < 0.001))
//	{
//		return true;
//	}
//	else
//	{
//		return false;
//	}
//}
//
//void Fragment::FindBoundaryCurve()
//{
//	varray<Quadedge> BoundaryEdges, AllBoundaryEdge; 
//	varray<varray<Quadedge>> BoundaryEdgess;
//	int flag;
//
//	for (int i = 0; i < m_allsurface.size(); i++)  //��������
//	{
//		for (int j = 0; j < m_allsurface[i].Transformpoints.size(); j++) //��Ԫ����
//		{
//			for (int k = 0; k < m_allsurface[i].Transformpoints[j].size() - 1; k++) //��Ԫ���������
//			{
//				if (m_allsurface[i].Transformpoints[j][k].boundary == true && m_allsurface[i].Transformpoints[j][k + 1].boundary == true)
//				{
//					Quadedge bound;
//					ChangePoints(bound.NodePoint[0], m_allsurface[i].Transformpoints[j][k]);
//					ChangePoints(bound.NodePoint[1], m_allsurface[i].Transformpoints[j][k + 1]);
//					BoundaryEdges.push_back(bound);
//				}
//			}
//		}
//		BoundaryEdgess.push_back(BoundaryEdges);
//		BoundaryEdges.clear();
//	}
//
//	for (int i = 0; i < BoundaryEdgess.size(); i++)
//	{
//		for (int j = 0; j < BoundaryEdgess[i].size(); j++)
//		{
//			flag = 0;
//			for (int k = 0; k < BoundaryEdgess[i].size(); k++)
//			{
//				if (j == k)
//				{
//					continue;
//				}
//				if (CompareTwoEggs(BoundaryEdgess[i][j], BoundaryEdgess[i][k]))
//				{
//					flag++;
//				}
//			}
//			if (flag == 0)
//			{
//				AllBoundaryEdge.push_back(BoundaryEdgess[i][j]);
//			}
//		}
//		AllBoundaryEdges.push_back(AllBoundaryEdge);
//		AllBoundaryEdge.clear();
//	}
//}
//
//void Fragment::ClearHex()
//{
//	for (int i = 0; i < Patchnum; i++)
//	{
//		for (int j = 0; j < PatchHex[i].size(); j++)
//		{
//			PatchHex[i][j].Intersectpoints.clear();
//			PatchHex[i][j].Transformpoints.clear();
//		}
//	}
//	m_allIntersetPoints.clear();
//}
//
//void Fragment::FindBoundaryPoints()
//{
//	Section SmallHex;
//	varray<point3d> Intersectpoint, Transformpoint;
//
//	for (int i = 0; i < Patchnum; i++)//�������ģ�͵�Ƭ��
//	{
//		for (int j = 0; j < PatchHex[i].size(); j++)//���еȲ�������Ԫ
//		{
//			if (PatchHex[i][j].Transformpoints.size() > 0)
//			{
//				for (int k = 0; k < PatchHex[i][j].Transformpoints.size(); k++)
//				{
//					Intersectpoint.push_back(PatchHex[i][j].Intersectpoints[k]);
//					Transformpoint.push_back(PatchHex[i][j].Transformpoints[k]);
//				}
//				SmallHex.Intersectpoints.push_back(Intersectpoint);
//				SmallHex.Transformpoints.push_back(Transformpoint);
//				Intersectpoint.clear();
//				Transformpoint.clear();
//			}
//		}
//	}
//
//	m_allsurface.push_back(SmallHex);
//}
//
//void Fragment::CoordinateTransformation()
//{
//	
//	Eigen::Matrix4d RX, RY, T, TRAN, InverTRAN;
//	point3d pt, ptt;
//	Eigen::Vector4d OldPoint, NewPoint;
//	varray<double> Matrix;
//
//	if (_C != 0)
//	{
//		pt.x = 0, pt.y = 0, pt.z = -_D / _C;
//	}
//	if (_B != 0)
//	{
//		pt.x = 0, pt.y = -_D / _B, pt.z = 0;
//	}
//	if (_A != 0)
//	{
//		pt.x = -_D / _A, pt.y = 0, pt.z = 0;
//	}
//	if (_B == 0 && _C == 0)
//	{
//		cosrx = 1;
//		sinrx = 0;
//	}
//	else
//	{
//		cosrx = _C / (sqrt(pow(_B, 2) + pow(_C, 2)));
//		sinrx = _B / (sqrt(pow(_B, 2) + pow(_C, 2)));
//	}
//
//	cosry = (sqrt(pow(_B, 2) + pow(_C, 2))) / (sqrt(pow(_A, 2) + pow(_B, 2) + pow(_C, 2)));
//	sinry = -_A / (sqrt(pow(_A, 2) + pow(_B, 2) + pow(_C, 2)));
//
//	T << 1, 0, 0, -pt.x,
//		0, 1, 0, -pt.y,
//		0, 0, 1, -pt.z,
//		0, 0, 0, 1;
//
//	RX << 1, 0, 0, 0,
//		0, cosrx, -sinrx, 0,
//		0, sinrx, cosrx, 0,
//		0, 0, 0, 1;
//
//	RY << cosry, 0, sinry, 0,
//		0, 1, 0, 0,
//		-sinry, 0, cosry, 0,
//		0, 0, 0, 1;
//
//	TRAN = RY * RX*T;
//
//	for (int i = 0; i < Patchnum; i++)
//	{
//		for (int j = 0; j < PatchHex[i].size(); j++)
//		{
//			for (int k = 0; k < PatchHex[i][j].Intersectpoints.size(); k++)
//			{
//				OldPoint << PatchHex[i][j].Intersectpoints[k].x, PatchHex[i][j].Intersectpoints[k].y, PatchHex[i][j].Intersectpoints[k].z, 1;
//				NewPoint = TRAN * OldPoint;
//				ChangePoints(ptt, PatchHex[i][j].Intersectpoints[k]);
//				ptt.x = NewPoint[0];
//				ptt.y = NewPoint[1];
//				ptt.z = NewPoint[2];
//				PatchHex[i][j].Transformpoints.push_back(ptt);
//			}
//		}
//	}
//	InverTRAN = TRAN.inverse();
//	for (int i = 0; i < 4; i++)
//	{
//		for (int j = 0; j < 4; j++)
//		{
//			Matrix.push_back(InverTRAN(i, j));
//		}
//	}
//	Tran.push_back(Matrix);
//}
//
//void Fragment::ChangePoints(point3d & p1, point3d & p2)
//{
//	p1.x = p2.x;
//	p1.y = p2.y;
//	p1.z = p2.z;
//	p1.u = p2.u;
//	p1.v = p2.v;
//	p1.w = p2.w;
//	p1.Gray = p2.Gray;
//	p1.boundary = p2.boundary;
//}
//
//double Fragment::Angle(point3d pt1, point3d pt2, point3d pt3)
//{
//	point3d pt, ptsin;
//	double product, productpt, productptt, angle, angleh, division;
//	pt.x = pt2.x - pt3.x;
//	pt.y = pt2.y - pt3.y;
//	pt.z = pt2.z - pt3.z;
//
//	ptsin.x = pt1.y*pt.z - pt1.z*pt.y;
//	ptsin.y = pt1.z*pt.x - pt1.x*pt.z;
//	ptsin.z = pt1.x*pt.y - pt1.y*pt.x;
//
//	product = pt1.x*pt.x + pt1.y*pt.y + pt1.z*pt.z;
//	productpt = sqrt(pow(pt1.x, 2) + pow(pt1.y, 2) + pow(pt1.z, 2));
//	productptt = sqrt(pow(pt.x, 2) + pow(pt.y, 2) + pow(pt.z, 2));
//
//	division = product / (productpt*productptt);
//
//	if (division > 1)
//	{
//		division = 1;
//	}
//	else if (division < -1)
//	{
//		division = -1;
//	}
//	angleh = acos(division);
//	angle = (angleh * 180) / 3.141592654;
//
//	if (ptsin.x + ptsin.y + ptsin.z >= 0)
//	{
//		return angle;
//	}
//	else
//	{
//		return 360 - angle;
//	}
//}

//void Fragment::BubbleSort(Hex & hexagon)
//{
//	point3d Centerpoint, firstcondition, ptt;
//	double ptangle;
//	varray<double> angle;
//
//	for (int i = 0; i < hexagon.Intersectpoints.size(); i++)
//	{
//		Centerpoint.x = Centerpoint.x + hexagon.Intersectpoints[i].x;
//		Centerpoint.y = Centerpoint.y + hexagon.Intersectpoints[i].y;
//		Centerpoint.z = Centerpoint.z + hexagon.Intersectpoints[i].z;
//	}
//
//	Centerpoint.x = Centerpoint.x / hexagon.Intersectpoints.size();
//	Centerpoint.y = Centerpoint.y / hexagon.Intersectpoints.size();
//	Centerpoint.z = Centerpoint.z / hexagon.Intersectpoints.size();
//
//	firstcondition.x = hexagon.Intersectpoints[0].x - Centerpoint.x;
//	firstcondition.y = hexagon.Intersectpoints[0].y - Centerpoint.y;
//	firstcondition.z = hexagon.Intersectpoints[0].z - Centerpoint.z;
//
//	for (int i = 0; i < hexagon.Intersectpoints.size(); i++)
//	{
//		angle.push_back(Angle(firstcondition, hexagon.Intersectpoints[i], Centerpoint));
//	}
//
//	for (int i = 0; i < hexagon.Intersectpoints.size(); i++)
//	{
//		for (int j = 0; j < hexagon.Intersectpoints.size() - 1 - i; j++)
//		{
//			if (angle[j] > angle[j + 1])
//			{
//				ptangle = angle[j];
//				angle[j] = angle[j + 1];
//				angle[j + 1] = ptangle;
//				ChangePoints(ptt, hexagon.Intersectpoints[j]);
//				ChangePoints(hexagon.Intersectpoints[j], hexagon.Intersectpoints[j + 1]);
//				ChangePoints(hexagon.Intersectpoints[j + 1], ptt);
//			}
//		}
//	}
//	ChangePoints(ptt, hexagon.Intersectpoints[0]);
//	hexagon.Intersectpoints.push_back(ptt);
//}

//void Fragment::SurfaceFitting()
//{
//	point3d firstcondition;
//	int num = 0;
//	for (int i = 0; i < Patchnum; i++)
//	{
//		for (int j = 0; j < PatchHex[i].size(); j++)
//		{
//			if (PatchHex[i][j].Intersectpoints.size() > 2)
//			{
//				BubbleSort(PatchHex[i][j]);
//			}
//		}
//	}
//}
//
//void Fragment::intersectionLinePlane(point3d p1, point3d p2, Hex & hexagon)
//{
//	point3d pt, pt1, ptl, ptr;
//	double num, den, n, distancept, distanceptl, distanceptr;
//	pt.x = p2.x - p1.x;
//	pt.y = p2.y - p1.y;
//	pt.z = p2.z - p1.z;
//	pt.u = p2.u - p1.u;
//	pt.v = p2.v - p1.v;
//	pt.w = p2.w - p1.w;
//	pt.Gray = p2.Gray - p1.Gray;
//
//	distancept = sqrt(pow(pt.x, 2) + pow(pt.y, 2) + pow(pt.z, 2));
//
//	num = _A * p1.x + _B * p1.y + _C * p1.z + _D; //��ĸ
//	den = -_A * pt.x - _B * pt.y - _C * pt.z;   //����
//
//	if (fabs(den) > 1e-5)
//	{
//		n = num / den;
//		pt1.x = p1.x + n * pt.x;
//		pt1.y = p1.y + n * pt.y;
//		pt1.z = p1.z + n * pt.z;
//		pt1.u = p1.u + n * pt.u;
//		pt1.v = p1.v + n * pt.v;
//		pt1.w = p1.w + n * pt.w;
//		pt1.Gray = p1.Gray + n * pt.Gray;
//
//		ptl.x = pt1.x - p1.x;
//		ptl.y = pt1.y - p1.y;
//		ptl.z = pt1.z - p1.z;
//
//		ptr.x = pt1.x - p2.x;
//		ptr.y = pt1.y - p2.y;
//		ptr.z = pt1.z - p2.z;
//
//		distanceptl = sqrt(pow(ptl.x, 2) + pow(ptl.y, 2) + pow(ptl.z, 2));
//		distanceptr = sqrt(pow(ptr.x, 2) + pow(ptr.y, 2) + pow(ptr.z, 2));
//
//		if ((distanceptl / distancept < 1 && distanceptr / distancept < 1) || distanceptl == 0 || distanceptr == 0)
//		{
//			if (pt1.u == 0 || pt1.u == 1 || pt1.v == 0 || pt1.v == 1 || pt1.w == 0 || pt1.w == 1)
//			{
//				pt1.boundary = true;
//			}
//			hexagon.Intersectpoints.push_back(pt1);
//			hexagon.Intersect = true;
//			m_allIntersetPoints.push_back(pt1);
//		}
//	}
//}

point3d Fragment::getShadow(point3d p, float D)//����ͶӰ��
{
	point3d shadow;
	float den = pow(_A, 2) + pow(_B, 2) + pow(_C, 2);
	shadow.x = (pow(_B, 2) + pow(_C, 2))*p.x - _A * (_B*p.y + _C * p.z + D);
	shadow.x /= den;
	shadow.y = (pow(_A, 2) + pow(_C, 2))*p.y - _B * (_A*p.x + _C * p.z + D);
	shadow.y /= den;
	shadow.z = (pow(_A, 2) + pow(_B, 2))*p.z - _C * (_A*p.x + _B * p.y + D);
	shadow.z /= den;
	return shadow;
}

bool Fragment::SetInOrOut(point3d p ,int n)//�жϵ��Ƿ�������,����n��ʾ�ڼ�������
{
	point3d shadow = getShadow(p, _begin + n * _each);//ȡ������һ�����������ͶӰ��

	int nCross = 0,nCross1=0;//�ֱ��¼��ߵĽ�����ұߵĽ���
	point3d lastp1 = { 0,0,0 };
	point3d lastp2 = { 0,0,0 };
	for (int i = 0; i != eachSurfacePoints_point3d[n].size(); i+=2)
	{
		point3d p1 = eachSurfacePoints_point3d[n][i];
		point3d p2 = eachSurfacePoints_point3d[n][i+1];//�ֱ������������һ���߶�


		if (lastp1 != p1||lastp2!=p2)
		{
			
			if (p.z < min(p1.z, p2.z))
				continue;
			if (p.z > max(p1.z, p2.z))
				continue;
			//�鿴�Ƿ��ڱ߽����
			if (fabs(p1.z - p2.z) < 0.001)
			{
				float minx = min(p1.x, p2.x);
				float maxx = max(p1.x, p2.x);
				if ((fabs(p.z - p1.z) < 0.001) && (fabs(p.z - p2.z) < 0.001) && (p.x >= minx && p.x <= maxx))
				{
					return true;
				}
			}
			//ͬ�����в鿴�Ƿ��ڱ߽����
			if (fabs(p1.x - p2.x) < 0.001)
			{
				float minx = min(p1.z, p2.z);
				float maxx = max(p1.z, p2.z);
				if ((fabs(p.x - p1.x) < 0.001) && (fabs(p.x - p2.x) < 0.001) && (p.z >= minx && p.z <= maxx))
				{
					return true;
				}
			}

			//�����ߺ��߶εĽ���
			float x = (double)(p.z - p1.z) * (double)(p2.x - p1.x) / (double)(p2.z - p1.z) + p1.x;
			//ͳ�Ƴ����ཻ���߶εĸ���	
			if (x > p.x)//���ұ��ཻ
			{
				if (lastp2 == p1)
				{
					lastp1 = p1;
					lastp2 = p2;
					continue;
				}
				lastp1 = p1;
				lastp2 = p2;
				nCross++;
			}
			else//������ཻ
			{
				nCross1++;
			}
		}
	}
	//ż�����ڶ����֮��
	//�������ڶ����֮��
	if (nCross1 > 0)//������Ҳ�н���
	{
		if (nCross > 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	if (nCross1 == 0)//������û�н���
	{
		if ((nCross % 2) == 1)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}

vector<float> Fragment::getSupportPoints()
{
	vector<float> SupportPoints;
	//�Ѷ�Ӧ�ĵ�ͶӰ��ͬһ��ƽ����
	for (int i = eachSurfacePoints.size()-1; i >0; --i)//�������Ƿ��Ž��У��ȴ����һ���濪ʼ
	{
		for (int j = 0; j != eachSurfacePoints[i].size(); j+=3)//����һ��ÿһ���㶼���м���
		{
			//���ͶӰ������ƽ�淶Χ�ڣ�ͶӰ�����һ��֧�ŵ㣬
			//ʣ�µ�û�����ڷ�Χ�ڵ�ͶӰ�����ͶӰ����һ�㣬ֱ�����ڻ�׼
			point3d temp = { eachSurfacePoints[i][j],eachSurfacePoints[i][j + 1],eachSurfacePoints[i][j + 2] };
			int tag = 0;//�����жϵ��Ƿ���ǰ��������ϣ�0����û�У������Ǿ����ڻ�׼����
			
			if (SetInOrOut(temp,i-1) == true)//�����������һ���������Ǿ�Ĭ�ϲ���֧����
			{
				tag = 1;
				continue;
			}
			else
			{
				//����������������һ�������棬��͵�����ǻ�׼��
				//�������е�����ÿһ�㶼����������в���бȶ�
				for (int k = i - 2; k > 0; --k)
				{
					if (SetInOrOut(temp, k) == true)//�����������������ϣ������������ڲ������Ǿͼ�¼����
					{
						SupportPoints.push_back(temp.x);
						SupportPoints.push_back(temp.y);
						SupportPoints.push_back(temp.z);
						point3d aa = getShadow(temp, _begin+k*_each);
						SupportPoints.push_back(aa.x);
						SupportPoints.push_back(aa.y);
						SupportPoints.push_back(aa.z);
						tag = 1;
						break;
					}
				}
			}
			
			if (tag == 0)//���ڻ�׼����
			{
				SupportPoints.push_back(temp.x);
				SupportPoints.push_back(temp.y);
				SupportPoints.push_back(temp.z);
				point3d aa=getShadow(temp, _begin);
				SupportPoints.push_back(aa.x);
				SupportPoints.push_back(aa.y);
				SupportPoints.push_back(aa.z);
			}
		}
	}

	return SupportPoints;
}

vector<vector<float>> Fragment::getPoints(varray < varray<varray<point3d>>> points)
{
	vector<float> temp;//����һ����ʱ���飬�����洢һ�����ϵ��е�
	if (_degree == 1)//�������ֻ��һ�������ǾͰ��������������_begin��_end���м�
	{
		_degree = 2;
		float arvge = (_end - _begin) / 2 + _begin;
		_end = arvge;
		_begin = arvge;
		_each = 1;
	}
	//�������һ�������ǾͰ��������з�
	
		for (int n = _begin; n <= _end; n += _each)//�ж��ٸ�����ƽ�ͽ��ж��ٴ�
		{
			for (int i = 0; i != points.size(); ++i)//points.size()������ĸ�������6����Ϊ�������壬ÿ���涼����һ��
			{
				for (int j = 0; j != points[i].size(); ++j)//ÿһ�������ִ��ںܶ��С��ϸ�ֵ�Ԫ��ÿһ��ϸ�ֵ�Ԫ������һ��
				{
					for (int k = 0; k != points[i][j].size(); k += 4)//ÿһ��С��Ԫ�ֱ��󽻵�
					{
						////0  3///
						//	�� ��
						////1��2///
						//ÿһ���ı��ζ��������˳����ɵ�
						intersectionLinePlane(points[i][j][k], points[i][j][k + 3], temp, n);//����0��3����ɵ�һ����
						intersectionLinePlane(points[i][j][k + 1], points[i][j][k + 2], temp, n);//1��2�����һ����
						intersectionLinePlane(points[i][j][k], points[i][j][k + 1], temp, n);//0��1���
						intersectionLinePlane(points[i][j][k + 2], points[i][j][k + 3], temp, n);//2��3���
					}
				}
			}
			//Ϊ���Ժ�ļ��㷽�㣬��Ҫ�������ת��Ϊpoint3d���д洢
			vector<point3d> temp_p;
			for (int i = 0; i != temp.size(); i+=3)
			{
				point3d q = { temp[i],temp[i + 1],temp[i + 2] };
				temp_p.push_back(q);
			}
			eachSurfacePoints_point3d.push_back(temp_p);

			eachSurfacePoints.push_back(temp);
			temp.clear();//�����ʱ���飬���������һ�μ���ѭ��
		}
	return eachSurfacePoints;
}


bool Fragment::intersectionLinePlane(point3d p1, point3d p2,vector<float> & temp,  float D)//��ý���
{
	point3d pt, pt1, ptl, ptr;
	double num, den, n, distancept, distanceptl, distanceptr;
	pt.x = p2.x - p1.x;
	pt.y = p2.y - p1.y;
	pt.z = p2.z - p1.z;
	pt.U = p2.U - p1.U;
	pt.V = p2.V - p1.V;
	pt.W = p2.W - p1.W;
	pt.Gray = p2.Gray - p1.Gray;

	distancept = sqrt(pow(pt.x, 2) + pow(pt.y, 2) + pow(pt.z, 2));

	num = _A * p1.x + _B * p1.y + _C * p1.z + D; //��ĸ
	den = -_A * pt.x - _B * pt.y - _C * pt.z;   //����

	if (fabs(den) > 1e-5)
	{
		n = num / den;
		pt1.x = p1.x + n * pt.x;
		pt1.y = p1.y + n * pt.y;
		pt1.z = p1.z + n * pt.z;
		pt1.U = p1.U + n * pt.U;
		pt1.V = p1.V + n * pt.V;
		pt1.W = p1.W + n * pt.W;
		pt1.Gray = p1.Gray + n * pt.Gray;

		ptl.x = pt1.x - p1.x;
		ptl.y = pt1.y - p1.y;
		ptl.z = pt1.z - p1.z;

		ptr.x = pt1.x - p2.x;
		ptr.y = pt1.y - p2.y;
		ptr.z = pt1.z - p2.z;

		distanceptl = sqrt(pow(ptl.x, 2) + pow(ptl.y, 2) + pow(ptl.z, 2));
		distanceptr = sqrt(pow(ptr.x, 2) + pow(ptr.y, 2) + pow(ptr.z, 2));

		if ((distanceptl / distancept < 1 && distanceptr / distancept < 1) || distanceptl == 0 || distanceptr == 0)
		{
			if (pt1.U == 0 || pt1.U == 1 || pt1.V == 0 || pt1.V == 1 || pt1.W == 0 || pt1.W == 1)
			{
				pt1.boundary = true;
			}
			temp.push_back(pt1.x);
			temp.push_back(pt1.y);
			temp.push_back(pt1.z);
		}
	}
	return true;
}

//void Fragment::CIntersecPoints()
//{
//	for (int i = 0; i < Patchnum; i++)
//	{
//		for (int j = 0; j < PatchHex[i].size(); j++)
//		{
//			for (int k = 0; k < 4; k++)
//			{
//				intersectionLinePlane(PatchHex[i][j].m_ptidx[k], PatchHex[i][j].m_ptidx[k + 4], PatchHex[i][j]);
//			}
//			for (int l = 0; l < 7; l = l + 2)
//			{
//				intersectionLinePlane(PatchHex[i][j].m_ptidx[l], PatchHex[i][j].m_ptidx[l + 1], PatchHex[i][j]);
//			}
//			intersectionLinePlane(PatchHex[i][j].m_ptidx[0], PatchHex[i][j].m_ptidx[2], PatchHex[i][j]);
//			intersectionLinePlane(PatchHex[i][j].m_ptidx[1], PatchHex[i][j].m_ptidx[3], PatchHex[i][j]);
//			intersectionLinePlane(PatchHex[i][j].m_ptidx[5], PatchHex[i][j].m_ptidx[7], PatchHex[i][j]);
//			intersectionLinePlane(PatchHex[i][j].m_ptidx[4], PatchHex[i][j].m_ptidx[6], PatchHex[i][j]);
//		}
//	}
//	for (int i = 0; i < Patchnum; i++)
//	{
//		for (int j = 0; j < PatchHex[i].size(); j++)
//		{
//			if (PatchHex[i][j].Intersectpoints.size() == 1 || PatchHex[i][j].Intersectpoints.size() == 2)
//			{
//				PatchHex[i][j].Intersectpoints.clear();
//			}
//		}
//	}
//}
//
//point3d Fragment::GetCubicPoint(int Unum, int Vnum, int Wnum, int uid, int vid, int wid, int uorder, int vorder, int worder, varray<double> Nu, varray<double> Nv, varray<double> Nw, varray<point3d> &bpt, double u, double v, double w)
//{
//	point3d val;
//	int num = 0;
//	varray<point3d> bracepoint;
//	int ii = 0, jj = 0, kk = 0, flag = 0;
//
//	for (int k = 0; k <= worder; k++)
//	{
//		kk = wid - worder + k;
//		for (int j = 0; j <= vorder; j++)
//		{
//			jj = vid - vorder + j;
//			for (int i = 0; i <= uorder; i++)
//			{
//				ii = uid - uorder + i;
//				val += alf_mult_point(Nu[i] * Nv[j] * Nw[k], bpt[kk*(Unum + 1)*(Vnum + 1) + jj * (Unum + 1) + ii]);
//			}
//		}
//	}
//
//	val.u = u; val.v = v; val.w = w;
//	return val;
//}
//
//void Fragment::GetAllCell(int Unum, int Vnum, int Wnum, int Uorder, int Vorder, int Worder, varray<double> &uknots, varray<double> &vknots, varray<double> &wknots, varray<point3d> &bpts)
//{
//	Hex ahex;
//	NurbsSurface temp;
//	double degree = 20.0f;//ϸ�־���
//	int k, j, i;
//	int u_id[2], v_id[2], w_id[2];//r,s,q  degeree
//	double u[2], v[2], w[2];
//	varray<double> Nu[2], Nv[2], Nw[2];
//	Nu[0].resize(Uorder + 1);
//	Nv[0].resize(Vorder + 1);
//	Nw[0].resize(Worder + 1);
//	Nu[1].resize(Uorder + 1);
//	Nv[1].resize(Vorder + 1);
//	Nw[1].resize(Worder + 1);
//	m_allHexs.clear();
//	for (k = 0; k < degree; k++)//uv��Ȳ���
//	{
//		w[0] = double(k) / degree;
//		w[1] = double(k + 1) / degree;
//		w_id[0] = Findnum(Worder, Wnum + 1, w[0], wknots);    //�ڵ�λ��
//		w_id[1] = Findnum(Worder, Wnum + 1, w[1], wknots);
//		temp.BasisFuns(w[0], w_id[0], Worder, wknots, Nw[0]);  //������ 
//		temp.BasisFuns(w[1], w_id[1], Worder, wknots, Nw[1]);  //������ 
//		for (j = 0; j < degree; j++)//vnum
//		{
//			v[0] = double(j) / degree;
//			v[1] = double(j + 1) / degree;
//			v_id[0] = Findnum(Vorder, Vnum + 1, v[0], vknots);
//			v_id[1] = Findnum(Vorder, Vnum + 1, v[1], vknots);
//			temp.BasisFuns(v[0], v_id[0], Vorder, vknots, Nv[0]);
//			temp.BasisFuns(v[1], v_id[1], Vorder, vknots, Nv[1]);
//			for (i = 0; i < degree; i++)
//			{
//				int p = 0;
//				u[0] = double(i) / degree;
//				u[1] = double(i + 1) / degree;
//				u_id[0] = Findnum(Uorder, Unum + 1, u[0], uknots);
//				u_id[1] = Findnum(Uorder, Unum + 1, u[1], uknots);
//				temp.BasisFuns(u[0], u_id[0], Uorder, uknots, Nu[0]);
//				temp.BasisFuns(u[1], u_id[1], Uorder, uknots, Nu[1]);
//				for (int l = 0; l < 2; l++) //u����
//				{
//					for (int m = 0; m < 2; m++) //v����
//					{
//						for (int n = 0; n < 2; n++) //w����
//						{
//							ahex.m_ptidx[p] = GetCubicPoint(Unum, Vnum, Wnum, u_id[l], v_id[m], w_id[n], Uorder, Vorder, Worder, Nu[l], Nv[m], Nw[n], bpts, u[l], v[m], w[n]);
//							//ahex.m_ptidx[p].Gray = GetGrayscale(Unum, Vnum, Wnum, u_id[l], v_id[m], w_id[n], Uorder, Vorder, Worder, Nu[l], Nv[m], Nw[n], bpts);
//							p++;
//						}
//					}
//				}
//				m_allHexs.push_back(ahex);
//			}
//		}
//	}
//	PatchHex.push_back(m_allHexs);
//	for (int i = 0; i < 2; i++)
//	{
//		Nu[i].clear();
//		Nv[i].clear();
//		Nw[i].clear();
//	}
//}
