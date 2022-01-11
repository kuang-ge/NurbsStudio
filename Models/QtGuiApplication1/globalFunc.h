//ȫ�ֺ���
#ifndef _globalFunc_H_
#define _globalFunc_H_

#include "varray.h"
#include "pointXd.h"
#include "../lib/eigen/Eigen/Dense"
#include <string>
#include "XVec.h"


using Eigen::Vector4d;
using Eigen::Matrix4d;
using std::string;

using namespace base;
/////////////////////////////////////////////���ı�    2018///////////////////////////////

template<class _T>
varray<_T> transVectorToVarray(const std::vector<_T> vector)
{
	varray<_T> varray;
	for (auto&i : vector) {
		varray.push_back(i);
	}
	return varray;
}


/*ά�Ƚ�һ
hightDimVarray����άvarray
lowDimVarray����άvarray
clear���Ƿ������άvarrayԭ������
*/
template<class _T>
void ReduceVarrayDim(const varray<varray<_T>>& hightDimVarray, varray<_T>& lowDimVarray, const bool clear)
{
	if (clear)
		lowDimVarray.clear();
	if (hightDimVarray.size() == 0)return;

	for (int i = 0; i < hightDimVarray.size(); ++i)
		for (int j = 0; j < hightDimVarray[i].size(); ++j)
			lowDimVarray.push_back(hightDimVarray[i][j]);
}

//��ά����ת��
template<class _T>
void Transpose(varray<varray<_T>>& Varr)
{
	varray<varray<_T>> tvarr;
	tvarr.resize(Varr[0].size());
	for (int i = 0; i < Varr.size(); ++i)
		for (int j = 0; j < Varr[i].size(); ++j)
			tvarr[j].push_back(Varr[i][j]);
	Varr = tvarr;
}

/*����ʽϵ��
Bin(x,y)=x!/(y!*(x-y)!)
*/
double Bin(const int x, const int y);

//��ʾ��淶����[a��b]
void Normalization(varray<point3d>& pts, const double a, const double b);

//��ʾ��淶����[a��b]
void Normalization(varray<varray<point3d>>& pts, const double a, const double b);

//����㼯����
void GetCenter(varray<point3d>& pts, point3d &center);

//�㼯����ƽ��
//mode=1 vecΪƽ���յ�,��ʱ��Ҫ����㼯����center
//mode=0 vecΪƽ������
template<class T>
void TransPts(varray<T> &pts, const T &vec, bool mode, const T &center = T())
{
	if (mode)
	{
		for (int i = 0; i < pts.size(); ++i)
			pts[i] = pts[i] + (vec - center);
	}
	else
	{
		for (int i = 0; i < pts.size(); ++i)
			pts[i] = pts[i] + vec;
	}
}

//���ʽת��
//point4dתVector4d
Vector4d P4dToV4d(const point4d P4d);

//���ʽת��
//Vector4dתpoint4d
point4d V4dToP4d(const Vector4d V4d);

Vector4d P4dToV4d(const Vec4 P4d);


//Vector4d��λ��
void Vector4dToUnit(Vector4d &v);

//Vector4d���
void AcrossB(const Vector4d& a, const Vector4d& b, Vector4d& c);

//����н�
double angleAB(const Vector4d& a, const Vector4d& b);

//�������߳���,LΪ���������
double CurveLength(const varray<point3d>& L);

//�ж�a�Ƿ���S��
template<class _T>
bool IsInSet(const _T a, const varray<_T>& S)
{
	for (int i = 0; i < S.size(); ++i)
	{
		if (a == S[i])
			return true;
	}
	return false;
}

//��������ϵobj��target�Ĺ��ɾ���
//Coordinates��X�ᣬY�ᣬZ��, ����ԭ��
void CalCoordinatesTransMat(const varray<point3d>& objCoordinates, const varray<point3d>& targetCoordinates, Matrix4d& mat);

//����ʸ���任����
void CalTransMat(const Vector4d& targetPts, const Vector4d& targetVec,
	const Vector4d& objPts, const Vector4d& objVec, Matrix4d& mat);
void CalTransMat(const point4d& targetPts, const point4d& targetVec,
	const point4d& objPts, const point4d& objVec, Matrix4d& mat);
void CalTransMat(const Vec4& targetPts, const Vec4& targetVec,
	const Vec4& objPts, const Vec4& objVec, Matrix4d& mat);

//����任
void TransByMat(varray<Vector4d>& varr, const Matrix4d& mat);
void TransByMat(varray<point4d>& varr, const Matrix4d& mat);
void TransByMat(varray<Vec4>& varr, const Matrix4d& mat);

//ȡ���ź���
int Sgn(double a);

/*����������㼰����
beginIdx:����Ϊ���Ԫ�صĵ�ǰ�±�
dir���Ƿ񰴵�ǰ˳�򴢴�
*/
template<class _T>
void SetBeginDir(varray<_T>& varr, const size_t beginIdx, const bool dir = true)
{
	varray<_T> res;
	varray<_T>::iterator it0, it;
	it = it0 = varr.begin() + beginIdx;
	do
	{
		res.push_back(*it);
		if (dir)
		{
			++it;
			if (it == varr.end())
				it = varr.begin();
		}
		else
		{
			--it;
			if (it == varr.begin())
			{
				res.push_back(*it);
				it = varr.end() - 1;
			}
		}
	} while (it != it0);
	varr = res;
}

//���Լ��
int GCD(int a, int b);

//����
template<class _T>
void KnotUnify(const varray<_T>& KnotsA, const varray<_T>& KnotsB, varray<_T>& NewKnots)
{
	NewKnots.clear();
	int Anum = 0, Bnum = 0;
	// ���Anum��Bnum��һ���������˾�ֹͣ
	while (Anum < KnotsA.size() && Bnum < KnotsB.size())
	{
		if (KnotsA[Anum] < KnotsB[Bnum])
		{
			NewKnots.push_back(KnotsA[Anum++]);
			continue;
		}
		else if (KnotsA[Anum] > KnotsB[Bnum])
		{
			NewKnots.push_back(KnotsB[Bnum++]);
			continue;
		}
		else
		{
			NewKnots.push_back(KnotsB[Bnum++]);
			Anum++;
		}
	}
	// ��� Anum< KnotsA.size()��˵��Bnum����ӽ�����NewKnots�ˣ�ֻʣ�±�Bnum���Anum��������
	while (Anum < KnotsA.size())
	{
		NewKnots.push_back(KnotsA[Anum++]);
	}
	while (Bnum < KnotsB.size())
	{
		NewKnots.push_back(KnotsB[Bnum++]);
	}
}
//�
//diffKnots=KnotsL-KnotsS
//KnotsL��С��С��KnotsS
template<class _T>
void KnotsDiff(const varray<_T>& KnotsL, const varray<_T>& KnotsS, varray<_T>& diffKnots)
{
	diffKnots.clear();
	if (KnotsL.size() <= KnotsS.size())
		return;

	int Ln = 0, Sn = 0;
	while (Ln < KnotsL.size() && Sn < KnotsS.size())
	{
		if (KnotsS[Sn] != KnotsL[Ln])
			diffKnots.push_back(KnotsL[Ln]);
		else
			++Sn;
		++Ln;
	}
	while (Ln < KnotsL.size())
	{
		diffKnots.push_back(KnotsL[Ln++]);
	}
}

//��pts������v1ͶӰ��ƽ��(v,p)
//v1Ĭ��Ϊ0��������ʱΪ��ͶӰ���൱��v1 = -v
//ƽ�з���0,���򷵻�1
int Project2Plane(const point3d& v, const point3d& p, const point3d& pts, point3d& res, const point3d& v1 = point3d());

//�㼯Cpts������v1ͶӰ��ƽ��(v,p)
//v1Ĭ��Ϊ0��������ʱΪ��ͶӰ���൱��v1 = -v
int Project2Plane(const point3d & v, const point3d & p,
	const varray<point4d>& Cpts, varray<point4d>& ProPts, const point3d & v1 = point3d());

#endif