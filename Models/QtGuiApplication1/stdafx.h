
// stdafx.h : ��׼ϵͳ�����ļ��İ����ļ���
// ���Ǿ���ʹ�õ��������ĵ�
// �ض�����Ŀ�İ����ļ�

#pragma once

#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN            // �� Windows ͷ���ų�����ʹ�õ�����
#endif


#include "lib3mf_implicit.hpp"
#include <thread>
#include <future>
#include "varray.h"
#include "XVec.h"
#include "pointXd.h"
#include <Windows.h>
#include <iostream>
#include <vector>
using namespace base;

struct threadParamSpline {
	varray<varray<Vec3>> l2;
	varray<varray<Vec3>> p3d2;
};



struct threadParamSplineVOL {
	varray<varray<varray<Vec3>>> q3;
	varray<varray<varray<Vec3>>> l3;
};
struct threadParam {
	varray<varray<point3d>> l2;
	varray<varray<point3d>> p3d2;
};

struct threadParamVOL {
	varray<varray<varray<point3d>>> q3;
	varray<varray<varray<point3d>>> l3;
};



