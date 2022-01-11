#pragma once
#pragma once
#include"varray.h"
#include"XVec.h"
#include "SplineSurface.h"
#include "SplineVolume.h"
using namespace base;
//�������۷���
//��Դ The Verdict Geometric Quality Library
//���д󲿷ַ���������paraview��ʵ��
//������ϣ��������Щָ������ۺ�����
//��Ϊ��һָ���������ײ���Ƭ��Ľ��
//�ο�Huang jin�� Evaluating Hex-mesh Quality Metrics via Correlation Analysis


class BaseInfo {
public:
	static BaseInfo* getInstancePtr() {
		if (instance == nullptr)
			instance = new BaseInfo();
		return instance;
	}
	BaseInfo(const BaseInfo&) = delete;
	BaseInfo& operator = (const BaseInfo&) = delete;
	varray<Vec3> pts;				//�ĸ��ǵ�
	varray<Vec3> Li;				//��Ԫ������
	varray<double> mLi;				//��Ԫ�߳���
	double Lmin;					//��С��
	double Lmax;					//����
	varray<Vec3> Di;				//�ԽǱ�
	varray<double> mDi;				//�ԽǱ߳���
	double Dmax;					//���ԽǱ�
	double Dmin;					//��С�ԽǱ�
	varray<Vec3> Xi;				//������
	varray<Vec3> Ni;				//������
	varray<Vec3> ni;				//��λ������
	Vec3 Nc;						//��������
	Vec3 nc;						//���ĵ�λ����
	varray<double> ai;				//�������
	varray<varray<Vec3>> Ai;		//3x3����
	varray<double> A2;				//Aiģ��ƽ�� 
	varray<double> adjA2;			
	varray<varray<Vec3>> Ai_;
	varray<double> alphi;
	varray<double> alphi_;
	///////////////�ı���ָ�����
	double getArea();					//�������
	double getAspectRatio();			//�����߱�
	double getCondition();
	double getDistortion();				//����Ť����
	double getEdgeRatio();				//�������
	double getJacobian();				//����ſɱ�ֵ
	void InitSufInfo(const varray<Vec4>& pts);
	void InitVolInfo(const varray<Vec4>& pts);
	void clearData();
	////////////�ı���ָ��
	double Area;					//���
	double AspectRatio;				//�ݺ��
	double Condition;
	double Distortion;				//Ť����
	double EdgeRatio;				//����
	double Jacobian;				//�ſɱ�
	////////////������ָ��
	double Diagonal;			//�Խ���[0.65,1]
	double Volume;				//���[0,+OO]
	double Stretch;				//��չ[0.25,1]
	double ScaledJacobian;		//�߽��ſɱ�[0.5,1]
	double Oddy;				//�ſɱȵ���һ�ֶ���[0,0.5]
	////////////������ָ�����
	double getDiagonal();
	double getVolume();
	double getStretch();
	double getScaledJacobian();

private:
	BaseInfo() {}

	static BaseInfo* instance;

};

//���ڼ��ģ����������
class BaseQuality
{
public:
	BaseQuality(int n);
	BaseQuality();
	~BaseQuality();
	varray<varray<varray<Vec4>>> vecs;		//ģ���ϵĵȲε�Ԫ	
	void virtual calQuality() = 0;			//ģ����������
	void virtual getSamplePts() = 0;		//���ģ�͵Ȳε�Ԫ
	int segmentNum;					//Ĭ�ϵ�ȡ������Ϊ10
};

class SufQuality : public BaseQuality {
public:
	SufQuality(varray<SplineSurface>& sf, int n = 5);
	~SufQuality();
	/////////////////////////////////
	//�û��������ֻ����������������//
	/////////////////////////////////
	void calQuality();
private:
	void getSamplePts();
	bool getBaseInfo(const varray<Vec4>& pt);
	bool compareData();
	void clearAll();
public:
	double minArea = INT_MAX;		//���������0
	double minJacobian = INT_MAX; //���������0
	double minEdgeRatio = INT_MAX;//>=1
	double maxEdgeRatio = INT_MIN;//<=1.3
	double minCondition = INT_MAX;//>=1
	double maxCondition = INT_MIN;//<=4
	double minAspectRatio = INT_MAX;//>=1
	double maxAspectRatio = INT_MIN;//<=1.3
	double minDistortion = INT_MAX;//>=0.5
	double maxDistortion = INT_MIN;//<=1


private:
	varray<SplineSurface> sf;		//����Ĵ���������
};

class VolQuality : public BaseQuality {
public:
	VolQuality(varray<SplineVolume>& sv, int n=10);
	~VolQuality();
	void calQuality();
	void getSamplePts();
	bool getBaseInfo(const varray<Vec4>& pt);
	bool compareData();
	void clearAll();

	double minDiagonal = INT_MAX;			//�Խ���[0.65,1]
	double maxDiagonal = INT_MIN;
	double minVolume = INT_MAX;				//���[0,+OO] 
	double maxVolume = INT_MIN;
	double minStretch = INT_MAX;				//��չ[0.25,1]
	double maxStretch = INT_MIN;
	double minScaledJacobian = INT_MAX;		//�߽��ſɱ�[0,1]
	double maxScaledJacobian = INT_MIN;		


private:
	varray<SplineVolume> sv;



};


