#pragma once

#include "Nurbs.h"
#include <list>
#include <stack>
#include <string>
#include <fstream>

using std::vector;
using std::unique_ptr;
using std::string;
using std::ifstream;
using std::ofstream;

/*�������ݽṹͨ����֮��Ĺ�ϵ���ӣ�������ɾ��ĵ����ɲ���*/
namespace YN {


	//ö�����飬ʶ����
	enum Direct { u0v0, u1v1, v0w0, v1w1, u0w0, u1w1 };
	struct Cube;

	//�������ݽṹ����
	struct CurveLine
	{
		YN::NurbsLine _lineData;
	};

	//�������ݽṹ����
	struct Surface
	{
		//������Ϣ����������Ϣ
		YN::NurbsSurface _surfaceData;//�������
		vector<int> _surfaceGloblControlPointID;//��Ŀ��Ƶ�ȫ�ֱ��
		//������Ϣ
		Surface * _nextSurface = NULL;//�ڽ���
		bool _haveNextSuf;//�Ƿ����ڽ���
		Cube * _pointOther = NULL;//�ڽ�Ƭ
		Cube * _pointThis = NULL;//��Ƭ
		int _pointThisNumber;//���ʶ�����
	};

	//�������ݽṹ����
	struct Cube
	{
		YN::NurbsVol _Vol_Data;
		vector<int> _volGloblControlPointID;//Ƭ�ڿ��Ƶ�ȫ�ֱ��
		Surface * _u0v0 = NULL;//0  ����Ϊ���ʶ�����,��Ӧö������Direct
		Surface * _u1v1 = NULL;//1
		Surface * _v0w0 = NULL;//2
		Surface * _v1w1 = NULL;//3
		Surface * _u0w0 = NULL;//4
		Surface * _u1w1 = NULL;//5
		bool _isOrdered = false;//�Ƿ񾭹�����
		bool _bordFirstSearch = true;
	};

	//�������ݽṹ
	class WingStruct
	{
	public:
		WingStruct();
		//��Ȩ��tag:false ��Ȩ��tag:true
		WingStruct(string path, bool tag);
		WingStruct(vector<NurbsVol> vols);
		~WingStruct();

		//��
		bool addVolToSurface(point4d target, YN::NurbsLine line, YN::NurbsSurface surface);
		//ɾ
		bool removeVol(point4d target);
		//�飬���ٲ�ѯ������Ŀ��㣬��ʼƬ�������ѯ·��Ƭ
		vector<Cube *> quickSearch(point4d, Cube *);
		//��
		void Modify(Cube * target, int direct, int u, int v, int w);
		//��������ı���
		vector<YN::NurbsSurface> getModelSurface();
		//�ļ��Ķ�ȡ,����tag���Ϊtrue��ʾ�ļ�����Ȩ�أ�tagΪfalse��ʾ�ļ�û��Ȩ��
		vector<Cube *> ReadVolumeTxt(string path,bool tag);
		//�ļ����������ֵ��ʩ����ģ���ϵ�����Լ��ʾ���
		vector<point4d>  OutputParaVolumeDataTxt(string filename, string Cp);
		//�ļ����ΪVTK
		void  OutputParaVolumeDataVTK(string filename);
	public:
		Cube * _root;//��ָ��
	private:
		//���ݽṹ���ܸ���
		WingStruct(const WingStruct&) = delete;
		WingStruct& operator=(const WingStruct&) = delete;
		//��ȡ��Ϣ
		void abstractInfo(vector<Cube *>&);
		//����cube����������Ϣ
		vector<Surface *> upDateSurface(vector<Cube *>&);
		vector<Surface *> upDateCubeSurface(Cube *);
		//NurbsVolת��ΪCube
		vector<Cube *> transNurbsVoltoCube(vector<YN::NurbsVol>);
		//����������Ϣ
		bool setTwinPatch(vector<Cube *>&, vector<Surface *>&);
		//���ò�������Ϣ
		bool Order(vector<Cube *>&);
		//��ȱ���
		vector<Cube *> bordFirstSearch();
		//��ȱ���,������ʼƬ������
		vector<Cube *> deepFirstSearch(Cube *, int);
		//����Ƭ����Ϊ��ͬUVW
		bool setNearCubeWithSameUVW(Cube *, Cube *, int index);
		//���ȫ�ֿ��Ƶ��id
		void setGloblControlPointId();
		//�жϵ��Ƿ���Ƭ�ڲ�,���ڲ����ر�ţ������ڲ�����-1
		int isPointInSurface(point4d, Cube *);
	};

};