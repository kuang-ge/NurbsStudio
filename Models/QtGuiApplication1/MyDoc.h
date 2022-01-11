#pragma once

#include "CNurbs.h"
#include "Fragment.h"
#include <qstring.h>
#include <qvector4d.h>
#include <qvector3d.h>

class MyDoc
{
public:
	//��������
	void cleanUp()
	{
		m_ShowData.clear();
		m_ShowLines.clear();
	}
	

	//������ʾ����
	void clearScrn()
	{
		VAO_L_ARRAY.clear();
		VAO_D_ARRAY.clear();
		VAO_CP_ARRAY.clear();
		VAO_FRAGMENT_ARRAY.clear();
		VAO_SUPPORT_ARRAY.clear();
		EACH_VAO_L_POINT_NUMBER.clear();
		EACH_VAO_CP_POINT_NUMBER.clear();
		EACH_VAO_D_POINT_NUMBER.clear();
		EACH_VAO_F_NUMBER.clear();
		EACH_VAO_SUPPORT_NUMBER.clear();
	}

	//��ȡ����ģʽָ��
	typedef std::shared_ptr<MyDoc> Ptr;
	MyDoc(const MyDoc&) = delete;
	MyDoc& operator=(const MyDoc&) = delete;
	static Ptr getInstance()
	{
		if (m_instance_ptr == nullptr)
		{
			std::lock_guard<std::mutex> lk(mutex);
			if (m_instance_ptr == nullptr)
			{
				m_instance_ptr = std::shared_ptr<MyDoc>(new MyDoc);
			}
		}
		return m_instance_ptr;
	}
	~MyDoc() {}
private:
	static Ptr m_instance_ptr;
	static std::mutex mutex;
	MyDoc() {}
public:
	
	QString path;//ģ��·��
	int modelType;//ģ�͵����ͣ�nurbsVOL��nurbsSurface�ȣ�

	//��ʼϸ�ֶ�
	int m_Unum = 20;
	int m_Vnum = 20;
	int m_Wnum = 20;

	int m_tolNum;//��Ԫ����
	int m_NUMPic = 6;//����Ƭ����Ƭ��Ŀ
	int m_CellIdx = -1;//��Ƭģ�͵ı��
	float m_dotSZ = 10.0;//��ʾʱ����Ĵ�С
	float m_CPTSdotSZ = 10.0;//��ʾʱ�����Ƶ�Ĵ�С
	float _lineWidth = 3.0;//�߿�
	bool m_showCtrlPts=false;//�Ƿ���ʾ���Ƶ�
	bool m_showCoordinates=true;//�Ƿ���ʾ����ϵ
	bool m_restView = false;//�Ƿ������ӽ�
	bool m_Fragment = false;//�Ƿ���Ҫ��Ƭ
	int mode = 7;

	int THREADNUMBERS = 8;//�̴߳�С������

	varray<varray<varray<point3d>>> m_ShowData;//3D��ʾ����
	varray<varray<varray<point3d>>> m_ShowLines;//������ʾ����
	varray<varray<point4d>> m_ShowData4D;//4D���Ƶ���ʾ����

	varray<varray<varray<Vec3>>> m_DisplayData;//3D��ʾ����
	varray<varray<varray<Vec3>>> m_DisplayLines;//������ʾ����
	varray<varray<Vec4>> m_DisplayData4D;//4D���Ƶ���ʾ����
	
	std::vector<float> showDta_vertices;//�����滭��������ۺ�Ϊһ��һά����
	std::vector<float> showLines_vertices;//�����滭�ߵ������ۺ�Ϊһ��һά����
	std::vector<float> showCp_vertices;//��������Ƶ��ۺ�Ϊһ��һά����
	std::vector<float> showFp_vertices;//��Ƭ��
	
	std::vector<unsigned int> VAO_L_ARRAY, VAO_D_ARRAY, VAO_CP_ARRAY, VAO_FRAGMENT_ARRAY, VAO_SUPPORT_ARRAY;//���涥�����黺�����ľ��
	std::vector<int> EACH_VAO_L_POINT_NUMBER, EACH_VAO_CP_POINT_NUMBER, EACH_VAO_D_POINT_NUMBER, EACH_VAO_F_NUMBER, EACH_VAO_SUPPORT_NUMBER;//����ÿһ���ߵĸ���

	std::vector<unsigned> VAO_D;//��Ⱦ�Ǿ���ģ�͵ľ��
	std::vector<int> VAO_D_SIZE;
	int m_GLMode = 7;//0=POINTS,3=LINE,7=FILL

	float m_fragmentA, m_fragmentB, m_fragmentC;//��Ƭ����
	int m_fragmentDegree;//���ٸ�����
	float m_fragmentBegin, m_fragmentEnd;//��ʼ�ͽ�����λ��
	
	float _mouseSenstive=0.01;//���������
	float _keyBordSensitive=10;//����������

	QVector3D _bordLineColor = { 0.1f,0.1f,0.1f };//�߽�����ɫ Ĭ�Ϻ�ɫ
	QVector3D _surfaceColor = { 0.0f,1.0f,0.0f };//��Ƭ��ɫ Ĭ����ɫ
	QVector3D _controlPoinrColor = { 1.0f,0.0f,0.0f };//���Ƶ���ɫ Ĭ�Ϻ�ɫ
	QVector3D _sliceSurfaceColor = { 1.0f,0.0f,0.0f};//��Ƭ����ɫ Ĭ�Ϻ�ɫ
	QVector3D _sliceSuppotrColor = { 0.0f,0.0f,1.0f };//��Ƭ֧��ɫ��ɫ Ĭ����ɫ
	QVector4D _lightColor = { 1.0f,1.0f,1.0f,1.0f };

	//������ĳ�ʼλ��
	float cameraX = 0;
	float cameraY = 0;
	float cameraZ = 80;
};
