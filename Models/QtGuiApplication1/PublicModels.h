#pragma once
#include "FeatureNetwork.h"

#include "SplineVolume.h"
#include "NurbsTrans.h"
#include "Option.h"

//�����ݽṹͷ�ļ�
#include "WingStruct.h"
#include "Nurbs.h"

//using namespace YN;

//����Ƕ�Բ��
class Circle0 {
public:
	Circle0() {};

	//����뾶�ͽǶ�
	Circle0(double r, double a) :r(r), a(a) {};

	Vec4 getP1() {
		Vec4 p01;
		p01.x = v0.x * l0;
		p01.y = v0.y * l0;
		p01.z = v0.z * l0;
		p01.w = 1;
		return p01;
	}
	Vec4 getP2() {
		Vec4 p02;
		p02.x = v1.x * l1;
		p02.y = v1.y * l1;
		p02.z = v1.z * l1;
		p02.w = w;
		return p02;
	}
	Vec4 getP3() {
		Vec4 p03;
		p03.x = v2.x * l2;
		p03.y = v2.y * l2;
		p03.z = v2.z * l2;
		p03.w = 1;
		return p03;
	}

	Spline getCircle() {
		Spline SL;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		SL.m_Knots = knots;
		SL.m_Degree = 2;
		p01.x = v0.x * l0;
		p01.y = v0.y * l0;
		p01.z = v0.z * l0;
		p01.w = 1;
		p02.x = v1.x * l1;
		p02.y = v1.y * l1;
		p02.z = v1.z * l1;
		p02.w = w;
		p03.x = v2.x * l2;
		p03.y = v2.y * l2;
		p03.z = v2.z * l2;
		p03.w = 1;
		SL.m_CtrlPts.push_back(p01);
		SL.m_CtrlPts.push_back(p02);
		SL.m_CtrlPts.push_back(p03);

		return SL;
	}

public:
	double r;//�뾶
	double a;//�Ƕ�
	Vec4 p01, p02, p03;//���Ƶ�
private:
	double a1 = a / 2;
	double a2 = a / 2;

	double l0 = r;
	double l1 = r / cos(a1);
	double l2 = r;
	double w = cos(a / 2);//Ȩ��
	Vec3 v1 = { 0,1,0 };//+y�ᵥλ����
	Vec3 v0 = v1.RotateZ(-a / 2);
	Vec3 v2 = v1.RotateZ(a / 2);
};

//����Ƕ�Բ��
class Annulus0 {
public:
	double r1;//�뾶
	double r2;
	double a;//�Ƕ�

	Annulus0() {};

	Annulus0(double r1, double r2, double a) :r1(r1), r2(r2), a(a) {}

	SplineSurface getSurface() {
		varray<Spline> SL1;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		SL1.resize(4);
		for (int i = 0;i < 4;i++) {
			SL1[i].m_Degree = 2;
			SL1[i].m_Knots = knots;
		}

		Circle0 c1(r1, a), c2(r2, a);
		/*Spline s1, s2;
		s1 = c1.getCircle();
		s2 = c2.getCircle();*/

		Vec4 p01 = c1.getP1();
		Vec4 p02 = c2.getP1();

		Vec4 p03 = c1.getP3();
		Vec4 p04 = c2.getP3();

		SL1[0].m_CtrlPts.push_back(c1.getP1());
		SL1[0].m_CtrlPts.push_back(c1.getP2());
		SL1[0].m_CtrlPts.push_back(c1.getP3());

		SL1[1].m_CtrlPts.push_back(p01);
		SL1[1].m_CtrlPts.push_back((p02 + p01) / 2);
		SL1[1].m_CtrlPts.push_back(p02);

		SL1[2].m_CtrlPts.push_back(c2.getP1());
		SL1[2].m_CtrlPts.push_back(c2.getP2());
		SL1[2].m_CtrlPts.push_back(c2.getP3());

		SL1[3].m_CtrlPts.push_back(p03);
		SL1[3].m_CtrlPts.push_back((p03 + p04) / 2);
		SL1[3].m_CtrlPts.push_back(p04);

		SplineSurface ss1;
		ss1.CoonsInterpolate(SL1);

		return ss1;
	}
};

//������
class Rectangle0 {
public:
	double l;//��
	double h;//��
	Vec4 v1, v2, v3, v4;

	Rectangle0() {};
	Rectangle0(double l,double h) :l(l),h(h){
		Vec4 p1 = { -l / 2,-h / 2,0,1 };
		Vec4 p2 = { l / 2,-h / 2,0,1 };
		Vec4 p3 = { -l / 2,h / 2,0,1 };
		Vec4 p4 = { l / 2,h / 2,0,1 };

		this->v1 = p1;
		this->v2 = p2;
		this->v3 = p3;
		this->v4 = p4;
	};

	//�ɳ��Ϳ�õ�������
	varray<Spline> getRectangle() {
		varray<Spline> SL;
		SL.resize(4);
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		for (int i = 0;i < SL.size();i++) {
			SL[i].m_Knots = knots;
			SL[i].m_Degree = 2;
		}
		SL[0].m_CtrlPts.push_back(v1);
		SL[0].m_CtrlPts.push_back((v1 + v2) / 2);
		SL[0].m_CtrlPts.push_back(v2);

		SL[1].m_CtrlPts.push_back(v1);
		SL[1].m_CtrlPts.push_back((v1 + v3) / 2);
		SL[1].m_CtrlPts.push_back(v3);

		SL[2].m_CtrlPts.push_back(v3);
		SL[2].m_CtrlPts.push_back((v3 + v4) / 2);
		SL[2].m_CtrlPts.push_back(v4);

		SL[3].m_CtrlPts.push_back(v2);
		SL[3].m_CtrlPts.push_back((v4 + v2) / 2);
		SL[3].m_CtrlPts.push_back(v4);

		return SL;
	}

	//���ĵ��γɳ����Σ��ʷ��㷨��ʹ�ã���һ�����򴫲Σ�
	varray<Spline> getRectangle(Vec4 v1, Vec4 v2, Vec4 v3, Vec4 v4) {
		varray<Spline> SL;
		SL.resize(4);
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		for (int i = 0;i < SL.size();i++) {
			SL[i].m_Knots = knots;
			SL[i].m_Degree = 2;
		}
		SL[0].m_CtrlPts.push_back(v1);
		SL[0].m_CtrlPts.push_back((v1 + v2) / 2);
		SL[0].m_CtrlPts.push_back(v2);

		SL[1].m_CtrlPts.push_back(v2);
		SL[1].m_CtrlPts.push_back((v2 + v3) / 2);
		SL[1].m_CtrlPts.push_back(v3);

		SL[2].m_CtrlPts.push_back(v3);
		SL[2].m_CtrlPts.push_back((v3 + v4) / 2);
		SL[2].m_CtrlPts.push_back(v4);

		SL[3].m_CtrlPts.push_back(v4);
		SL[3].m_CtrlPts.push_back((v4 + v1) / 2);
		SL[3].m_CtrlPts.push_back(v1);

		return SL;
	}

	SplineSurface getSurface() {
		varray<Spline> sl;
		sl = getRectangle();
		SplineSurface ss;
		ss.CoonsInterpolate(sl);
		return ss;
	}

};

//��(�����ı��ʷ�)
class Spline0 {
public:
	Vec4 v1, v2;

	Spline0() {};
	Spline0(Vec4 v1,Vec4 v2) {
		this->v1 = v1;
		this->v2 = v2;
	};

	Spline getSpline() {
		Spline SL;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		SL.m_Knots = knots;
		SL.m_Degree = 2;

		SL.m_CtrlPts.push_back(v1);
		SL.m_CtrlPts.push_back((v1 + v2) / 2);
		SL.m_CtrlPts.push_back(v2);

		return SL;
	}

	//�߶Σ����������������꣩
	Spline getSpline(Vec4 v1,Vec4 v2) {
		Spline SL;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		SL.m_Knots = knots;
		SL.m_Degree = 2;

		SL.m_CtrlPts.push_back(v1);
		SL.m_CtrlPts.push_back((v1 + v2) / 2);
		SL.m_CtrlPts.push_back(v2);

		return SL;
	}

	//Բ�������������������� Ĭ�Ͼ�ʮ�ȣ�
	Spline getArcSpline(Vec4 v1, Vec4 v2, Vec4 v3) {
		Spline SL;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		SL.m_Knots = knots;
		SL.m_Degree = 2;

		SL.m_CtrlPts.push_back(v1);
		SL.m_CtrlPts.push_back(v2);
		SL.m_CtrlPts.push_back(v3);

		return SL;
	}

	//Բ��������Բ�뾶��
	varray<Spline> arcSplines(double r) {
		Model_Solution m;
		varray<Spline> SL;
		double w = cos(PI / 4);
		Spline sl;
		Vec4 v1 = { -r,0,0,1 };
		Vec4 v2 = { -r,r,0,w };
		Vec4 v3 = { 0,r,0,1 };
		sl = getArcSpline(v1, v2, v3);
		SL.push_back(sl);
		m.Rolate(sl, PI / 2, 3);
		SL.push_back(sl);
		m.Rolate(sl, PI / 2, 3);
		SL.push_back(sl);
		m.Rolate(sl, PI / 2, 3);
		SL.push_back(sl);

		return SL;

	}


};

//������״��Ŀǰ�У������ı��Ρ�Բ��
class RandomModel {
private:
	//�����ı��ε��ĸ�����
	Vec4 v1, v2, v3, v4,v5,v6;

	//Բ�İ뾶
	double r;
	double w = cos(PI / 4);

	//�ı���+��Բ
	double a, b;//����С��Ϊ��
public:
	RandomModel() {};

	//˳ʱ�����������
	RandomModel(Vec4 v1, Vec4 v2, Vec4 v3, Vec4 v4) :v1(v1), v2(v2), v3(v3), v4(v4){};

	//Բ r:�뾶
	RandomModel(double r) :r(r) {};

	//���Ϳ�
	RandomModel(double a, double b) :a(a), b(b) {};

	//�ı���
	SplineSurface getSurface() {
		varray<Spline> SL1;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		SL1.resize(4);
		for (int i = 0; i < 4; i++) {
			SL1[i].m_Degree = 2;
			SL1[i].m_Knots = knots;
		}
		SL1[0].m_CtrlPts.push_back(v1);
		SL1[0].m_CtrlPts.push_back((v1 + v4) / 2);
		SL1[0].m_CtrlPts.push_back(v4);
		SL1[1].m_CtrlPts.push_back(v1);
		SL1[1].m_CtrlPts.push_back((v1 + v2) / 2);
		SL1[1].m_CtrlPts.push_back(v2);
		SL1[2].m_CtrlPts.push_back(v2);
		SL1[2].m_CtrlPts.push_back((v2 + v3) / 2);
		SL1[2].m_CtrlPts.push_back(v3);
		SL1[3].m_CtrlPts.push_back(v4);
		SL1[3].m_CtrlPts.push_back((v4 + v3) / 2);
		SL1[3].m_CtrlPts.push_back(v3);
		SplineSurface ss1;
		ss1.CoonsInterpolate(SL1);
		
		return ss1;
	}

	//Բ
	varray<SplineSurface> getArcSurface() {
		varray<SplineSurface> SS;
		varray<Spline> SL1,SL2;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		SL1.resize(4);
		SL2.resize(4);
		for (int i = 0; i < 4; i++) {
			SL1[i].m_Degree = 2;
			SL1[i].m_Knots = knots;
			SL2[i].m_Degree = 2;
			SL2[i].m_Knots = knots;
		}

		Vec4 v1 = { -r,0,0,1 };
		Vec4 v2 = { 0,r,0,1 };
		Vec4 v3 = { r,0,0,1 };
		Vec4 v4 = { 0,-r,0,1 };

		Vec4 p1 = { v1.x,v4.y,0,w };
		Vec4 p2 = { v1.x,v2.y,0,w };
		Vec4 p3 = { v3.x,v2.y,0,w };
		Vec4 p4 = { v3.x,v4.y,0,w };

		Vec4 v5 = { -r/2,0,0,1 };
		Vec4 v6 = { 0,r/2,0,1 };
		Vec4 v7 = { r/2,0,0,1 };
		Vec4 v8 = { 0,-r/2,0,1 };

		SL1[0].m_CtrlPts.push_back(v1);
		SL1[0].m_CtrlPts.push_back((v1 + v5) / 2);
		SL1[0].m_CtrlPts.push_back(v5);
		SL1[1].m_CtrlPts.push_back(v1);
		SL1[1].m_CtrlPts.push_back(p2);
		SL1[1].m_CtrlPts.push_back(v2);
		SL1[2].m_CtrlPts.push_back(v2);
		SL1[2].m_CtrlPts.push_back((v2 + v6) / 2);
		SL1[2].m_CtrlPts.push_back(v6);
		SL1[3].m_CtrlPts.push_back(v5);
		SL1[3].m_CtrlPts.push_back((v6 + v5) / 2);
		SL1[3].m_CtrlPts.push_back(v6);
		SplineSurface ss1,temp;
		ss1.CoonsInterpolate(SL1);
		SS.push_back(ss1);
		temp = ss1;
		Model_Solution m;
		m.Rolate(temp, PI / 2, 3);
		SS.push_back(temp);
		m.Rolate(temp, PI / 2, 3);
		SS.push_back(temp);
		m.Rolate(temp, PI / 2, 3);
		SS.push_back(temp);

		SL2[0].m_CtrlPts.push_back(v5);
		SL2[0].m_CtrlPts.push_back((v8 + v5) / 2);
		SL2[0].m_CtrlPts.push_back(v8);
		SL2[1].m_CtrlPts.push_back(v5);
		SL2[1].m_CtrlPts.push_back((v6 + v5) / 2);
		SL2[1].m_CtrlPts.push_back(v6);
		SL2[2].m_CtrlPts.push_back(v6);
		SL2[2].m_CtrlPts.push_back((v7 + v6) / 2);
		SL2[2].m_CtrlPts.push_back(v7);
		SL2[3].m_CtrlPts.push_back(v8);
		SL2[3].m_CtrlPts.push_back((v8 + v7) / 2);
		SL2[3].m_CtrlPts.push_back(v7);
		SplineSurface ss2;
		ss2.CoonsInterpolate(SL2);
		SS.push_back(ss2);



		return SS;
	}

	//�ı���+�ڿ�Բ
	varray<SplineSurface> getRecArcSurface();

};

//�ʷֺ��� outerΪ��������innerΪ��������addlinesΪ�����ߣ� genusΪÿ���������Ƿ�Ϊ�ף� allSurfΪ�ʷֽ��
//void quad(varray<Spline>& outer, varray<varray<Spline>>& inner, varray<Spline>& addlines,
//			varray<bool>& genus, varray<SplineSurface>& allSurf) {
//
//	varray<Spline>allLines;
//	//�ֱ��ȡ�ⲿ�������ߣ��ڲ��������ߺͶ��⸨����
//	//���и����ߵ������ǽ��׺ͱ߽���������
//	varray<SplineSurface> NS;
//	varray<varray<Spline>> surf;					//��ÿ�������ֱ�װ������
//	surf.push_back(outer);
//	for (auto& suf : inner)
//		surf.push_back(suf);
//	varray<varray<int>>seg;							//��ÿ�������Ŀɷָ��Լ�¼������0Ϊ���ɷ֣�1Ϊ�ɷָ�
//	seg.push_back(varray<int>(outer.size(), 1));    //Ĭ������Ϊ���ɷָ�
//	for (int i = 0; i < inner.size(); i++)
//		seg.push_back(varray<int>(inner[i].size(), 1));
//
//	SfCtainTreeNode* root = CreateSurfContainTree(surf, addlines, seg, genus);//���������������
//	 
//	QuadWithContainTree(root);						//�ʷ�
//
//	//���²����ǴӰ�������ȡ�����棬���Ե�����װ�ɺ���
//	queue<SfCtainTreeNode*> q;
//	q.push(root);
//	allLines.clear();
//	while (!q.empty())
//	{
//		varray<SplineSurface> tmpsf;
//		SfCtainTreeNode*cur = q.front();
//		q.pop();
//		//���뵱ǰ�ڵ���������
//		for (auto i : cur->quadPolNumber)
//		{
//			for (auto j : i) {
//				allLines.push_back(cur->allLines[j]);
//			}
//		}
//		//���뵱ǰ�ڵ���������
//		cur->GetSurfs(tmpsf);
//		for (auto& s : tmpsf) {
//			allSurf.push_back(s);
//		}
//
//		//�Ľ���Ȩ��
//		for (int i = 0; i < allSurf.size(); i++) {
//			for (int j = 0; j < allSurf[i].m_CtrlPts.size(); j++) {
//				if (allSurf[i].m_CtrlPts[j].w > 1) {
//					allSurf[i].m_CtrlPts[j].w = 1;
//				}	
//			}
//		}
//
//		list<SfCtainTreeNode*>::iterator it = cur->childs.begin();
//		for (; it != cur->childs.end(); it++) {
//			q.push(*it);
//		}
//	}
//}




//�������ܺ�����
class PublicSolution {

public:
	/**
	 *	ʩ���غɺ����ӿ�
	 *	str1�����ɵ���ģ���ļ�·������.txt��
	 *	str2: ����ķ����ļ�·��������.txt��
	 */
	void setWCandWF(const string str1, const string str2)
	{
		RWGeometric rwg;
		Model_Solution M;
		varray<SplineVolume> NVs;
		varray<varray<SplineSurface>> NFs;
		varray<SplineSurface> NF;
		TestBolcks::pList plist;

		rwg.ReadSplineVolume(str1, NVs);
		NFs = M.GetSurfaces(NVs);
		for (auto& SF : NFs) {
			for (auto& i : SF) {
				NF.push_back(i);
			}
		}
		rwg.WriteSplineSurface("E:\\kuang_models\\WcWfFile\\newSurface.txt", NF);
		string str,modelName;
		cout << "ɢ���ļ��Ѿ����ɣ��뵽E:\\kuang_models\\WcWfFile\\newSurface.txt�鿴......" << endl;
		cout << "������ģ�����ƣ�";
		cin >> modelName;
		cout << "�Ƿ��Ѽ�¼ʩ��Լ��������������Ϣ��";

		varray<int> wc;
		varray<int> wf;
		wc.clear();
		wf.clear();
		while (1) {
			cin >> str;
			if (str == "y" || str == "Y") {
				ifstream wcWfFile;
				string temp = "";
				int flag = 0;
				wcWfFile.open("E:\\kuang_models\\WcWfFile\\WcWf.txt", ios::binary);
				if (!wcWfFile) {
					cout << "�ļ���ʧ��!" << endl;
				}
				while (!wcWfFile.eof()) {
					wcWfFile >> temp;
					if (temp == "WC") {
						flag = 1;
						continue;
					}
					else if (temp == "WF") {
						flag = 2;
						continue;
					}
					switch (flag) {
					case 1: {
						int c = atoi(temp.c_str());
						wc.push_back(c);
						break;
					}
					case 2: {
						int f = atoi(temp.c_str());
						wf.push_back(f);
						break;
					}
					default:
						break;
					}
				}
				cout << "Լ���������" << wc.size() << endl;
				cout << "�غ��������" << wf.size() << endl;
				wcWfFile.close();

				//��¼ģ��Լ�����غɵ�����
				ofstream wcWfLog;
				wcWfLog.open("E:\\kuang_models\\WcWfFile\\WcWfLog.txt", ios::binary | ios::app);
				wcWfLog << modelName << endl;
				wcWfLog << "WC ";
				for (auto& i : wc) {
					wcWfLog << i << " ";
				}
				wcWfLog << endl;
				wcWfLog << "WF ";
				for (auto& i : wf) {
					wcWfLog << i << " ";
				}
				wcWfLog << endl;
				
				plist.OutputParaVolumeDataTxt(NVs, str2);

				//ʩ��Լ������
				varray<SplineVolume> NVS;
				rwg.ReadSplineVolume(str1, NVS);
				TestBolcks::pList plt;
				varray<varray<SplineSurface>> NS = M.GetSurfaces(NVS);
				varray<SplineSurface> NSf;
				int ssNum = 0;
				for (int i = 0; i < NS.size(); i++)
				{
					for (int j = 0; j < NS[i].size(); j++)
					{
						NSf.push_back(NS[i][j]);
						ssNum++;
					}
				}
				cout << "Ƭ����" << ssNum << endl;
				rwg.WriteSplineSurface(str2 + "������.txt", NSf);

				ofstream wcWf;
				wcWf.open(str2+ "ctrlptsIdx.txt", ios::binary|ios::app);

				varray<SplineSurface> WC;
				for (int i = 0; i < wc.size(); i++) {
					WC.push_back(NSf[wc[i]]);
				}
				varray<varray<int>>WCidx = plt.getfaceidx(NVS, WC);
				//plt.showdata();
				int cnum = 0;
				cout << "WC" << " ";
				wcWf << "WC" << " ";
				for (int i = 0; i < WCidx.size(); i++)
				{
					for (int j = 0; j < WCidx[i].size(); j++)
					{
						cnum++;
						cout << WCidx[i][j] << " ";
						wcWf << WCidx[i][j] << " ";
					}
				}
				cout << endl;
				wcWf << endl;
				varray<SplineSurface> WF;
				for (int i = 0; i < wf.size(); i++) {
					WF.push_back(NSf[wf[i]]);
				}
				varray<varray<int>>WFidx = plt.getfaceidx(NVS, WF);
				int fnum = 0;
				cout << "WF" << " ";
				wcWf << "WF" << " ";
				for (int i = 0; i < WFidx.size(); i++)
				{
					for (int j = 0; j < WFidx[i].size(); j++)
					{
						fnum++;
						cout << WFidx[i][j] << " ";
						wcWf << WFidx[i][j] << " ";
					}
				}
				wcWf.close();
				cout << endl;
				cout << "Լ�����Ƶ�����" << cnum << endl;
				cout << "�غɿ��Ƶ�����" << fnum << endl;

				wcWfLog << "Լ���������" << wc.size() << endl;
				wcWfLog << "�غ��������" << wf.size() << endl;
				wcWfLog << "Լ�����Ƶ�����" << cnum << endl;
				wcWfLog << "�غɿ��Ƶ�����" << fnum << endl;
				wcWfLog.close();
				

				break;
			}
			else
				cout << "�������룡���������룺";
		}

	}

	//�������
	SplineVolume loft(SplineSurface ss1, SplineSurface ss2) {
		varray<Vec4> vec;
		for (int i = 0; i < ss1.m_CtrlPts.size(); i++) {
			vec.push_back((ss1.m_CtrlPts[i] + ss2.m_CtrlPts[i]) / 2);
			vec[i].w = 1;
		}
		SplineVolume SV;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		SV.m_uNum = 3;//���Ƶ����
		SV.m_vNum = 3;
		SV.m_wNum = 3;
		SV.m_uKnots = knots;
		SV.m_uDegree = 2;
		SV.m_vKnots = knots;
		SV.m_vDegree = 2;
		SV.m_wKnots = knots;
		SV.m_wDegree = 2;

		for (int i = 0; i < ss1.m_CtrlPts.size(); i++) {
			SV.m_CtrlPts.push_back(ss1.m_CtrlPts[i]);
		}
		for (int i = 0; i < vec.size(); i++) {
			SV.m_CtrlPts.push_back(vec[i]);
		}
		for (int i = 0; i < ss2.m_CtrlPts.size(); i++) {
			SV.m_CtrlPts.push_back(ss2.m_CtrlPts[i]);
		}

		return SV;

	}

	//�������
	varray<SplineVolume> loft(varray<SplineSurface> ss1, varray<SplineSurface> ss2) {
		varray<SplineVolume> SV;
		for (int i = 0; i < ss1.size(); i++) {
			SplineVolume sv;
			sv = loft(ss1[i], ss2[i]);
			SV.push_back(sv);
		}
		return SV;
	}

	//��������ת
	void reverse(varray<Vec4> v1) {
		stack<Vec4> st;
		for (int i = 0; i < v1.size(); i++) {
			st.push(v1[i]);
		}
		for (int i = 0; i < st.size(); i++) {
			v1.push_back(st.top());
			st.pop();
		}
		
	}

	/**
	*	����Զ��ʷ֣� Spline����(xyƽ��) 
	*	axis:������x,y:(1,2)
	*	�Ƿ����ԭͼ�� 1���� 2����
	*	
	*/
	//�㾵��
	varray<Vec4> mirror(varray<Vec4> v1,int axis) {
		varray<Vec4> V;
		
		if (axis == 1) {
			for (auto& v : v1) {
				v.y = -v.y;
			}
		}
		if (axis == 2) {
			for (auto& v : v1) {
				v.x = -v.x;
			}
		}
		for (auto& v : v1) {
			V.push_back(v);
		}
			
		return V;
	}

	varray<Spline> mirror(Spline sl,int axis,int choice) {
		varray<Spline> SL;
		if (choice == 1) {
			SL.push_back(sl);
			
		}
		if (axis == 1) {
			for (int i = 0; i < sl.m_CtrlPts.size(); i++) {
				if (sl.m_CtrlPts[i].y != 0) {
					sl.m_CtrlPts[i].y = -sl.m_CtrlPts[i].y;
				}

			}
			reverse(sl.m_CtrlPts);
			SL.push_back(sl);
		}
		if (axis == 2) {
			for (int i = 0; i < sl.m_CtrlPts.size(); i++) {
				if (sl.m_CtrlPts[i].x != 0) {
					sl.m_CtrlPts[i].x = -sl.m_CtrlPts[i].x;
				}
			}
			reverse(sl.m_CtrlPts);
			SL.push_back(sl);
		}
		
		
		
		return SL;
	}
	
	varray<Spline> mirror(varray<Spline> sl, int axis,int choice){
		varray<Spline> SL;
		if (axis == 1) {
			for (int i = 0; i < sl.size(); i++) {
				varray<Spline> SL1;
				SL1 = mirror(sl[i], axis, choice);
				for (int j = 0; j < SL1.size(); j++) {
					SL.push_back(SL1[j]);
				}
			}
		}
		if (axis == 2) {
			for (int i = 0; i < sl.size(); i++) {
				varray<Spline> SL1;
				SL1 = mirror(sl[i], axis, choice);
				for (int j = 0; j < SL1.size(); j++) {
					SL.push_back(SL1[j]);
				}
			}
		}
		
		return SL;
	}

	varray<SplineSurface> mirror(SplineSurface& ss, int axis) {
		varray<SplineSurface> SS;
		SS.push_back(ss);
		ss.m_CtrlPts = mirror(ss.m_CtrlPts, axis);
		SS.push_back(ss);

		return SS;
	}

	varray<SplineSurface> mirror(varray<SplineSurface>& ss, int axis) {
		varray<SplineSurface> SS,SS1;
		for (auto& s : ss) {
			SS1 = mirror(s, axis);
			for (auto& s1 : SS1) {
				SS.push_back(s1);
			}
		}
		return SS;
	}

	varray<SplineVolume> mirror(SplineVolume& sv, int axis) {
		varray<SplineVolume> SV;
		SV.push_back(sv);
		sv.m_CtrlPts = mirror(sv.m_CtrlPts, axis);
		SV.push_back(sv);

		return SV;
	}

	varray<SplineVolume> mirror(varray<SplineVolume>& sv, int axis) {
		varray<SplineVolume> SV, SV1;
		for (auto& s : sv) {
			SV1 = mirror(s, axis);
			for (auto& s1 : SV1) {
				SV.push_back(s1);
			}
		}
		return SV;
	}

	//���VTK��ʽ
	void outPutVTK(varray<SplineVolume>& SV,string path) {
		varray<NurbsVol> NV;
		NV = NurbsTrans::SplinevolsToCvols(SV);
		CPolyParaVolume cp;       //���vtk�ļ��������
		cp = NV;
		cp.OutputParaVolumeDataVTK(path);
	}

	//����������ļ�,���üӺ�׺
	void outPutTXT(varray<SplineSurface>& ss,string path) {
		TestBolcks::pList plist;
		plist.OutputParaSurfaceDataTxt(ss, path);
	}

	void outPutTXT(varray<SplineVolume>& sv, string path) {
		TestBolcks::pList plist;
		plist.OutputParaVolumeDataTxt(sv, path);
	}
	
	//�ʷֳ����������
	void quadAdjust(varray<SplineSurface>& ss) {
		varray<SplineSurface> SS;

		//�Ľ�1��Ȩ�ش���1�Ķ�����Ϊ1
		for (int i = 0; i < ss.size(); i++) {
			for (int j = 0; j < ss[i].m_CtrlPts.size(); j++) {
				if (ss[i].m_CtrlPts[j].w > 1) {
					ss[i].m_CtrlPts[j].w = 1;
				}
			}
		}

		varray<Vec4> v;//���Ȩ��Ϊcos(PI/4)�ĵ�
		for (int i = 0; i < ss.size(); i++) {
			for (int j = 0; j < ss[i].m_CtrlPts.size(); j++) {
				if (abs(ss[i].m_CtrlPts[j].w - cos(PI / 4))<1e-4) {
					v.push_back(ss[i].m_CtrlPts[j]);
				}
			}
		}
		SS = ss;
		//�Ľ�2����v�е��غϵĵ�Ȩ������Ϊcos(PI/4)
		for (int i = 0; i < SS.size(); i++) {
			for (int j = 0; j < SS[i].m_CtrlPts.size(); j++) {
				for (auto&v1 : v) {
					if (abs(SS[i].m_CtrlPts[j].x - v1.x) < 1e-4&&
						abs(SS[i].m_CtrlPts[j].y - v1.y) < 1e-4&&
						abs(SS[i].m_CtrlPts[j].z - v1.z) < 1e-4) {
						SS[i].m_CtrlPts[j].w = cos(PI / 4);
					}
				}
			}
		}

		ss.clear();
		ss = SS;
	}


	//�����м���Ȩ��
	void calWeight(Spline &sl) {
		
		Vec4 v1 = sl.m_CtrlPts[0];
		Vec4 v2 = sl.m_CtrlPts[1];
		Vec4 v3 = sl.m_CtrlPts[2];

		
		Vec4 p1 = v2 - v1;
		Vec4 p2 = v3 - v1;
		double dotP = p1.Dot(p2);
		double MA = p1.Magnitude();
		double MB = p2.Magnitude();
		double w = dotP / (MA*MB);
		
		sl.m_CtrlPts[0].w = 1;
		sl.m_CtrlPts[1].w = w;
		sl.m_CtrlPts[2].w = 1;
	}
	static double distance(Vec4 p1, Vec4 p2)
	{
		return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
	}

	//Ĭ�ϱ߽綼���ڵ�һ���� �ҳ�����ԭ������ĵ�
	Vec4 calMinPoint(varray<Spline> &sl) {
		double dist,temp;
		Vec4 vecDistMin;//����ԭ������ĵ�
		Vec4 P0 = { 0,-50,0,0 };//�ο���
		temp = distance(sl[0].m_CtrlPts[0], P0);
		vecDistMin = sl[0].m_CtrlPts[0];
		for (auto&i : sl) {
			int cptLenth = i.m_CtrlPts.size();//���Ƶ����
			dist = distance(i.m_CtrlPts[0], P0);
			if (dist < temp) {
				temp = dist;
				vecDistMin = i.m_CtrlPts[0];
			}
			dist = distance(i.m_CtrlPts[cptLenth-1], P0);
			if (dist < temp) {
				temp = dist;
				vecDistMin = i.m_CtrlPts[cptLenth-1];
			}
		}
		cout << vecDistMin.x << vecDistMin.y << vecDistMin.z << endl;
		return vecDistMin;

	}

	//�Ա߽��������
	void sortEdg(varray<Spline> &sl) {
		Vec4 vMin = calMinPoint(sl);
		varray<Spline> SL;//����ź���ı߽���
		varray<Spline> SL1;
		Spline temp,Line4;
		//�ҵ��������С���йص������߽���
		for (auto&i : sl) {
			if (JudgeTwoPointsCoincide(vMin, i.m_CtrlPts[0])) {
				SL.push_back(i);
			}
			else if (JudgeTwoPointsCoincide(vMin, i.m_CtrlPts[2])) {
				Vec4 temp;
				temp = i.m_CtrlPts[0];
				i.m_CtrlPts[0] = i.m_CtrlPts[2];
				i.m_CtrlPts[2] = temp;
				SL.push_back(i);
			}
			else {
				SL1.push_back(i);
			}
		}
		//ȷ���߽�˳��
		Vec4 v1 = SL[0].m_CtrlPts[2] - SL[0].m_CtrlPts[0];
		Vec4 v2 = SL[1].m_CtrlPts[2] - SL[1].m_CtrlPts[0];
		v1.Normalize();
		v2.Normalize();
		
		Vec4 cross = v1.Cross(v2);//���������
		if (cross.z < 0)
		{
			temp = SL[0];
			SL[0] = SL[1];
			SL[1] = temp;
		}

		for (auto&i : SL1) {
			if (JudgeTwoPointsCoincide(SL[1].m_CtrlPts[2], i.m_CtrlPts[0])) {
				SL.push_back(i);
			}
			else if (JudgeTwoPointsCoincide(SL[1].m_CtrlPts[2], i.m_CtrlPts[2])) {
				Vec4 temp;
				temp = i.m_CtrlPts[0];
				i.m_CtrlPts[0] = i.m_CtrlPts[2];
				i.m_CtrlPts[2] = temp;
				SL.push_back(i);
			}
			else if (JudgeTwoPointsCoincide(SL[0].m_CtrlPts[2], i.m_CtrlPts[0])) {
				Line4 = i;
			}
			else if (JudgeTwoPointsCoincide(SL[0].m_CtrlPts[2], i.m_CtrlPts[2])) {
				Vec4 temp;
				temp = i.m_CtrlPts[0];
				i.m_CtrlPts[0] = i.m_CtrlPts[2];
				i.m_CtrlPts[1] = i.m_CtrlPts[1];
				i.m_CtrlPts[2] = temp;
				Line4 = i;
			}
			
		}
		SL.push_back(Line4);
		//����Ȩ��
		for (auto&i : SL) {
			calWeight(i);
		}
		sl = SL;

		
		

		/*cout << SL[0].m_CtrlPts[0].x << SL[0].m_CtrlPts[0].y << SL[0].m_CtrlPts[0].z << endl;
		cout << SL[0].m_CtrlPts[2].x << SL[0].m_CtrlPts[2].y << SL[0].m_CtrlPts[2].z << endl;
		cout << SL[1].m_CtrlPts[0].x << SL[1].m_CtrlPts[0].y << SL[1].m_CtrlPts[0].z << endl;
		cout << SL[1].m_CtrlPts[2].x << SL[1].m_CtrlPts[2].y << SL[1].m_CtrlPts[2].z << endl;
		cout << SL[2].m_CtrlPts[0].x << SL[2].m_CtrlPts[0].y << SL[2].m_CtrlPts[0].z << endl;
		cout << SL[2].m_CtrlPts[2].x << SL[2].m_CtrlPts[2].y << SL[2].m_CtrlPts[2].z << endl;
		cout << SL[3].m_CtrlPts[0].x << SL[3].m_CtrlPts[0].y << SL[3].m_CtrlPts[0].z << endl;
		cout << SL[3].m_CtrlPts[2].x << SL[3].m_CtrlPts[2].y << SL[3].m_CtrlPts[2].z << endl;*/
	}

	void quadAdjustagain(varray<SplineSurface>& ss) {
		varray<SplineSurface> SS;
		for (auto& i : ss) {
			SplineSurface ss;
			varray<Spline> sl;
			i.GetEdgeLines(sl);
			sortEdg(sl);
			ss.CoonsInterpolate(sl);
			SS.push_back(ss);
		}
		ss.clear();
		for (auto&i : SS) {
			ss.push_back(i);
		}

		
	}
	//�ʷֺ��� outerΪ��������innerΪ��������addlinesΪ�����ߣ� genusΪÿ���������Ƿ�Ϊ�ף� allSurfΪ�ʷֽ��
	void quad(varray<Spline>& outer, varray<varray<Spline>>& inner, varray<Spline>& addlines,
		varray<bool>& genus, varray<SplineSurface>& allSurf) {

		varray<Spline>allLines;
		//�ֱ��ȡ�ⲿ�������ߣ��ڲ��������ߺͶ��⸨����
		//���и����ߵ������ǽ��׺ͱ߽���������
		varray<SplineSurface> NS;
		varray<varray<Spline>> surf;					//��ÿ�������ֱ�װ������
		surf.push_back(outer);
		for (auto& suf : inner)
			surf.push_back(suf);
		varray<varray<int>>seg;							//��ÿ�������Ŀɷָ��Լ�¼������0Ϊ���ɷ֣�1Ϊ�ɷָ�
		seg.push_back(varray<int>(outer.size(), 1));    //Ĭ������Ϊ���ɷָ�
		for (int i = 0; i < inner.size(); i++)
			seg.push_back(varray<int>(inner[i].size(), 1));

		SfCtainTreeNode* root = CreateSurfContainTree(surf, addlines, seg, genus);//���������������

		QuadWithContainTree(root);						//�ʷ�

		//���²����ǴӰ�������ȡ�����棬���Ե�����װ�ɺ���
		queue<SfCtainTreeNode*> q;
		q.push(root);
		allLines.clear();
		while (!q.empty())
		{
			varray<SplineSurface> tmpsf;
			SfCtainTreeNode*cur = q.front();
			q.pop();
			//���뵱ǰ�ڵ���������
			for (auto i : cur->quadPolNumber)
			{
				for (auto j : i) {
					allLines.push_back(cur->allLines[j]);
				}
			}
			//���뵱ǰ�ڵ���������
			cur->GetSurfs(tmpsf);
			for (auto& s : tmpsf) {
				allSurf.push_back(s);
			}

			
			
			list<SfCtainTreeNode*>::iterator it = cur->childs.begin();
			for (; it != cur->childs.end(); it++) {
				q.push(*it);
			}
		}
		//�ʷֵ���
		//quadAdjust(allSurf);
		quadAdjustagain(allSurf);
	}

	

	//��ģ��ȥ��
	static void removeRepeatVols(varray<SplineVolume>& Vols) {
		varray<varray<Vec4>> Rmaps;
		int len = Vols.size();
		int j = 0;
		for (int i = 0; i < len; i++) {
			varray<Vec4> cpts;
			cpts.push_back(Vols[i].GetVolPoint(0, 0, 0));
			cpts.push_back(Vols[i].GetVolPoint(0, 1, 0));
			cpts.push_back(Vols[i].GetVolPoint(1, 0, 0));
			cpts.push_back(Vols[i].GetVolPoint(0, 1, 1));
			bool isFind = false;
			for (auto& var : Rmaps) {
				if (distance(cpts[0], var[0]) < 0.001
					&&distance(cpts[1], var[1]) < 0.001
					&&distance(cpts[2], var[2]) < 0.001
					&&distance(cpts[3], var[3]) < 0.001) {
					isFind = true;
					break;
				}
			}
			if (isFind) {
				continue;
			}
			else {
				Rmaps.push_back(cpts);
				Vols[j++] = Vols[i];
			}

		}
		Vols.resize(j);
	}

	/*
	 *	ת���ӿڣ��������д�����ݽṹ���
	 */
	void Transform(Vec4& v,point4d& p) {
		p.x = v.x;
		p.y = v.y;
		p.z = v.z;
		p.w = v.w;
	}

	void Transform(point4d& p, Vec4& v) {
		v.x = p.x;
		v.y = p.y;
		v.z = p.z;
		v.w = p.w;
	}

	void Transform(varray<Vec4>& v, vector<point4d>& p) {
		for (int i = 0; i < v.size(); i++) {
			Transform(v[i], p[i]);
		}
	}

	void Transform(vector<point4d>& p, varray<Vec4>& v) {
		v.resize(p.size());
		for (int i = 0; i < p.size(); i++) {
			Transform(p[i], v[i]);
		}
	}

	void Transform(SplineSurface &SS, YN::NurbsSurface &NS) {
		NS._u_Degree = SS.m_uDegree;
		NS._v_Degree = SS.m_vDegree;
		NS._u_Num = SS.m_uNum;
		NS._v_Num = SS.m_vNum;
		vector<point4d> p;
		Transform(SS.m_CtrlPts, p);
		auto p1 = make_shared<vector<point4d>>(p);
		NS._ControlPts = p1;
		
		for (auto& i : SS.m_uKnots) {
			NS._u_Knots->push_back(i);
		}
		for (auto& i : SS.m_vKnots) {
			NS._v_Knots->push_back(i);
		}
	
	}

	void Transform(YN::NurbsSurface &NS, SplineSurface &SS) {
		SS.m_uDegree = NS._u_Degree;
		SS.m_vDegree = NS._v_Degree;
		SS.m_uNum = NS._u_Num;
		SS.m_vNum = NS._v_Num;
		vector<double> v = *NS._u_Knots;//����ָ�룬ָ����
		for (auto&i : v) {
			SS.m_uKnots.push_back(i);
		}
		
		v = *NS._v_Knots;
		for (auto&i : v) {
			SS.m_vKnots.push_back(i);
		}
		
		Transform(*NS._ControlPts, SS.m_CtrlPts);
	}

	void Transform(varray<SplineSurface> &SS, vector<YN::NurbsSurface> &NS) {
		for (int i = 0; i < SS.size(); i++) {
			Transform(SS[i], NS[i]);
		}
	}

	void Transform(vector<YN::NurbsSurface> &NS, varray<SplineSurface> &SS) {
		SS.resize(NS.size());
		for (int i = 0; i < NS.size(); i++) {
			Transform(NS[i], SS[i]);
		}
	}

	void Transform(SplineVolume &SV, YN::NurbsVol &NV) {
		NV._u_Degree = SV.m_uDegree;
		NV._v_Degree = SV.m_vDegree;
		NV._w_Degree = SV.m_wDegree;
		NV._u_Num = SV.m_uNum;
		NV._v_Num = SV.m_vNum;
		NV._w_Num = SV.m_wNum;
		for (auto& i : SV.m_uKnots) {
			NV._u_Knots->push_back(i);
		}
		for (auto& i : SV.m_vKnots) {
			NV._v_Knots->push_back(i);
		}
		for (auto& i : SV.m_wKnots) {
			NV._w_Knots->push_back(i);
		}
		vector<point4d> p;
		Transform(SV.m_CtrlPts, p);
		for (auto& i : p) {
			NV._ControlPts->push_back(i);
		}
	}

	void Transform(varray<SplineVolume> &SV, vector<YN::NurbsVol> &NV) {
		for (int i = 0; i < SV.size(); i++) {
			Transform(SV[i], NV[i]);
		}
	}
};
	




//--------------------ģ��--------------------//

//���������
class CarPart {

public:
	
	//����ߴ����
	double r = 1.5;//�ն�����
	double angle = 102.53*PI / 180;//б�߽Ƕ�
	double angle1 = (90 - 75.07)*PI / 180;//������ƫ�ƽǶ�1
	double angle2 = (90 - 71.07)*PI / 180;//������ƫ�ƽǶ�2,��������ĸ�Ϊ�����ĽǶ�
	double angle3 = (90 - 70.07)*PI / 180;//������ƫ�ƽǶ�3
	double angle4 = (90 - 70.07)*PI / 180;//������ƫ�ƽǶ�4
	double angle5 = (90 - 69.07)*PI / 180;//������ƫ�ƽǶ�5
	double angle6 = (90 - 73.07)*PI / 180;//������ƫ�ƽǶ�6
	double Y_L = 58.22*cos(angle1);//��Сƽ��֮�����
	double Y_L1 = 24.21;//�ն���������
	double Y_L2;//����б��
	double Y_L3 = 40;//�ױ߳���
	double Y_L4=2;//�����������߳���
	double Y_L5 = 15;//���泤��
	double Y_H = 5;//�������Ӳ��ֿ��
	double Y_H1 = 10, Y_H3 = 13;//�����ն��߶�
	double Y_H2;

	//���������Ӳ���
	double Y_L7 = 17.6;//����ױ��м䲿�ֳ�
	double Y_L8 = 37.5;//���岿������ͼ�߳�
	double Y_L9 = 30;//���岿������ͼ��
	
	//���Ӵ��滻����
	double Y_L6 = (Y_L3 / 2 - cos(PI - angle)*Y_L)*2;//��Բ�����Ӳ�����Χ���
	 
	//���ƽ��ƽ�ƾ���
	double s1 = 5-5*sqrt(2)/4;
	double s2 = 5 + 5 * sqrt(2) / 4;
	double s3 = 20;
	double s4 = 30;
	//����Բ���뾶
	double r1 = ((s3 - s2) > Y_L6) ? Y_L6/4 : (s3 - s2)/4, r2 = Y_L4 / 4;
	//Բ���������ȡ�������ȡ�Բ������
	double h1 = sin(angle3)*Y_H, h2 = 10;

	Model_Solution m;
	PublicSolution ps;
	RWGeometric rwg;

	void initialize() {
		Y_H2 = Y_L2 - 3 * Y_H - Y_H1 - Y_H3 - Y_H / 2;//�������ն��߶�
	}

	//��Բ�����ӵ��滻����
	varray<SplineVolume> part() {
		varray<SplineSurface> SS,SS1;
		//��
		RandomModel rm((s3-s2), Y_L4);
		SS = rm.getRecArcSurface();
		//��
		RandomModel rm1((s3 - s2),Y_L6 );
		SS1 = rm1.getRecArcSurface();

		m.Trans(SS1, (sin(PI - angle)*Y_L / cos(angle3)-(Y_L-Y_H))*cos(angle3), 3);
		m.Trans(SS1, (sin(PI - angle)*Y_L / cos(angle3) - (Y_L - Y_H))*sin(angle3), -2);
		varray<SplineVolume> SV;
		SV = ps.loft(SS, SS1);

		
		rwg.WriteSplineVolume("E:\\kuang_models\\newvolume.txt", SV);
		return SV;
	}

	//Բ�������Ӳ���
	varray<SplineVolume> cylinderLoft() {

		varray<SplineVolume> SV, temp;
		SplineVolume tmp;
		PublicSolution ps;

		//��
		RandomModel rm(r2);
		varray<SplineSurface> SS;
		SS = rm.getArcSurface();
		//��ת45��
		m.Rolate(SS, PI / 4, 3);
		//��
		RandomModel rm1(r1);
		varray<SplineSurface> SS1;
		SS1 = rm1.getArcSurface();
		m.Rolate(SS1, PI / 4, 3);

		m.Trans(SS1, (sin(PI - angle)*Y_L / cos(angle3) - (Y_L - Y_H))*cos(angle3), 3);
		m.Trans(SS1, (sin(PI - angle)*Y_L / cos(angle3) - (Y_L - Y_H))*sin(angle3), -2);
		SV = ps.loft(SS, SS1);

		temp = m.CreatSweepVol(SS1, h2, 3);
		
		for (int i = 0; i < temp.size(); i++) {
			SV.push_back(temp[i]);
		}
		
		temp = part();
		for (int i = 0; i < temp.size(); i++) {
			SV.push_back(temp[i]);
		}
		m.Rolate(SV, -PI / 2, 1);
		m.Trans(SV, s2 + (s3 - s2) / 2, -3);
		m.Trans(SV, cos(angle3)*(Y_L  - Y_H), 2);
		m.Trans(SV, sin(angle3)*(Y_L  - Y_H), 3);


		//rwg.WriteSplineSurface("E:\\kuang_models\\newsurface.txt", SS);
		rwg.WriteSplineVolume("E:\\kuang_models\\newVolume1.txt", SV);
		return SV;
	}

	//��---����������������
	void carPartLine(varray<Spline> &outLine, varray<varray<Spline>> &inLine, varray<Spline> &addLine) {
		initialize();
		Vec4 v1 = { Y_L1 / 2,Y_H / 2,0,1 };
		Vec4 v2 = { -Y_L1 / 2,Y_H / 2,0,1 };
		Vec4 v3 = { -v2.x - Y_H1 / tan(PI - angle),Y_H / 2 + Y_H1,0,1 };
		Vec4 v4 = { v2.x + Y_H1 / tan(PI - angle),Y_H / 2 + Y_H1,0,1 };


		Vec4 v6 = { v2.x + (Y_H1 + Y_H) / tan(PI - angle),v2.y + Y_H1 + Y_H,0,1 };
		Vec4 v5 = { -v6.x,v6.y,0,1 };
		Vec4 v7 = { v6.x + Y_H2 / tan(PI - angle),v6.y + Y_H2 ,0,1 };
		Vec4 v8 = { -v7.x,v7.y,0,1 };

		Vec4 v10 = { -Y_L4/2,Y_L-Y_H,0,1 };
		Vec4 v9 = { v10.x - Y_H3 / tan(PI - angle),v10.y - Y_H3,0,1 };
		Vec4 v11 = { -v10.x,v10.y ,0,1 };
		Vec4 v12 = { -v9.x,v9.y,0,1 };

		//����������
		Vec4 p1 = { -Y_L3 / 2,0,0,1 };
		Vec4 p2 = { -Y_L3 / 2 + cos(PI - angle)*Y_L,sin(PI - angle)*Y_L2,0,1 };
		Vec4 p3 = { -p2.x,p2.y,0,1 };
		Vec4 p4 = { -p1.x,p1.y,0,1 };

		Vec4 p5 = { -Y_L7 / 2,0,0,1 };
		Vec4 p6 = { Y_L7 / 2,0,0,1 };

		varray<Spline> SL, inner1, inner2, inner3;//�洢������������

		Spline sl, temp;
		Spline sl_1, temp_1;
		Spline sl_2, temp_2;

		//�·�Բ��
		Circle0 cc1(r, angle);
		sl = cc1.getCircle();
		m.Rolate(sl, PI - angle / 2, 3);
		sl_2 = sl;
		m.Trans(sl, r + Y_H / 2, 2);
		m.Trans(sl, v2.x, 1);
		inner1.push_back(sl);
		sl_1 = sl;
		m.Trans(sl_1, Y_H1 + Y_H, 2);
		m.Trans(sl_1, (Y_H1 + Y_H) / tan(PI - angle), 1);
		inner2.push_back(sl_1);
		
		m.Trans(sl_2, Y_L-Y_H-Y_H3+r, 2);
		m.Trans(sl_2, Y_L4/2+Y_H3/tan(PI-angle), -1);
		inner3.push_back(sl_2);

		Circle0 cc2(r, PI - angle);
		temp = cc2.getCircle();
		m.Rolate(temp, (PI - angle) / 2, 3);
		temp_2 = temp;
		m.Trans(temp, Y_H1 + Y_H / 2 - r, 2);
		m.Trans(temp, v4.x, 1);
		inner1.push_back(temp);
		temp_1 = temp;
		m.Trans(temp_1, Y_H2 + Y_H, 2);
		m.Trans(temp_1, (Y_H2 + Y_H) / tan(PI - angle), 1);
		inner2.push_back(temp_1);
		
		m.Trans(temp_2, Y_L-Y_H-r, 2);
		m.Trans(temp_2, Y_L4/2, -1);
		inner3.push_back(temp_2);

		Spline0 sl3(sl.m_CtrlPts[0], temp.m_CtrlPts[2]);
		sl = sl3.getSpline();
		inner1.push_back(sl);

		Spline0 sl4(sl_1.m_CtrlPts[0], temp_1.m_CtrlPts[2]);
		sl = sl4.getSpline();
		inner2.push_back(sl);

		Spline0 sl5(sl_2.m_CtrlPts[0], temp_2.m_CtrlPts[2]);
		sl = sl5.getSpline();
		inner3.push_back(sl);

		//����
		inner1 = ps.mirror(inner1, 2, 1);
		inner2 = ps.mirror(inner2, 2,1);
		inner3 = ps.mirror(inner3, 2,1);

		Spline0 sl1(v1, v2);
		sl = sl1.getSpline();
		inner1.push_back(sl);
		Spline0 sl2(v4, v3);
		sl = sl2.getSpline();
		inner1.push_back(sl);
		Spline0 sl11(v5, v6);
		sl = sl11.getSpline();
		inner2.push_back(sl);
		Spline0 sl12(v7, v8);
		sl = sl12.getSpline();
		inner2.push_back(sl);
		Spline0 sl13(v12, v9);
		sl = sl13.getSpline();
		inner3.push_back(sl);
		Spline0 sl14(v10, v11);
		sl = sl14.getSpline();
		inner3.push_back(sl);

		//��������
		inLine.push_back(inner1);
		inLine.push_back(inner2);
		inLine.push_back(inner3);

		varray<Spline> inSL;
		for (auto& i : inLine) {
			for (auto& j : i) {
				inSL.push_back(j);
			}
		}


		//������
		Spline0 nl1(p1, p2);
		sl = nl1.getSpline();
		SL.push_back(sl);

		/*Spline0 nl2(p2, v10);
		sl = nl2.getSpline();
		SL.push_back(sl);

		Spline0 nl3(v11, p3);
		sl = nl3.getSpline();
		SL.push_back(sl);*/

		Spline0 nl3(p2, p3);
		sl = nl3.getSpline();
		SL.push_back(sl);


		Spline0 nl4(p3, p4);
		sl = nl4.getSpline();
		SL.push_back(sl);

		Spline0 nl5(p4, p6);
		sl = nl5.getSpline();
		SL.push_back(sl);

		Spline0 nl12(p6, p5);
		sl = nl12.getSpline();
		SL.push_back(sl);

		Spline0 nl13(p5, p1);
		sl = nl13.getSpline();
		SL.push_back(sl);

		outLine = SL;
		SL.clear();

		Spline0 nl6(v4, v6);
		sl = nl6.getSpline();
		SL.push_back(sl);

		Spline0 nl7(v7, v9);
		sl = nl7.getSpline();
		SL.push_back(sl);

		Spline0 nl8(p2, v10);
		sl = nl8.getSpline();
		SL.push_back(sl);

		Spline0 nl9(v11, p3);
		sl = nl9.getSpline();
		SL.push_back(sl);

		Spline0 nl10(v2, p5);
		sl = nl10.getSpline();
		SL.push_back(sl);

		Spline0 nl11(v1, p6);
		sl = nl11.getSpline();
		SL.push_back(sl);
		addLine = SL;

		for (auto&i : addLine) {
			inSL.push_back(i);
		}
		for (auto&i : outLine) {
			inSL.push_back(i);
		}
		
		//rwg.WriteSpline("E:\\kuang_models\\CarpartinSpline.txt", SL);
		rwg.WriteSpline("E:\\kuang_models\\CarpartoutSpline.txt", outLine);
		rwg.WriteSpline("E:\\kuang_models\\CarpartinSpline.txt", inSL);
		rwg.WriteSpline("E:\\kuang_models\\CarpartaddSpline.txt", addLine);
	}



	//�������
	void carPart() {
		Model_Solution m;
		varray<Spline> outLine, outLine1, outLine2, outLine3, outLine4, outLine5;
		varray<varray<Spline>> inLine, inLine1, inLine2, inLine3, inLine4,inLine5;
		varray<Spline> addLine, addLine1, addLine2, addLine3, addLine4, addLine5;
		
		Y_L2 = Y_L / cos(angle2);
		carPartLine(outLine1, inLine1, addLine1);
		Y_L2 = Y_L / cos(angle3);
		carPartLine(outLine2, inLine2, addLine2);
		Y_L2 = Y_L / cos(angle4);
		carPartLine(outLine3, inLine3, addLine3);
		Y_L2 = Y_L / cos(angle5);
		carPartLine(outLine4, inLine4, addLine4);
		Y_L2 = Y_L / cos(angle6);
		carPartLine(outLine5, inLine5, addLine5);
		varray<bool> genus, genus1, genus2, genus3, genus4, genus5;
		varray<SplineSurface> allSurf, allSurf1, allSurf2, allSurf3, allSurf4, allSurf5;//����ʷֽ��
		
		/*genus1.resize(4);
		genus1[0] = false;
		genus1[1] = true;
		genus1[2] = true;
		genus1[3] = true;*/
		genus2.resize(4);
		genus2[0] = false;
		genus2[1] = true;
		genus2[2] = true;
		genus2[3] = true;
		/*genus3.resize(4);
		genus3[0] = false;
		genus3[1] = true;
		genus3[2] = true;
		genus3[3] = true;*/
		/*genus4.resize(4);
		genus4[0] = false;
		genus4[1] = true;
		genus4[2] = true;
		genus4[3] = true;*/
		/*genus5.resize(4);
		genus5[0] = false;
		genus5[1] = true;
		genus5[2] = true;
		genus5[3] = true;*/
		
		//ת�Ƶ���һ����
		for (auto&i : outLine2) {
			m.Trans(i, 50, 1);
		}
		for (auto&i : inLine2) {
			m.Trans(i, 50, 1);
		}
		for (auto&i : addLine2) {
			m.Trans(i,50, 1);
		}
		
		//ps.quad(outLine1, inLine1, addLine1, genus1, allSurf1);
		ps.quad(outLine2, inLine2, addLine2, genus2, allSurf2);
		//ps.quad(outLine3, inLine3, addLine3, genus3, allSurf3);
		//ps.quad(outLine4, inLine4, addLine4, genus4, allSurf4);
		//ps.quad(outLine5, inLine5, addLine5, genus5, allSurf5);
		for (auto&i : allSurf2) {
			m.Trans(i, 50, -1);
		}
		//rwg.WriteSplineSurface("E:\\kuang_models\\Carpartsurf.txt", allSurf2);

		

		/*m.Rolate(allSurf1, angle2, 1);
		m.Trans(allSurf1, s1, -3);*/
		m.Rolate(allSurf2, angle3, 1);
		m.Trans(allSurf2, s2, -3);
		allSurf3 = allSurf2;
		allSurf1 = allSurf2;
		m.Trans(allSurf3, s3-s2, -3);
		m.Trans(allSurf1, s2-s1, 3);
		allSurf5 = allSurf1;
		m.Trans(allSurf5, s1, 3);
		allSurf4 = allSurf3;
		m.Trans(allSurf4, s4-s3, -3);
		
		varray<SplineSurface> SS;
		
		for (int j = 0; j < allSurf1.size(); j++) {
			allSurf.push_back(allSurf1[j]);
		}
		for (int j = 0; j < allSurf2.size(); j++) {
			allSurf.push_back(allSurf2[j]);
		}
		for (int j = 0; j < allSurf3.size(); j++) {
			allSurf.push_back(allSurf3[j]);
		}
		for (int j = 0; j < allSurf4.size(); j++) {
			allSurf.push_back(allSurf4[j]);
		}
		for (int j = 0; j < allSurf4.size(); j++) {
			allSurf.push_back(allSurf5[j]);
		}

		
		varray<Spline> SL;
		for (int i = 0; i < inLine.size(); i++) {
			for (int j = 0; j < inLine[i].size(); j++) {
				SL.push_back(inLine[i][j]);
			}
		}
		varray<SplineVolume> SV,SV1,SV2,SV3, SV4,SV5;
		SV1 = m.CreatSweepVol(allSurf5, s1, -3);
		SV2 = m.CreatSweepVol(allSurf1, s2-s1, -3);
		SV3 = m.CreatSweepVol(allSurf2, s3-s2, -3);
		SV4 = m.CreatSweepVol(allSurf3, s4-s3, -3);
		/*SV1 = ps.loft(allSurf5, allSurf1);
		SV2 = ps.loft(allSurf1, allSurf2);
		SV3 = ps.loft(allSurf2, allSurf3);
		SV4 = ps.loft(allSurf3, allSurf4);*/
		for (auto& i : SV1) {
			SV.push_back(i);
		}
		for (auto& i : SV2) {
			SV.push_back(i);
		}
		for (auto& i : SV3) {
			SV.push_back(i);
		}
		
		for (auto& i : SV4) {
			SV.push_back(i);
		}
		for (auto&i : SV) {
			i.OrderCtrlPts(i);
		}
		/*SV[2].OrderCtrlPts(SV[2]);
		SV[3].OrderCtrlPts(SV[3]);
		SV[4].OrderCtrlPts(SV[4]);
		SV[5].OrderCtrlPts(SV[5]);
		SV[16].OrderCtrlPts(SV[16]);
		SV[17].OrderCtrlPts(SV[17]);
		SV[19].OrderCtrlPts(SV[19]);
		SV[20].OrderCtrlPts(SV[20]);
		SV[21].OrderCtrlPts(SV[21]);

		SV[26].OrderCtrlPts(SV[26]);
		SV[27].OrderCtrlPts(SV[27]);
		SV[28].OrderCtrlPts(SV[28]);
		SV[29].OrderCtrlPts(SV[29]);
		SV[40].OrderCtrlPts(SV[40]);
		SV[41].OrderCtrlPts(SV[41]);
		SV[43].OrderCtrlPts(SV[43]);
		SV[44].OrderCtrlPts(SV[44]);
		SV[45].OrderCtrlPts(SV[45]);
		
		SV[50].OrderCtrlPts(SV[50]);
		SV[51].OrderCtrlPts(SV[51]);
		SV[52].OrderCtrlPts(SV[52]);
		SV[53].OrderCtrlPts(SV[53]);
		SV[64].OrderCtrlPts(SV[64]);
		SV[65].OrderCtrlPts(SV[65]);
		SV[67].OrderCtrlPts(SV[67]);
		SV[68].OrderCtrlPts(SV[68]);
		SV[69].OrderCtrlPts(SV[69]);
		
		SV[74].OrderCtrlPts(SV[74]);
		SV[75].OrderCtrlPts(SV[75]);
		SV[76].OrderCtrlPts(SV[76]);
		SV[77].OrderCtrlPts(SV[77]);
		SV[88].OrderCtrlPts(SV[88]);
		SV[89].OrderCtrlPts(SV[89]);
		SV[91].OrderCtrlPts(SV[91]);
		SV[92].OrderCtrlPts(SV[92]);
		SV[93].OrderCtrlPts(SV[93]);*/

		SV.erase(SV.begin()+49);
		SV5 = cylinderLoft();
		for (int i = 0; i < SV5.size(); i++) {
			SV.push_back(SV5[i]);
		}
		SV[106].OrderCtrlPts(SV[106]);
		SV[108].OrderCtrlPts(SV[108]);
		for (auto& i : SV) {
			i.OrderCtrlPts(i);
		}
		cout << SV.size() << endl;
		m.Trans(SV, Y_L8 / 2, 2);

		cout << "��������......" << endl;
		
		varray<SplineVolume> SVtemp;
		m.MirrorVols(SV, SVtemp, 2);
		for (auto& i : SVtemp) {
			i.OrderCtrlPts(i);
		}
		for (auto& i : SVtemp) {
			SV.push_back(i);
		}
		
		
		//rwg.WriteSplineSurface("E:\\kuang_models\\CarpartSurface.txt", allSurf);
		//rwg.WriteSplineVolume("E:\\kuang_models\\CarpartVolume.txt", SV);
		//rwg.WriteSpline("E:\\kuang_models\\CarpartoutSpline.txt", outLine);
		//rwg.WriteSpline("E:\\kuang_models\\CarpartinSpline.txt", SL);
		//rwg.WriteSpline("E:\\kuang_models\\CarpartaddSpline.txt", addLine);

		m.Rolate(SV, PI / 2, 3);
		m.Trans(SV, Y_L3 / 2, 2);
		m.Trans(SV, Y_L8 / 2, 1);
		m.Trans(SV, Y_L9, 3);

		SVtemp.clear();
		//��ȡ���岿��
		rwg.ReadSplineVolume("E:\\kuang_models\\kuang\\carpart\\AllcarModel.txt", SVtemp);
		for (auto& i : SVtemp) {
			SV.push_back(i);
		}
		
		/*YN::WingStruct ws("E:\\kuang_models\\kuang\\carpartctrlpts.txt", true);
		vector<YN::NurbsSurface> YNSS;
		varray<SplineSurface> YNSS1;
		YNSS = ws.getModelSurface();
		ps.Transform(YNSS, YNSS1);*/
		
		/*for (int i = 0; i < SV.size(); i++) {
			SV[i].Knots_Refine_Num(1);
		}*/
		ps.outPutVTK(SV, "E:\\kuang_models\\CarpartVolume.vtk");
		rwg.WriteSplineVolume("E:\\kuang_models\\CarpartVolume.txt", SV);
		//rwg.WriteNurbsSurface("E:\\kuang_models\\CarpartVolume.txt", YNSS1);
		//rwg.WriteSplineSurface("E:\\kuang_models\\CarpartVolume.txt", YNSS1);
		//rwg.WriteNurbsVol("E:\\kuang_models\\CarpartVolume.txt", NVtemp);
	}
};


//����
class MachineTool {
private:
	//����
	double x = (1.1*1.1*11.9) / 92;

	//y�ᴲ��
	double L1 = 80 * x, L2 = 25 * x, L3 = 15 * x, H1 = 40 * x, H2 = 12 * x;//����ߴ�
	double M1 = 200 * x, M2 = 50 * x, M3 = 18 * x, M4 = 18 * x, M5 = 15 * x, M6 = 15 * x;//�����ȵ�

	//����
	double a1 = 180 * x, b1 = 100 * x, b2 = 40 * x, c1 = 240 * x;
	double a2 = 2 * M3 / 3;

	//֧�ż�
	double c2 = c1 / 30;//���쳤��

	//���Ӳ��ֲ���
	double u = 1.1*1.1*11.9;//��
	double v = 10;//��
	double w = 12;//��

	Model_Solution m;
	RWGeometric rwg;

public:
	//y�ᴲ��
	varray<SplineVolume> machineTool1() {
		/*
			����ʹ��PublicModels�е�RandomModel���������ı���
			������Ը��ģ��Լ��ٴ����ظ���
		*/
		varray<Spline> SL1;
		varray<Spline> SL2;
		varray<Spline> SL3;
		varray<Spline> SL4;
		varray<Spline> SL5;
		varray<Spline> SL6;
		varray<Spline> SL7;
		varray<Spline> SL8;
		varray<Spline> SL9;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		SL1.resize(4);
		SL2.resize(4);
		SL3.resize(4);
		SL4.resize(4);
		SL5.resize(4);
		SL6.resize(4);
		SL7.resize(4);
		SL8.resize(4);
		SL9.resize(4);
		for (int i = 0; i < 4; i++) {
			SL1[i].m_Degree = 2;
			SL1[i].m_Knots = knots;
			SL2[i].m_Degree = 2;
			SL2[i].m_Knots = knots;
			SL3[i].m_Degree = 2;
			SL3[i].m_Knots = knots;
			SL4[i].m_Degree = 2;
			SL4[i].m_Knots = knots;
			SL5[i].m_Degree = 2;
			SL5[i].m_Knots = knots;
			SL6[i].m_Degree = 2;
			SL6[i].m_Knots = knots;
			SL7[i].m_Degree = 2;
			SL7[i].m_Knots = knots;
			SL8[i].m_Degree = 2;
			SL8[i].m_Knots = knots;
			SL9[i].m_Degree = 2;
			SL9[i].m_Knots = knots;
		}
		//����������������
		Vec4 v1 = { -L1 / 2,-H1 / 2,0,1 };
		Vec4 v2 = { -L1 / 2,0,0,1 };
		Vec4 v3 = { -L1 / 4,0,0,1 };
		Vec4 v4 = { -L1 / 4,H2,0,1 };
		Vec4 v5 = { -L2 / 2 - L3,H2,0,1 };
		Vec4 v6 = { -L2 / 2 - L3,H1 / 2,0,1 };
		Vec4 v7 = { -L2 / 2 ,H1 / 2,0,1 };
		Vec4 v8 = { -L2 / 2 ,0,0,1 };
		Vec4 v9 = { L2 / 2 ,0,0,1 };
		Vec4 v10 = { L2 / 2 ,H1 / 2,0,1 };
		Vec4 v11 = { L2 / 2 + L3,H1 / 2,0,1 };
		Vec4 v12 = { L2 / 2 + L3,H2,0,1 };
		Vec4 v13 = { L1 / 4,H2,0,1 };
		Vec4 v14 = { L1 / 4,0,0,1 };
		Vec4 v15 = { L1 / 2,0,0,1 };
		Vec4 v16 = { L1 / 2,-H1 / 2,0,1 };
		Vec4 v17 = { L2 / 2,-H1 / 2,0,1 };
		Vec4 v18 = { -L2 / 2,-H1 / 2,0,1 };
		Vec4 v19 = (v1 + v18) / 2;
		Vec4 v20 = (v17 + v16) / 2;

		SL1[0].m_CtrlPts.push_back(v1);
		SL1[0].m_CtrlPts.push_back((v1 + v19) / 2);
		SL1[0].m_CtrlPts.push_back(v19);
		SL1[1].m_CtrlPts.push_back(v1);
		SL1[1].m_CtrlPts.push_back((v1 + v2) / 2);
		SL1[1].m_CtrlPts.push_back(v2);
		SL1[2].m_CtrlPts.push_back(v2);
		SL1[2].m_CtrlPts.push_back((v2 + v3) / 2);
		SL1[2].m_CtrlPts.push_back(v3);
		SL1[3].m_CtrlPts.push_back(v19);
		SL1[3].m_CtrlPts.push_back((v19 + v3) / 2);
		SL1[3].m_CtrlPts.push_back(v3);
		SplineSurface ss1;
		ss1.CoonsInterpolate(SL1);

		SL2[0].m_CtrlPts.push_back(v19);
		SL2[0].m_CtrlPts.push_back((v18 + v19) / 2);
		SL2[0].m_CtrlPts.push_back(v18);
		SL2[1].m_CtrlPts.push_back(v19);
		SL2[1].m_CtrlPts.push_back((v19 + v3) / 2);
		SL2[1].m_CtrlPts.push_back(v3);
		SL2[2].m_CtrlPts.push_back(v3);
		SL2[2].m_CtrlPts.push_back((v8 + v3) / 2);
		SL2[2].m_CtrlPts.push_back(v8);
		SL2[3].m_CtrlPts.push_back(v18);
		SL2[3].m_CtrlPts.push_back((v18 + v8) / 2);
		SL2[3].m_CtrlPts.push_back(v8);
		SplineSurface ss2;
		ss2.CoonsInterpolate(SL2);

		SL3[0].m_CtrlPts.push_back(v18);
		SL3[0].m_CtrlPts.push_back((v18 + v17) / 2);
		SL3[0].m_CtrlPts.push_back(v17);
		SL3[1].m_CtrlPts.push_back(v18);
		SL3[1].m_CtrlPts.push_back((v18 + v8) / 2);
		SL3[1].m_CtrlPts.push_back(v8);
		SL3[2].m_CtrlPts.push_back(v8);
		SL3[2].m_CtrlPts.push_back((v8 + v9) / 2);
		SL3[2].m_CtrlPts.push_back(v9);
		SL3[3].m_CtrlPts.push_back(v17);
		SL3[3].m_CtrlPts.push_back((v17 + v9) / 2);
		SL3[3].m_CtrlPts.push_back(v9);
		SplineSurface ss3;
		ss3.CoonsInterpolate(SL3);

		SL4[0].m_CtrlPts.push_back(v17);
		SL4[0].m_CtrlPts.push_back((v17 + v20) / 2);
		SL4[0].m_CtrlPts.push_back(v20);
		SL4[1].m_CtrlPts.push_back(v17);
		SL4[1].m_CtrlPts.push_back((v17 + v9) / 2);
		SL4[1].m_CtrlPts.push_back(v9);
		SL4[2].m_CtrlPts.push_back(v9);
		SL4[2].m_CtrlPts.push_back((v9 + v14) / 2);
		SL4[2].m_CtrlPts.push_back(v14);
		SL4[3].m_CtrlPts.push_back(v20);
		SL4[3].m_CtrlPts.push_back((v20 + v14) / 2);
		SL4[3].m_CtrlPts.push_back(v14);
		SplineSurface ss4;
		ss4.CoonsInterpolate(SL4);

		SL5[0].m_CtrlPts.push_back(v20);
		SL5[0].m_CtrlPts.push_back((v20 + v16) / 2);
		SL5[0].m_CtrlPts.push_back(v16);
		SL5[1].m_CtrlPts.push_back(v20);
		SL5[1].m_CtrlPts.push_back((v20 + v14) / 2);
		SL5[1].m_CtrlPts.push_back(v14);
		SL5[2].m_CtrlPts.push_back(v14);
		SL5[2].m_CtrlPts.push_back((v14 + v15) / 2);
		SL5[2].m_CtrlPts.push_back(v15);
		SL5[3].m_CtrlPts.push_back(v16);
		SL5[3].m_CtrlPts.push_back((v16 + v15) / 2);
		SL5[3].m_CtrlPts.push_back(v15);
		SplineSurface ss5;
		ss5.CoonsInterpolate(SL5);

		SL6[0].m_CtrlPts.push_back(v3);
		SL6[0].m_CtrlPts.push_back((v3 + v8) / 2);
		SL6[0].m_CtrlPts.push_back(v8);
		SL6[1].m_CtrlPts.push_back(v3);
		SL6[1].m_CtrlPts.push_back((v3 + v4) / 2);
		SL6[1].m_CtrlPts.push_back(v4);
		SL6[2].m_CtrlPts.push_back(v4);
		SL6[2].m_CtrlPts.push_back((v4 + v7) / 2);
		SL6[2].m_CtrlPts.push_back(v7);
		SL6[3].m_CtrlPts.push_back(v8);
		SL6[3].m_CtrlPts.push_back((v8 + v7) / 2);
		SL6[3].m_CtrlPts.push_back(v7);
		SplineSurface ss6;
		ss6.CoonsInterpolate(SL6);

		SL7[0].m_CtrlPts.push_back(v5);
		SL7[0].m_CtrlPts.push_back((v5 + v4) / 2);
		SL7[0].m_CtrlPts.push_back(v4);
		SL7[1].m_CtrlPts.push_back(v5);
		SL7[1].m_CtrlPts.push_back((v5 + v6) / 2);
		SL7[1].m_CtrlPts.push_back(v6);
		SL7[2].m_CtrlPts.push_back(v6);
		SL7[2].m_CtrlPts.push_back((v6 + v7) / 2);
		SL7[2].m_CtrlPts.push_back(v7);
		SL7[3].m_CtrlPts.push_back(v4);
		SL7[3].m_CtrlPts.push_back((v7 + v4) / 2);
		SL7[3].m_CtrlPts.push_back(v7);
		SplineSurface ss7;
		ss7.CoonsInterpolate(SL7);

		SL8[0].m_CtrlPts.push_back(v9);
		SL8[0].m_CtrlPts.push_back((v14 + v9) / 2);
		SL8[0].m_CtrlPts.push_back(v14);
		SL8[1].m_CtrlPts.push_back(v9);
		SL8[1].m_CtrlPts.push_back((v9 + v10) / 2);
		SL8[1].m_CtrlPts.push_back(v10);
		SL8[2].m_CtrlPts.push_back(v10);
		SL8[2].m_CtrlPts.push_back((v10 + v13) / 2);
		SL8[2].m_CtrlPts.push_back(v13);
		SL8[3].m_CtrlPts.push_back(v14);
		SL8[3].m_CtrlPts.push_back((v14 + v13) / 2);
		SL8[3].m_CtrlPts.push_back(v13);
		SplineSurface ss8;
		ss8.CoonsInterpolate(SL8);

		SL9[0].m_CtrlPts.push_back(v13);
		SL9[0].m_CtrlPts.push_back((v13 + v12) / 2);
		SL9[0].m_CtrlPts.push_back(v12);
		SL9[1].m_CtrlPts.push_back(v13);
		SL9[1].m_CtrlPts.push_back((v10 + v13) / 2);
		SL9[1].m_CtrlPts.push_back(v10);
		SL9[2].m_CtrlPts.push_back(v10);
		SL9[2].m_CtrlPts.push_back((v10 + v11) / 2);
		SL9[2].m_CtrlPts.push_back(v11);
		SL9[3].m_CtrlPts.push_back(v12);
		SL9[3].m_CtrlPts.push_back((v12 + v11) / 2);
		SL9[3].m_CtrlPts.push_back(v11);
		SplineSurface ss9;
		ss9.CoonsInterpolate(SL9);

		varray<SplineSurface> SS;
		SS.push_back(ss1);
		SS.push_back(ss2);
		SS.push_back(ss3);
		SS.push_back(ss4);
		SS.push_back(ss5);
		SS.push_back(ss6);
		SS.push_back(ss7);
		SS.push_back(ss8);
		SS.push_back(ss9);

		Model_Solution m;
		double M7 = M1 / 2 - M2 / 2 - M3 - M4 - M5;
		varray<SplineVolume> SV;
		varray<SplineVolume> SV1;
		varray<SplineVolume> SV2;
		varray<SplineVolume> SV3;
		varray<SplineVolume> SV4;
		varray<SplineVolume> SV5;
		varray<SplineVolume> SV6;
		varray<SplineVolume> SV7;
		varray<SplineVolume> SV8;
		varray<SplineVolume> SV9;
		SV1 = m.CreatSweepVol(SS, M7, 3);
		m.Trans(SS, M7, 3);
		SV2 = m.CreatSweepVol(SS, M5, 3);
		m.Trans(SS, M5, 3);
		SV3 = m.CreatSweepVol(SS, M4, 3);
		m.Trans(SS, M4, 3);
		SV4 = m.CreatSweepVol(SS, M3, 3);
		m.Trans(SS, M3, 3);
		SV5 = m.CreatSweepVol(SS, M2, 3);
		m.Trans(SS, M2, 3);
		SV6 = m.CreatSweepVol(SS, M3, 3);
		m.Trans(SS, M3, 3);
		SV7 = m.CreatSweepVol(SS, M4, 3);
		m.Trans(SS, M4, 3);
		SV8 = m.CreatSweepVol(SS, M5, 3);
		m.Trans(SS, M5, 3);
		SV9 = m.CreatSweepVol(SS, M7, 3);
		m.Trans(SS, M7, 3);
		for (int i = 0; i < 9; i++) {
			SV.push_back(SV1[i]);
			SV.push_back(SV2[i]);
			SV.push_back(SV3[i]);
			SV.push_back(SV4[i]);
			SV.push_back(SV5[i]);
			SV.push_back(SV6[i]);
			SV.push_back(SV7[i]);
			SV.push_back(SV8[i]);
			SV.push_back(SV9[i]);
		}

		////rwg.WriteSplineSurface("E:\\kuang_models\\machineToolSurface.txt", SS);
		//rwg.WriteSplineVolume("E:\\kuang_models\\machineToolVolume.txt", SV);

		return SV;
	}
	varray<SplineVolume> machineTool2() {
		/*
			����ʹ��PublicModels�е�RandomModel���������ı���
			������Ը��ģ��Լ��ٴ����ظ���
		*/
		varray<Spline> SL1;
		varray<Spline> SL2;
		varray<Spline> SL3;
		varray<Spline> SL4;
		varray<Spline> SL5;
		varray<Spline> SL6;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		SL1.resize(4);
		SL2.resize(4);
		SL3.resize(4);
		SL4.resize(4);
		SL5.resize(4);
		SL6.resize(4);

		for (int i = 0; i < 4; i++) {
			SL1[i].m_Degree = 2;
			SL1[i].m_Knots = knots;
			SL2[i].m_Degree = 2;
			SL2[i].m_Knots = knots;
			SL3[i].m_Degree = 2;
			SL3[i].m_Knots = knots;
			SL4[i].m_Degree = 2;
			SL4[i].m_Knots = knots;
			SL5[i].m_Degree = 2;
			SL5[i].m_Knots = knots;
			SL6[i].m_Degree = 2;
			SL6[i].m_Knots = knots;
		}
		//����������������
		Vec4 v1 = { -M2 / 2 - M3 - M4 - M5,-M6,0,1 };
		Vec4 v2 = { -M2 / 2 - M3 - M4 - M5,0,0,1 };
		Vec4 v3 = { -M2 / 2 - M3 - M4,0,0,1 };
		Vec4 v4 = { -M2 / 2 - M3 - M4,-M6,0,1 };
		Vec4 v5 = { -M2 / 2 - M3,-M6 + a2,0,1 };
		Vec4 v6 = { -M2 / 2 - M3,0,0,1 };
		Vec4 v7 = { -M2 / 2 ,0,0,1 };
		Vec4 v8 = { -M2 / 2 - 2 * M3 / 3 ,-M6 + a2,0,1 };
		Vec4 v9 = { -M2 / 2 - 2 * M3 / 3 ,-M6,0,1 };
		Vec4 v10 = { -M2 / 2 ,-M6 ,0,1 };
		Vec4 v11 = { M2 / 2 ,-M6 ,0,1 };
		Vec4 v12 = { M2 / 2 ,0,0,1 };
		Vec4 v13 = { M2 / 2 + M3,0,0,1 };
		Vec4 v14 = { M2 / 2 + M3,-M6 + a2,0,1 };
		Vec4 v15 = { M2 / 2 + 2 * M3 / 3 ,-M6 + a2,0,1 };
		Vec4 v16 = { M2 / 2 + 2 * M3 / 3 ,-M6,0,1 };
		Vec4 v17 = { M2 / 2 + M3 + M4,-M6,0,1 };
		Vec4 v18 = { M2 / 2 + M3 + M4,0,0,1 };
		Vec4 v19 = { M2 / 2 + M3 + M4 + M5,0,0,1 };
		Vec4 v20 = { M2 / 2 + M3 + M4 + M5,-M6,0,1 };



		SL1[0].m_CtrlPts.push_back(v1);
		SL1[0].m_CtrlPts.push_back((v1 + v4) / 2);
		SL1[0].m_CtrlPts.push_back(v4);
		SL1[1].m_CtrlPts.push_back(v1);
		SL1[1].m_CtrlPts.push_back((v1 + v2) / 2);
		SL1[1].m_CtrlPts.push_back(v2);
		SL1[2].m_CtrlPts.push_back(v2);
		SL1[2].m_CtrlPts.push_back((v2 + v3) / 2);
		SL1[2].m_CtrlPts.push_back(v3);
		SL1[3].m_CtrlPts.push_back(v4);
		SL1[3].m_CtrlPts.push_back((v4 + v3) / 2);
		SL1[3].m_CtrlPts.push_back(v3);
		SplineSurface ss1;
		ss1.CoonsInterpolate(SL1);

		SL2[0].m_CtrlPts.push_back(v5);
		SL2[0].m_CtrlPts.push_back((v5 + v8) / 2);
		SL2[0].m_CtrlPts.push_back(v8);
		SL2[1].m_CtrlPts.push_back(v5);
		SL2[1].m_CtrlPts.push_back((v5 + v6) / 2);
		SL2[1].m_CtrlPts.push_back(v6);
		SL2[2].m_CtrlPts.push_back(v6);
		SL2[2].m_CtrlPts.push_back((v6 + v7) / 2);
		SL2[2].m_CtrlPts.push_back(v7);
		SL2[3].m_CtrlPts.push_back(v8);
		SL2[3].m_CtrlPts.push_back((v7 + v8) / 2);
		SL2[3].m_CtrlPts.push_back(v7);
		SplineSurface ss2;
		ss2.CoonsInterpolate(SL2);

		SL3[0].m_CtrlPts.push_back(v9);
		SL3[0].m_CtrlPts.push_back((v9 + v10) / 2);
		SL3[0].m_CtrlPts.push_back(v10);
		SL3[1].m_CtrlPts.push_back(v9);
		SL3[1].m_CtrlPts.push_back((v9 + v8) / 2);
		SL3[1].m_CtrlPts.push_back(v8);
		SL3[2].m_CtrlPts.push_back(v8);
		SL3[2].m_CtrlPts.push_back((v8 + v7) / 2);
		SL3[2].m_CtrlPts.push_back(v7);
		SL3[3].m_CtrlPts.push_back(v10);
		SL3[3].m_CtrlPts.push_back((v10 + v7) / 2);
		SL3[3].m_CtrlPts.push_back(v7);
		SplineSurface ss3;
		ss3.CoonsInterpolate(SL3);

		SL4[0].m_CtrlPts.push_back(v11);
		SL4[0].m_CtrlPts.push_back((v11 + v16) / 2);
		SL4[0].m_CtrlPts.push_back(v16);
		SL4[1].m_CtrlPts.push_back(v11);
		SL4[1].m_CtrlPts.push_back((v11 + v12) / 2);
		SL4[1].m_CtrlPts.push_back(v12);
		SL4[2].m_CtrlPts.push_back(v12);
		SL4[2].m_CtrlPts.push_back((v12 + v15) / 2);
		SL4[2].m_CtrlPts.push_back(v15);
		SL4[3].m_CtrlPts.push_back(v16);
		SL4[3].m_CtrlPts.push_back((v16 + v15) / 2);
		SL4[3].m_CtrlPts.push_back(v15);
		SplineSurface ss4;
		ss4.CoonsInterpolate(SL4);

		SL5[0].m_CtrlPts.push_back(v15);
		SL5[0].m_CtrlPts.push_back((v15 + v14) / 2);
		SL5[0].m_CtrlPts.push_back(v14);
		SL5[1].m_CtrlPts.push_back(v15);
		SL5[1].m_CtrlPts.push_back((v15 + v12) / 2);
		SL5[1].m_CtrlPts.push_back(v12);
		SL5[2].m_CtrlPts.push_back(v12);
		SL5[2].m_CtrlPts.push_back((v12 + v13) / 2);
		SL5[2].m_CtrlPts.push_back(v13);
		SL5[3].m_CtrlPts.push_back(v14);
		SL5[3].m_CtrlPts.push_back((v14 + v13) / 2);
		SL5[3].m_CtrlPts.push_back(v13);
		SplineSurface ss5;
		ss5.CoonsInterpolate(SL5);

		SL6[0].m_CtrlPts.push_back(v17);
		SL6[0].m_CtrlPts.push_back((v17 + v20) / 2);
		SL6[0].m_CtrlPts.push_back(v20);
		SL6[1].m_CtrlPts.push_back(v17);
		SL6[1].m_CtrlPts.push_back((v17 + v18) / 2);
		SL6[1].m_CtrlPts.push_back(v18);
		SL6[2].m_CtrlPts.push_back(v18);
		SL6[2].m_CtrlPts.push_back((v18 + v19) / 2);
		SL6[2].m_CtrlPts.push_back(v19);
		SL6[3].m_CtrlPts.push_back(v20);
		SL6[3].m_CtrlPts.push_back((v20 + v19) / 2);
		SL6[3].m_CtrlPts.push_back(v19);
		SplineSurface ss6;
		ss6.CoonsInterpolate(SL6);

		varray<SplineSurface> SS;
		SS.push_back(ss1);
		SS.push_back(ss2);
		SS.push_back(ss3);
		SS.push_back(ss4);
		SS.push_back(ss5);
		SS.push_back(ss6);

		Model_Solution m;
		varray<SplineVolume> SV;
		SV = m.CreatSweepVol(SS, L1, 3);

		////rwg.WriteSplineSurface("E:\\kuang_models\\machineToolSurface.txt", SS);
		//rwg.WriteSplineVolume("E:\\kuang_models\\machineToolVolume.txt", SV);

		return SV;
	}

	//֧�ż�
	varray<SplineVolume> machineTool4() {
		double d1 = b2 / 2, d2 = b2 / 3;
		double d3 = 2 * a2 / 3;

		Vec4 v1 = { -d2,0,0,1 };
		Vec4 v2 = { -d2,d1 / 4,0,1 };
		Vec4 v3 = { -d2,d1 / 2,0,1 };
		Vec4 v4 = { -d2,d1,0,1 };
		Vec4 v5 = { 0,b2,0,1 };
		Vec4 v6 = { -d2 / 3,d1 / 2,0,1 };
		Vec4 v7 = { -d2 / 3,d1 / 4,0,1 };
		Vec4 v8 = { 0,0,0,1 };

		varray<SplineSurface> SS;
		SplineSurface ss;
		SplineSurface ss1;
		RandomModel rm1(v1, v2, v7, v8);
		ss = rm1.getSurface();
		SS.push_back(ss);

		RandomModel rm2(v2, v3, v6, v7);
		ss = rm2.getSurface();
		ss1 = ss;
		SS.push_back(ss);

		RandomModel rm3(v3, v4, v5, v6);
		ss = rm3.getSurface();
		SS.push_back(ss);

		RandomModel rm4(v7, v6, v5, v8);
		ss = rm4.getSurface();
		SS.push_back(ss);

		varray<SplineVolume> SV;
		SV = m.CreatSweepVol(SS, c2, 3);

		varray<SplineVolume> SV1;
		SV1 = SV;
		m.Trans(SV1, c2 * 4, 3);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}

		SS.clear();
		SS.push_back(ss1);
		varray<SplineVolume> SV2;
		SV2 = m.CreatSweepVol(SS, c2 * 3, 3);
		m.Trans(SV2, c2, 3);
		for (int i = 0; i < SV2.size(); i++) {
			SV.push_back(SV2[i]);
		}

		SV1 = SV;
		SV.clear();
		m.Trans(SV1, a1 / 2, -1);
		m.Trans(SV1, 5 * c2, -3);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		m.Trans(SV1, 3 * c1 / 4 - 5 * c2, -3);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = SV;
		m.Rolate(SV1, PI, 2);
		m.Trans(SV1, 3 * c1 / 4, -3);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		m.Trans(SV, 3 * c1 / 4, 3);

		rwg.WriteSplineVolume("E:\\kuang_models\\supportpartVolume.txt", SV);

		return SV;
	}

	//����
	varray<SplineVolume> machineTool3() {
		/*
			����ʹ��PublicModels�е�RandomModel���������ı���
			������Ը��ģ��Լ��ٴ����ظ���
		*/
		varray<Spline> SL1;
		varray<Spline> SL2;
		varray<double> knots;
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(0);
		knots.push_back(1);
		knots.push_back(1);
		knots.push_back(1);
		SL1.resize(4);
		SL2.resize(4);
		for (int i = 0; i < 4; i++) {
			SL1[i].m_Degree = 2;
			SL1[i].m_Knots = knots;
			SL2[i].m_Degree = 2;
			SL2[i].m_Knots = knots;
		}
		//������������
		Vec4 v1 = { -a1 / 2 ,0,0,1 };
		Vec4 v2 = { -a1 / 2 ,b2,0,1 };
		Vec4 v3 = { -M2 / 2 - a2 * 3,b1 - a2 * 2,0,1 };
		Vec4 v4 = { -M2 / 2 - a2 * 3,0,0,1 };
		Vec4 v5 = { M2 / 2 + a2 * 3,0,0,1 };
		Vec4 v6 = { M2 / 2 + a2 * 3,b1 - a2 * 2,0,1 };
		Vec4 v7 = { a1 / 2 ,b2,0,1 };
		Vec4 v8 = { a1 / 2 ,0,0,1 };


		Vec4 v9 = { -M2 / 2 - a2 * 2,0,0,1 };
		Vec4 v10 = { -M2 / 2 - a2 * 2 ,a2,0,1 };
		Vec4 v11 = { -M2 / 2 - a2 * 3 ,a2,0,1 };
		Vec4 v12 = { -M2 / 2 - a2 * 3,M6,0,1 };
		Vec4 v13 = { -M2 / 2 - M3,M6,0,1 };
		Vec4 v14 = { -M2 / 2 - M3,a2,0,1 };
		Vec4 v15 = { -M2 / 2 - a2,a2,0,1 };
		Vec4 v16 = { -M2 / 2 - a2,0,0,1 };
		Vec4 v17 = { M2 / 2 + a2,0,0,1 };
		Vec4 v18 = { M2 / 2 + a2,a2,0,1 };
		Vec4 v19 = { M2 / 2 + M3,a2,0,1 };
		Vec4 v20 = { M2 / 2 + M3,M6,0,1 };
		Vec4 v21 = { M2 / 2 + a2 * 3,M6,0,1 };
		Vec4 v22 = { M2 / 2 + a2 * 3 ,a2,0,1 };
		Vec4 v23 = { M2 / 2 + a2 * 2 ,a2,0,1 };
		Vec4 v24 = { M2 / 2 + a2 * 2,0,0,1 };

		//������ײ㲿��
		varray<SplineSurface> SS;
		varray<SplineSurface> SS3;//��������������Ƭ
		varray<SplineSurface> ss1;
		varray<SplineSurface> ss2;
		Rectangle0 rt(a2, a2);
		SplineSurface ss;
		ss = rt.getSurface();
		m.Trans(ss, M2 / 2 + a2 + a2 / 2, -1);
		m.Trans(ss, b1 + a2 / 2, 2);
		ss1.push_back(ss);
		m.Trans(ss, a2, 1);
		ss1.push_back(ss);
		m.Trans(ss, M2 + a2, 1);
		ss1.push_back(ss);
		m.Trans(ss, a2, 1);
		ss1.push_back(ss);
		m.Trans(ss1, a2, -2);
		for (int i = 0; i < ss1.size(); i++) {
			SS.push_back(ss1[i]);
		}
		ss2 = ss1;
		m.Trans(ss1, a2, -2);
		for (int i = 0; i < ss1.size(); i++) {
			SS.push_back(ss1[i]);
		}

		Rectangle0 rt1(a2, b1 - a2 * 2);
		ss = rt1.getSurface();
		m.Trans(ss, M2 / 2 + 5 * a2 / 2, -1);
		m.Trans(ss, (b1 - a2 * 2) / 2, 2);
		SS.push_back(ss);
		SS3.push_back(ss);
		m.Trans(ss, a2, 1);
		SS.push_back(ss);
		SS3.push_back(ss);
		m.Trans(ss, a2, 1);
		SS.push_back(ss);
		SS3.push_back(ss);
		m.Trans(ss, a2 + M2, 1);
		SS.push_back(ss);
		SS3.push_back(ss);
		m.Trans(ss, a2, 1);
		SS.push_back(ss);
		SS3.push_back(ss);
		m.Trans(ss, a2, 1);
		SS.push_back(ss);
		SS3.push_back(ss);

		Rectangle0 rt2(M2, b1 - a2 * 2);
		ss = rt2.getSurface();
		m.Trans(ss, (b1 - a2 * 2) / 2, 2);
		SS.push_back(ss);
		SS3.push_back(ss);

		Rectangle0 rt3(M2, a2);
		ss = rt3.getSurface();
		m.Trans(ss, b1 - a2 * 2 + a2 / 2, 2);
		SS.push_back(ss);
		SplineSurface s;//������
		s = ss;

		SL1[0].m_CtrlPts.push_back(v1);
		SL1[0].m_CtrlPts.push_back((v1 + v4) / 2);
		SL1[0].m_CtrlPts.push_back(v4);
		SL1[1].m_CtrlPts.push_back(v1);
		SL1[1].m_CtrlPts.push_back((v1 + v2) / 2);
		SL1[1].m_CtrlPts.push_back(v2);
		SL1[2].m_CtrlPts.push_back(v2);
		SL1[2].m_CtrlPts.push_back((v2 + v3) / 2);
		SL1[2].m_CtrlPts.push_back(v3);
		SL1[3].m_CtrlPts.push_back(v4);
		SL1[3].m_CtrlPts.push_back((v4 + v3) / 2);
		SL1[3].m_CtrlPts.push_back(v3);
		SplineSurface ss3;
		ss3.CoonsInterpolate(SL1);
		SS.push_back(ss3);
		SS3.push_back(ss3);

		SL2[0].m_CtrlPts.push_back(v5);
		SL2[0].m_CtrlPts.push_back((v5 + v8) / 2);
		SL2[0].m_CtrlPts.push_back(v8);
		SL2[1].m_CtrlPts.push_back(v5);
		SL2[1].m_CtrlPts.push_back((v5 + v6) / 2);
		SL2[1].m_CtrlPts.push_back(v6);
		SL2[2].m_CtrlPts.push_back(v6);
		SL2[2].m_CtrlPts.push_back((v6 + v7) / 2);
		SL2[2].m_CtrlPts.push_back(v7);
		SL2[3].m_CtrlPts.push_back(v8);
		SL2[3].m_CtrlPts.push_back((v8 + v7) / 2);
		SL2[3].m_CtrlPts.push_back(v7);
		SplineSurface ss4;
		ss4.CoonsInterpolate(SL2);
		SS.push_back(ss4);
		SS3.push_back(ss4);

		//�м��οղ���
		varray<SplineSurface> SS1;
		SS1.push_back(ss2[1]);
		SS1.push_back(ss2[2]);
		m.Trans(SS1, a2, 2);
		m.Trans(s, a2, 2);
		SS1.push_back(s);
		m.Trans(s, a2, 2);
		SS1.push_back(s);

		Rectangle0 rt4(M3 / 3, M6 - a2);
		ss = rt4.getSurface();
		m.Trans(ss, b1 + a2 + (M6 - a2) / 2, 2);
		m.Trans(ss, M2 / 2 + a2 + M3 / 6, -1);
		SS1.push_back(ss);
		m.Trans(ss, 2 * (M2 / 2 + a2 + M3 / 6), 1);
		SS1.push_back(ss);

		Rectangle0 rt5(a2, M6 - a2);
		ss = rt5.getSurface();
		m.Trans(ss, b1 + a2 + (M6 - a2) / 2, 2);
		m.Trans(ss, M2 / 2 + a2 / 2, -1);
		SS1.push_back(ss);
		m.Trans(ss, 2 * (M2 / 2 + a2 / 2), 1);
		SS1.push_back(ss);

		Rectangle0 rt6(M2, M6 - a2);
		ss = rt6.getSurface();
		m.Trans(ss, b1 + a2 + (M6 - a2) / 2, 2);
		SS1.push_back(ss);

		//����
		varray<SplineSurface> SS2;
		SplineSurface ss5;

		RandomModel rm(v17, v18, v19, v24);
		ss5 = rm.getSurface();
		SS2.push_back(ss5);
		RandomModel rm1(v19, v20, v23, v24);
		ss5 = rm1.getSurface();
		SS2.push_back(ss5);
		RandomModel rm2(v23, v20, v21, v22);
		ss5 = rm2.getSurface();
		SS2.push_back(ss5);
		RandomModel rm3(v9, v14, v15, v16);
		ss5 = rm3.getSurface();
		SS2.push_back(ss5);
		RandomModel rm4(v9, v10, v13, v14);
		ss5 = rm4.getSurface();
		SS2.push_back(ss5);
		RandomModel rm5(v11, v12, v13, v10);
		ss5 = rm5.getSurface();
		SS2.push_back(ss5);

		////�������䲿�� ��ʱ��Ҫ
		//Vec4 p1 = { a1 / 2 + a1 / 10,0,0,1 };
		//Vec4 p2 = { a1 / 2 + a1 / 10,0,0,1 };
		//Vec4 p1 = { a1 / 2 + a1 / 10,0,0,1 };
		//Vec4 p1 = { a1 / 2 + a1 / 10,0,0,1 };
		//RandomModel();



		varray<SplineVolume> SV;
		varray<SplineVolume> SV1;
		SV = m.CreatSweepVol(SS, c2, 3);
		m.Trans(SS, c2, 3);
		SV1 = m.CreatSweepVol(SS, c2 * 3, 3);
		m.Trans(SS, c2 * 3, 3);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = m.CreatSweepVol(SS, c2, 3);
		m.Trans(SS, c2, 3);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = m.CreatSweepVol(SS, c1 * 3 / 4 - 10 * c2, 3);
		m.Trans(SS, c1 * 3 / 4 - 10 * c2, 3);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = m.CreatSweepVol(SS, c2, 3);
		m.Trans(SS, c2, 3);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = m.CreatSweepVol(SS, c2 * 3, 3);
		m.Trans(SS, c2 * 3, 3);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = m.CreatSweepVol(SS, c2, 3);
		m.Trans(SS, c2, 3);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}

		SV1 = m.CreatSweepVol(SS, v, 3);//���������Ӳ���
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = m.CreatSweepVol(ss2, c1 / 10, -3);//�������
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		m.Trans(SS1, 3 * c1 / 4, 3);
		SV1 = m.CreatSweepVol(SS1, v, 3);//���������Ӳ���
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}

		SV1 = m.CreatSweepVol(SS2, c2, 3);
		m.Trans(SS2, c2, 3);
		m.Trans(SV1, b1, 2);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = m.CreatSweepVol(SS2, 3 * c2, 3);
		m.Trans(SS2, 3 * c2, 3);
		m.Trans(SV1, b1, 2);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = m.CreatSweepVol(SS2, c2, 3);
		m.Trans(SS2, c2, 3);
		m.Trans(SV1, b1, 2);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = m.CreatSweepVol(SS2, 3 * c1 / 4 - 10 * c2, 3);
		m.Trans(SS2, 3 * c1 / 4 - 10 * c2, 3);
		m.Trans(SV1, b1, 2);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = m.CreatSweepVol(SS2, c2, 3);
		m.Trans(SS2, c2, 3);
		m.Trans(SV1, b1, 2);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = m.CreatSweepVol(SS2, 3 * c2, 3);
		m.Trans(SS2, c2 * 3, 3);
		m.Trans(SV1, b1, 2);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		SV1 = m.CreatSweepVol(SS2, c2, 3);
		m.Trans(SS2, c2, 3);
		m.Trans(SV1, b1, 2);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}

		SV1 = machineTool4();
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}

		/*
		rwg.WriteSplineSurface("E:\\kuang_models\\newSurface.txt", SS);
		rwg.WriteSplineSurface("E:\\kuang_models\\newSurface1.txt", SS1);
		rwg.WriteSplineSurface("E:\\kuang_models\\newSurface2.txt", SS2);*/
		rwg.WriteSplineVolume("E:\\kuang_models\\dizuovolume.txt", SV);

		return SV;
	}

	//y�ᴲ��
	varray<SplineVolume> getVolume1() {
		varray<SplineVolume> SV;
		varray<SplineVolume> SV1;
		varray<SplineVolume> SV2;
		SV1 = machineTool1();
		SV2 = machineTool2();
		m.Rolate(SV2, PI / 2, 2);
		m.Trans(SV2, H1 / 2, -2);
		m.Trans(SV2, M1 / 2, 3);
		m.Trans(SV2, L1 / 2, -1);
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		for (int i = 0; i < SV2.size(); i++) {
			SV.push_back(SV2[i]);
		}

		rwg.WriteSplineVolume("E:\\kuang_models\\yzhouchuangshenVolume.txt", SV);

		return SV;

	}

	//y�ᴲ��,������֧�ż�
	varray<SplineVolume> getVolume2() {
		varray<SplineVolume> SV;
		varray<SplineVolume> SV1;
		SV = getVolume1();
		m.Trans(SV, M1 / 2, -3);
		m.Rolate(SV, PI / 2, 2);
		m.Trans(SV, b1 + M6 + H1 / 2, 2);
		m.Trans(SV, c1 / 2, 3);
		SV1 = machineTool3();
		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}

		m.Rolate(SV, PI, 2);
		m.Trans(SV, b1 + M6 + w, -2);
		m.Trans(SV, c1 * 3 / 4, 3);
		m.Trans(SV, 11.0 / 2, 1);
		rwg.WriteSplineVolume("E:\\kuang_models\\bothmachinevolume.txt", SV);

		return SV;
	}

	//y�ᴲ������������
	varray<SplineVolume> getVolume3() {
		varray<SplineVolume> SV;
		varray<SplineVolume> SV1;
		MachineTool mt;
		SV = mt.getVolume2();

		rwg.ReadSplineVolume("E:\\kuang_models\\buningyuan\\MachineTool\\jichuang.txt", SV1);

		for (int i = 0; i < SV1.size(); i++) {
			SV.push_back(SV1[i]);
		}
		rwg.WriteSplineVolume("E:\\kuang_models\\buningyuan\\MachineTool\\machinetool.txt", SV);
		return SV;
	}

	//�ʷֳ�������ģ�ͣ����м����ԣ�y�ᴲ�����
	void testQuad() {
		//�����
		Vec4 v1 = { 0,0,0,1 };
		Vec4 v2 = { 0,H1 / 2,0,1 };
		Vec4 v3 = { L1 / 4,H1 / 2,0,1 };
		Vec4 v4 = { L1 / 4,H1 / 2 + H2,0,1 };
		Vec4 v5 = { L1 / 2 - L2 / 2 - L3 ,H1 / 2 + H2,0,1 };
		Vec4 v6 = { L1 / 2 - L2 / 2 - L3,H1,0,1 };
		Vec4 v7 = { L1 / 2 - L2 / 2,H1,0,1 };
		Vec4 v8 = { L1 / 2 - L2 / 2,H1 / 2,0,1 };
		Vec4 v9 = { L1 / 2 + L2 / 2 ,H1 / 2,0,1 };
		Vec4 v10 = { L1 / 2 + L2 / 2,H1,0,1 };
		Vec4 v11 = { L1 / 2 + L2 / 2 + L3,H1,0,1 };
		Vec4 v12 = { L1 / 2 + L2 / 2 + L3,H1 - (H1 / 2 - H2),0,1 };
		Vec4 v13 = { L1 / 2 + L1 / 4,H1 - (H1 / 2 - H2),0,1 };
		Vec4 v14 = { L1 / 2 + L1 / 4,H1 / 2,0,1 };
		Vec4 v15 = { L1,H1 / 2,0,1 };
		Vec4 v16 = { L1,0,0,1 };


		varray<Spline> sps;
		Spline0 sl;
		Spline sp;
		sp = sl.getSpline(v1, v2);
		sps.push_back(sp);

		sp = sl.getSpline(v2, v3);
		sps.push_back(sp);

		sp = sl.getSpline(v3, v4);
		sps.push_back(sp);

		sp = sl.getSpline(v4, v5);
		sps.push_back(sp);

		sp = sl.getSpline(v5, v6);
		sps.push_back(sp);

		sp = sl.getSpline(v6, v7);
		sps.push_back(sp);


		sp = sl.getSpline(v7, v8);
		sps.push_back(sp);

		sp = sl.getSpline(v8, v9);
		sps.push_back(sp);

		sp = sl.getSpline(v9, v10);
		sps.push_back(sp);

		sp = sl.getSpline(v10, v11);
		sps.push_back(sp);

		sp = sl.getSpline(v11, v12);
		sps.push_back(sp);

		sp = sl.getSpline(v12, v13);
		sps.push_back(sp);

		sp = sl.getSpline(v13, v14);
		sps.push_back(sp);

		sp = sl.getSpline(v14, v15);
		sps.push_back(sp);

		sp = sl.getSpline(v15, v16);
		sps.push_back(sp);

		sp = sl.getSpline(v16, v1);
		sps.push_back(sp);

		RWGeometric rwg;
		rwg.WriteSpline("E:\\kuang_models\\outSpline.txt", sps);


		varray<Spline>outer;//�������������
		varray<Spline>inner1;//�������������
		varray<varray<Spline>> inner;
		varray<Spline>addlines, allLines;//������ ��������������������Ϊ�����
		varray<bool> genus;
		varray<SplineSurface> allSurf;//����ʷֽ��

		rwg.ReadSpline("E:\\kuang_models\\outSpline.txt", outer);
		genus.resize(1);
		genus[0] = false;//��֮ǰ�õ���push_back,���������
		PublicSolution ps;
		ps.quad(outer, inner, addlines, genus, allSurf);

		rwg.WriteSplineSurface("E:\\kuang_models\\newSurface.txt", allSurf);

	}


};

class Part {
public:
	Model_Solution m;
	PublicSolution ps;
	RWGeometric rwg;

	double m_r = 0.25;
	double m_s1 = 5;
	double m_s2 = 0.6;
	double m_s3 = 1;

	void maoPartOne() {
		Vec4 v1 = { 0,m_r,0,1 };
		Vec4 v2 = { 0,3 * m_s2+m_r,0,1 };
		Vec4 v3 = { m_s1,3 * m_s2 + m_r,0,1 };
		Vec4 v4 = { m_s1,m_r,0,1 };
		Vec4 v5 = { m_s3*2,3 * m_s2 + m_r,0,1 };
		Vec4 v6 = { m_s3 * 4,3 * m_s2 + m_r,0,1 };
		Vec4 v7 = { m_s3*3,m_r,0,1 };
		Vec4 v8 = { m_s3,m_r,0,1 };
		Vec4 v9 = { m_s3 * 3,3 * m_s2 + m_r,0,1 };

		Vec4 p1 = { m_s3,m_s2,0,1 };
		Vec4 p2 = { m_s3 * 2,2 * m_s2 + m_r+m_r,0,1 };
		Vec4 p3 = { m_s3 * 3,m_s2,0,1 };
		Vec4 p4 = { m_s3 * 4,2 * m_s2 + m_r+m_r,0,1 };
		Vec4 p5 = { m_s3 * 2,2 * m_s2,0,1 };
		Vec4 p6 = { m_s3+m_r,m_s2+m_r,0,1 };
		Vec4 p7 = { m_s3 * 3,m_s2 + m_r + m_r,0,1 };
		

		varray<Spline> outSpline;
		varray<Spline> subSpline;
		varray<varray<Spline>> inSpline;
		varray<Spline> addSpline;
		varray<SplineSurface> allSurf;
		varray<Spline> allSpline;

		Spline0 sl(v1, v2);
		outSpline.push_back(sl.getSpline());

		Spline0 sl1(v2, v5);
		outSpline.push_back(sl1.getSpline());
		
		Spline0 sl2(v5, v9);
		outSpline.push_back(sl2.getSpline());
		Spline0 SL5(v9, v6);
		outSpline.push_back(SL5.getSpline());

		Spline0 sl3(v6, v3);
		outSpline.push_back(sl3.getSpline());

		Spline0 SL(v3, v4);
		outSpline.push_back(SL.getSpline());

		Spline0 SL1(v4, v7);
		outSpline.push_back(SL1.getSpline());

		Spline0 SL2(v7, v8);
		outSpline.push_back(SL2.getSpline());

		Spline0 SL4(v8, v1);
		outSpline.push_back(SL4.getSpline());

		Spline0 sl4;
		subSpline = sl4.arcSplines(m_r);
		m.Trans(subSpline, m_s3, 1);
		m.Trans(subSpline, m_s2 + m_r, 2);
		inSpline.push_back(subSpline);

		m.Trans(subSpline, m_s3, 1);
		m.Trans(subSpline, m_s2, 2);
		inSpline.push_back(subSpline);

		m.Trans(subSpline, m_s3, 1);
		m.Trans(subSpline, m_s2, -2);
		inSpline.push_back(subSpline);

		m.Trans(subSpline, m_s3, 1);
		m.Trans(subSpline, m_s2, 2);
		inSpline.push_back(subSpline);

		Spline0 sl5(p1, v8);
		addSpline.push_back(sl5.getSpline());
		Spline0 sl6(p2, v5);
		addSpline.push_back(sl6.getSpline());
		Spline0 sl7(p3, v7);
		addSpline.push_back(sl7.getSpline());
		Spline0 sl8(p4, v6);
		addSpline.push_back(sl8.getSpline());
		Spline0 sl9(p5, p6);
		addSpline.push_back(sl9.getSpline());
		Spline0 sl10(p7, v9);
		addSpline.push_back(sl10.getSpline());

		for (auto&i : outSpline) {
			allSpline.push_back(i);
		}
		for (auto&i : inSpline) {
			for (auto&j : i) {
				allSpline.push_back(j);
			}
		}
		for (auto&i : addSpline) {
			allSpline.push_back(i);
		}

		varray<bool> genus;
		genus.resize(5);
		genus[0] = false;
		genus[1] = true;
		genus[2] = true;
		genus[3] = true;
		genus[4] = true;
		rwg.WriteSpline("E:\\kuang_models\\maoSplines.txt", allSpline);
		ps.quad(outSpline, inSpline, addSpline, genus, allSurf);
		//varray<SplineVolume> SV;
		//SV = m.CreatSweepVol(allSurf, 1, 3);
		ps.outPutTXT(allSurf, "E:\\kuang_models\\maoSurf");
		
		//rwg.WriteSplineSurface("E:\\kuang_models\\maoSurfaces.txt", allSurf);
		//rwg.WriteSplineVolume("E:\\kuang_models\\maoSPlineVolume.txt", SV);
		//ps.outPutVTK(SV, "E:\\kuang_models\\maoSPlineVolume.vtk");
	}
};





	