#include "QtGuiApplication1.h"
#include <QtWidgets/QApplication>
#include "MyDoc.h"
//#include <queue>
#include "FeatureNetwork.h"
#include "PolyIGA.h"
#include "quadPart.h"
#include "MeshQuality.h"

#include "kuang.h"
#include "PublicModels.h"




//两个静态成员的初始化
MyDoc::Ptr MyDoc::m_instance_ptr = nullptr;
std::mutex MyDoc::mutex;




void meshquality_test()
{
	RWGeometric rwg;
	//varray<SplineSurface> sf;
	//rwg.ReadSplineSurface("D:\\r\\sections.txt", sf);
	//SufQuality sq(sf);
	//sq.calQuality();
	varray<SplineVolume> sv;
	rwg.ReadSplineVolume("D:\\r\\fuld_vols.txt", sv);
	VolQuality vq(sv);
	vq.calQuality();
	
}

void fuld() {
	//RWGeometric rwg;
	//varray<SplineSurface> sf;
	//rwg.ReadSplineSurface("D:\\r\\sections.txt", sf);
	//varray<SplineVolume> sv;
	//Model_Solution M;
	//sv = M.CreatSweepVol(sf, 80, 3);
	//rwg.WriteSplineVolume("D:\\r\\fuld_vols.txt", sv);
	//varray<SplineVolume> targets;
	//varray<SplineVolume> res;
	//rwg.ReadSplineVolume("D:\\r\\rod_new2.txt", targets);
	//for (int i = 0; i < targets.size(); i+=6) {
	//	auto pt = targets[i].m_CtrlPts[6];
	//	auto tempV = sv;
	//	M.Trans(tempV, pt.x, 1);
	//	M.Trans(tempV, pt.y, 2);
	//	for (auto& s : tempV) {
	//		res.push_back(s);
	//	}
	//}
	//rwg.WriteSplineVolume("D:\\r\\rod_new3.txt", res);

	////52,53
	RWGeometric rwg;
	varray<SplineVolume> sv;
	rwg.ReadSplineVolume("D:\\r\\rod_new3.txt", sv);

	TestBolcks::pList plist;
	plist.getVolidx(sv);
	varray<int>WC;
	varray<int>WF;
	Model_Solution M;
	auto sf = M.GetSurfaces(sv);
	varray<SplineSurface> sufs;
	for (auto& s1 : sf) {
		for (auto& s2 : s1) {
			sufs.push_back(s2);
		}
	}

	plist.OutputParaVolumeDataTxt(sv, "D:\\r\\row_new3");
	cout << "WC ";
	for (int i = 52; i < sufs.size(); i+=54) {
		for (auto& ptwc : sufs[i].m_CtrlPts) {
			cout << plist.ptMap[ptwc] << " ";
		}
	}

	cout << "WF ";
	for (int i = 53; i < sufs.size(); i += 54) {
		for (auto& ptwc : sufs[i].m_CtrlPts) {
			cout << plist.ptMap[ptwc] << " ";
		}
	}



	
	
}



int main(int argc, char *argv[])
{
	/*Vec4 v1 = { 2,2,0,1 };
	Vec4 v2 = { 4,2,0,1 };
	Vec4 v3 = { 2,4,0,1 }; 
	Vec4 v4 = { 4,4,0,1 };

	varray<Spline> SL;
	Spline0 sl(v1,v2);
	SL.push_back(sl.getSpline());
	Spline0 sl1(v1, v3);
	SL.push_back(sl1.getSpline());
	Spline0 sl2(v3, v4);
	SL.push_back(sl2.getSpline());
	Spline0 sl3(v2, v4);
	SL.push_back(sl3.getSpline());
	varray<SplineSurface> SS;
	SplineSurface ss;
	ss.CoonsInterpolate(SL);
	SS.push_back(ss);
	PublicSolution ps;
	RWGeometric rwg;
	rwg.WriteSplineSurface("E:\\kuang_models\\asplinesurface.txt", SS);
	for (auto&i : SS) {
		varray<Spline> sl;
		i.GetEdgeLines(sl);
		rwg.WriteSpline("E:\\kuang_models\\asplines.txt", sl);
		ps.sortEdg(sl);
		rwg.WriteSpline("E:\\kuang_models\\asplines1.txt", sl);
	}*/
	


	Part pt;
	pt.maoPartOne();


	/*cout << "WC ";
	for (int i = 262; i < 304; i++) {
		cout << i * 6 + 4 << " ";
	}*/
	/*CarPart cp;
	cp.carPart();*/
	/*RWGeometric rwg;
	varray<SplineVolume> SV;
	rwg.ReadSplineVolume("E:\\kuang_models\\CarpartVolume.txt",SV);
	for (auto& i : SV) {
		i.Knots_Refine_Num(1);
	}
	rwg.WriteSplineVolume("E:\\kuang_models\\CarpartVolumeRef.txt", SV);
	PublicSolution ps;
	ps.setWCandWF("E:\\kuang_models\\CarpartVolumeRef.txt", "E:\\kuang_models\\CarpartVolume");*/
	
	
	

	QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
	QApplication a(argc, argv);
	QtGuiApplication1 w;
	w.show();
	return a.exec();
}
