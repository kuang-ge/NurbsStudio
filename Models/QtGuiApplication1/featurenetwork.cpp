#include "FeatureNetwork.h"

Feature_Line::Feature_Line(const Spline & nl, bool seg ):m_seg(seg)
{
	this->m_CtrlPts = nl.m_CtrlPts;
	this->m_Degree = nl.m_Degree;
	this->m_Knots = nl.m_Knots;
}

Feature_Line::Feature_Line(const Vec3 & p1, const Vec3 & p2, bool seg) :m_seg(seg)
{
	this->m_Degree = 1;
	this->m_Knots.push_back(0);
	this->m_Knots.push_back(0);
	this->m_Knots.push_back(1);
	this->m_Knots.push_back(1);
	this->m_CtrlPts.push_back(Vec4(p1.x,p1.y,p1.z, 1));
	this->m_CtrlPts.push_back(Vec4(p2.x, p2.y,p2.z, 1));
	this->DegreeElevate(2);
}

Spline Feature_Line::GetLine()
{
	return *this;
}



bool Model_Solution::JudgeTwoPointsCoincide(const Vec3 & p1, const Vec3 & p2)
{
	if ((abs(p1.x - p2.x) < 1e-4) && (abs(p1.y - p2.y) < 1e-4) && (abs(p1.z - p2.z) < 1e-4))
		return true;
	else
		return false;
}

void Model_Solution::MovePoint(Vec3 & p, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	if (mode == 1) {
		//Xƽ��
		p.x += dis;
	}
	else if (mode == 2)
		p.y += dis;
	else
		p.z += dis;
}

void Model_Solution::MovePoint(Vec4 & p, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	if (mode == 1) {
		//Xƽ��
		p.x += dis;
	}
	else if (mode == 2)
		p.y += dis;
	else
		p.z += dis;
}

void Model_Solution::MovePoints(varray<Vec3>& ps, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& p : ps) {
		if (mode == 1) {
			//Xƽ��
			p.x += dis;
		}
		else if (mode == 2)
			p.y += dis;
		else
			p.z += dis;
	}	
}

void Model_Solution::MovePoints(varray<Vec4>& ps, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& p : ps) {
		if (mode == 1) {
			//Xƽ��
			p.x += dis;
		}
		else if (mode == 2)
			p.y += dis;
		else
			p.z += dis;
	}
}

void Model_Solution::MoveLine(Spline & nl, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& p : nl.m_CtrlPts) {
		if (mode == 1) {
			//Xƽ��
			p.x += dis;
		}
		else if (mode == 2)
			p.y += dis;
		else
			p.z += dis;
	}
}

void Model_Solution::MoveLine(Feature_Line & fl, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& p : fl.m_CtrlPts) {
		if (mode == 1) {
			//Xƽ��
			p.x += dis;
		}
		else if (mode == 2)
			p.y += dis;
		else if (mode == 3)
			p.z += dis;
	}
}

void Model_Solution::MoveLine(Spline & nl, const Vec3 & p1, const Vec3 & p2)
{
	double dx, dy, dz;
	dx = p2.x - p1.x;
	dy = p2.y - p1.y;
	dz = p2.z - p1.z;
	
	if (abs(dx - 0.0) > 1e-4) {
		MoveLine(nl, dx, 1);
	}
	if (abs(dy - 0.0) > 1e-4) {
		MoveLine(nl, dy, 2);
	}
	if (abs(dz - 0.0) > 1e-4) {
		MoveLine(nl, dz, 3);
	}
}

void Model_Solution::MoveLine(Feature_Line & fl, const Vec3 & p1, const Vec3 & p2)
{
	double dx, dy, dz;
	dx = p2.x - p1.x;
	dy = p2.y - p1.y;
	dz = p2.z - p1.z;

	if (abs(dx - 0.0) > 1e-4) {
		MoveLine(fl, dx, 1);
	}
	if (abs(dy - 0.0) > 1e-4) {
		MoveLine(fl, dy, 2);
	}
	if (abs(dz - 0.0) > 1e-4) {
		MoveLine(fl, dz, 3);
	}
}

void Model_Solution::MoveLines(varray<Spline>& nls, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& l : nls) {
		MovePoints(l.m_CtrlPts, dis, mode);
	}
}

void Model_Solution::MoveLines(varray<Feature_Line>& fls, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& l : fls) {
		MovePoints(l.m_CtrlPts, dis, mode);
	}
}

void Model_Solution::MoveLines(varray<Spline>& nls, const Vec3 & p1, const Vec3 & p2)
{
	double dx, dy, dz;
	dx = p2.x - p1.x;
	dy = p2.y - p1.y;
	dz = p2.z - p1.z;

	if (abs(dx - 0.0) > 1e-4) {
		MoveLines(nls, dx, 1);
	}
	if (abs(dy - 0.0) > 1e-4) {
		MoveLines(nls, dy, 2);
	}
	if (abs(dz - 0.0) > 1e-4) {
		MoveLines(nls, dz, 3);
	}
}

void Model_Solution::MoveLines(varray<Feature_Line>& fls, const Vec3 & p1, const Vec3 & p2)
{
	double dx, dy, dz;
	dx = p2.x - p1.x;
	dy = p2.y - p1.y;
	dz = p2.z - p1.z;

	if (abs(dx - 0.0) > 1e-4) {
		MoveLines(fls, dx, 1);
	}
	if (abs(dy - 0.0) > 1e-4) {
		MoveLines(fls, dy, 2);
	}
	if (abs(dz - 0.0) > 1e-4) {
		MoveLines(fls, dz, 3);
	}
}

void Model_Solution::MoveSurface(SplineSurface & sf, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	MovePoints(sf.m_CtrlPts, dis, mode);
}

void Model_Solution::MoveSurfaces(varray<SplineSurface>& sfs, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& s : sfs) {
		MovePoints(s.m_CtrlPts, dis, mode);
	}
}

void Model_Solution::MoveSurface(SplineSurface & sf, const Vec3 & p1, const Vec3 & p2)
{
	double dx, dy, dz;
	dx = p2.x - p1.x;
	dy = p2.y - p1.y;
	dz = p2.z - p1.z;

	if (abs(dx - 0.0) > 1e-4) {
		MoveSurface(sf, dx, 1);
	}
	if (abs(dy - 0.0) > 1e-4) {
		MoveSurface(sf, dy, 2);
	}
	if (abs(dz - 0.0) > 1e-4) {
		MoveSurface(sf, dz, 3);
	}
}

void Model_Solution::MoveSurfaces(varray<SplineSurface>& sfs, const Vec3 & p1, const Vec3 & p2)
{
	double dx, dy, dz;
	dx = p2.x - p1.x;
	dy = p2.y - p1.y;
	dz = p2.z - p1.z;

	if (abs(dx - 0.0) > 1e-4) {
		MoveSurfaces(sfs, dx, 1);
	}
	if (abs(dy - 0.0) > 1e-4) {
		MoveSurfaces(sfs, dy, 2);
	}
	if (abs(dz - 0.0) > 1e-4) {
		MoveSurfaces(sfs, dz, 3);
	}
}


void Model_Solution::MoveSurface(Feature_Surface & fs, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	MovePoints(fs.m_edgePoints, dis, mode);
	MoveLines(fs.m_edges, dis, mode);
	if (fs.m_iscomplete) {
		MoveSurfaces(fs.m_ns, dis, mode);
	}	
}

void Model_Solution::MoveSurfaces(varray<Feature_Surface>& fss, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& s : fss) {
		MoveSurface(s, dis, mode);
	}
}

void Model_Solution::MoveSurface(Feature_Surface & fs, const Vec3 & p1, const Vec3 & p2)
{
	double dx, dy, dz;
	dx = p2.x - p1.x;
	dy = p2.y - p1.y;
	dz = p2.z - p1.z;

	if (abs(dx - 0.0) > 1e-4) {
		MoveSurface(fs, dx, 1);
	}
	if (abs(dy - 0.0) > 1e-4) {
		MoveSurface(fs, dy, 2);
	}
	if (abs(dz - 0.0) > 1e-4) {
		MoveSurface(fs, dz, 3);
	}
}

void Model_Solution::MoveSurfaces(varray<Feature_Surface>& fss, const Vec3& p1, const Vec3& p2)
{
	double dx, dy, dz;
	dx = p2.x - p1.x;
	dy = p2.y - p1.y;
	dz = p2.z - p1.z;

	if (abs(dx - 0.0) > 1e-4) {
		MoveSurfaces(fss, dx, 1);
	}
	if (abs(dy - 0.0) > 1e-4) {
		MoveSurfaces(fss, dy, 2);
	}
	if (abs(dz - 0.0) > 1e-4) {
		MoveSurfaces(fss, dz, 3);
	}
}

void Model_Solution::MoveVol(SplineVolume & v, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	MovePoints(v.m_CtrlPts, dis, mode);
}

void Model_Solution::MoveVols(varray<SplineVolume>& vs, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& v : vs) {
		MoveVol(v, dis, mode);
	}
}

void Model_Solution::MoveVol(Feature_Vol & fv, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	MoveSurfaces(fv.m_surf, dis, mode);
	MoveLines(fv.m_path, dis, mode);
	if (fv.Iscomplete()) {
		MoveVols(fv.m_fv, dis, mode);
	}
}

void Model_Solution::MoveVols(varray<Feature_Vol>& fvs, double dis, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& v : fvs) {
		MoveVol(v, dis, mode);
	}
}

void Model_Solution::MirrorPoint(const Vec3 & p1, Vec3 & p2, int mode)
{
	p2 = p1;
	if (mode == 1) {
		//YOZ����
		p2.x = -(p2.x);
	}
	if (mode == 2) {
		//XOZ����
		p2.y = -(p2.y);
	}
	if (mode == 3) {
		//XOY�᾵��
		p2.z = -(p2.z);
	}
}

void Model_Solution::MirrorPoint(const Vec4 & p1, Vec4 & p2, int mode)
{
	p2 = p1;
	if (mode == 1) {
		//YOZ�᾵��
		p2.x = -(p2.x);
	}
	if (mode == 2) {
		//XOZ�᾵��
		p2.y = -(p2.y);
	}
	if (mode == 3) {
		//XOY�᾵��
		p2.z = -(p2.z);
	}
}

void Model_Solution::MirrorPoints(const varray<Vec3>& p1, varray<Vec3>& p2, int mode)
{
	p2 = p1;
	for (auto& p : p2) {
		if (mode == 1) {
			//YOZ����
			p.x = -(p.x);
		}
		if (mode == 2) {
			//XOZ����
			p.y = -(p.y);
		}
		if (mode == 3) {
			//XOY����
			p.z = -(p.z);
		}
	}
}

void Model_Solution::MirrorPoints(const varray<Vec4>& p1, varray<Vec4>& p2, int mode)
{
	p2 = p1;
	for (auto& p : p2) {
		if (mode == 1) {
			//YOZ�᾵��
			p.x = -(p.x);
		}
		if (mode == 2) {
			//XOZ�᾵��
			p.y = -(p.y);
		}
		if (mode == 3) {
			//XOY�᾵��
			p.z = -(p.z);
		}
	}
}

void Model_Solution::MirrorLine(const Spline & l1, Spline & l2, int mode)
{
	l2 = l1;
	MirrorPoints(l1.m_CtrlPts, l2.m_CtrlPts, mode);
}

void Model_Solution::MirrorLines(const varray<Spline>& l1, varray<Spline>& l2, int mode)
{
	l2.clear();
	l2 = l1;
	for (int i = 0; i < l2.size(); ++i) {
		MirrorLine(l1[i], l2[i], mode);
	}
}

void Model_Solution::MirrorLine(const Feature_Line & fl1, Feature_Line & fl2, int mode)
{
	fl2 = fl1;
	MirrorPoints(fl1.m_CtrlPts, fl2.m_CtrlPts, mode);
}

void Model_Solution::MirrorLines(const varray<Feature_Line>& fl1, varray<Feature_Line>& fl2, int mode)
{
	fl2.clear();
	fl2 = fl1;
	for (int i = 0; i < fl2.size(); ++i) {
		MirrorLine(fl1[i], fl2[i], mode);
	}
}


void Model_Solution::MirrorSuface(const SplineSurface & s1, SplineSurface & s2, int mode)
{
	s2 = s1;
	MirrorPoints(s1.m_CtrlPts, s2.m_CtrlPts, mode);
}

void Model_Solution::MirrorSufaces(const varray<SplineSurface>& s1, varray<SplineSurface>& s2, int mode)
{
	s2.clear();
	s2.resize(s1.size());
	for (int i = 0; i < s2.size(); ++i) {
		MirrorSuface(s1[i], s2[i], mode);
	}
}

void Model_Solution::MirrorSuface(const Feature_Surface & fs1, Feature_Surface & fs2, int mode)
{
	fs2 = fs1;
	MirrorPoints(fs1.m_edgePoints, fs2.m_edgePoints, mode);
	MirrorLines(fs1.m_edges, fs2.m_edges, mode);
	MirrorSufaces(fs1.m_ns, fs2.m_ns, mode);
}

void Model_Solution::MirrorSufaces(const varray<Feature_Surface>& fs1, varray<Feature_Surface>& fs2, int mode)
{
	fs2.clear();
	fs2.resize(fs1.size());
	for (int i = 0; i < fs2.size(); ++i) {
		MirrorSuface(fs1[i],fs2[i], mode);
	}
}

void Model_Solution::MirrorVol(const SplineVolume & v1, SplineVolume & v2, int mode)
{
	v2 = v1;
	MirrorPoints(v1.m_CtrlPts, v2.m_CtrlPts, mode);
}

void Model_Solution::MirrorVols(const varray<SplineVolume>& v1, varray<SplineVolume>& v2, int mode)
{
	v2.clear();
	v2.resize(v1.size());
	for (int i = 0; i < v2.size(); ++i) {
		MirrorVol(v1[i], v2[i], mode);
	}
}

void Model_Solution::MirrorVol(const Feature_Vol & fv1, Feature_Vol & fv2, int mode)
{
	fv2 = fv1;
	MirrorLines(fv1.m_path, fv2.m_path, mode);
	MirrorSufaces(fv1.m_surf, fv2.m_surf, mode);
}

void Model_Solution::MirrorVols(const varray<Feature_Vol>& fv1, varray<Feature_Vol>& fv2, int mode)
{
	fv2.clear();
	fv2.resize(fv1.size());
	for (int i = 0; i < fv2.size(); ++i) {
		MirrorVol(fv1[i], fv2[i], mode);
	}
}

void Model_Solution::RotatePoints(varray<Vec3>& ps, double ang, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	if (mode == 1) {
		//��x����ת
		for (auto& p : ps) {
			p = p.RotateX(ang);
		}
	}
	else if (mode == 2) {
		//��Y����ת
		for (auto& p : ps) {
			p = p.RotateY(ang);
		}
	}
	else if (mode == 3) {
		//��Z����ת
		for (auto& p : ps) {
			p = p.RotateZ(ang);
		}
	}
}

void Model_Solution::RotatePoints(varray<Vec4>& ps, double ang, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	if (mode == 1) {
		//��x����ת
		for (auto& p : ps) {
			p = p.RotateX(ang);
		}
	}
	else if (mode == 2) {
		//��Y����ת
		for (auto& p : ps) {
			p = p.RotateY(ang);
		}
	}
	else if (mode == 3) {
		//��Z����ת
		for (auto& p : ps) {
			p = p.RotateZ(ang);
		}
	}
}

void Model_Solution::RotateLine(Spline& nl, double ang, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	RotatePoints(nl.m_CtrlPts, ang, mode);
}

void Model_Solution::RotateLine(Feature_Line & fl, double ang, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	RotatePoints(fl.m_CtrlPts, ang, mode);
}

void Model_Solution::RotateLines(varray<Spline>& nl, double ang, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& l : nl) {
		RotatePoints(l.m_CtrlPts, ang, mode);
	}
}

void Model_Solution::RotateLines(varray<Feature_Line>& fl, double ang, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& l : fl) {
		RotatePoints(l.m_CtrlPts, ang, mode);
	}
}

void Model_Solution::RotateSurface(SplineSurface & ns, double ang, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	RotatePoints(ns.m_CtrlPts, ang, mode);
}

void Model_Solution::RotateSurface(Feature_Surface & fs, double ang, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	RotatePoints(fs.m_edgePoints, ang, mode);
	RotateLines(fs.m_edges, ang, mode);
	if (fs.m_iscomplete) {
		RotateSurfaces(fs.m_ns, ang, mode);
	}
}

void Model_Solution::RotateSurfaces(varray<SplineSurface>& ns, double ang, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& s : ns) {
		RotatePoints(s.m_CtrlPts, ang, mode);
	}
}

void Model_Solution::RotateSurfaces(varray<Feature_Surface>& fs, double ang, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& s : fs) {
		RotateSurface(s, ang, mode);
	}
}

void Model_Solution::RotateVol(Feature_Vol & fv, double ang, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	RotateLines(fv.m_path, ang, mode);
	RotateSurfaces(fv.m_surf, ang, mode);
}

void Model_Solution::RotateVols(varray<Feature_Vol>& fvs, double ang, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3);
	for (auto& v : fvs) {
		RotateVol(v, ang, mode);
	}
}

void Model_Solution::RoArrayPoint(const Vec3 & p, varray<Vec3>& points, int num, int mode)
{
	double angle = PI * 2.0 / num;
	Vec3 tmp = p;
	points.clear();
	points.resize(num);
	points[0] = p;
	if (mode == 3) {
		//��Z������
		for (int i = 1; i < num; ++i) {
			tmp = tmp.RotateZ(angle);
			points[i] = tmp;
		}
	}
}

void Model_Solution::RoArrayPoint(const Vec4 & p, varray<Vec4>& points, int num, int mode)
{
	double angle = PI * 2.0 / num;
	Vec4 tmp = p;
	points.clear();
	points.resize(num);
	points[0] = p;
	if (mode == 3) {
		//��Z������
		for (int i = 1; i < num; ++i) {
			tmp = tmp.RotateZ(angle);
			points[i] = tmp;
		}
	}
}

void Model_Solution::RoArrayLine(const Spline & nl, varray<Spline>& nls, int num, int mode)
{
	double angle = PI * 2.0 / num;
	Spline tmp = nl;
	nls.clear();
	nls.resize(num);
	nls[0] = nl;
	if (mode == 3) {
		//��Z������
		for (int i = 1; i < num; ++i) {
			RotateLine(tmp, angle, mode);
			nls[i] = tmp;
		}
	}
}

void Model_Solution::RoArrayLine(const Feature_Line & fl, varray<Feature_Line>& fls, int num, int mode)
{
	double angle = PI * 2.0 / num;
	Feature_Line tmp = fl;
	fls.clear();
	fls.resize(num);
	fls[0] = fl;
	if (mode == 3) {
		//��Z������
		for (int i = 1; i < num; ++i) {
			RotateLine(tmp, angle, mode);
			fls[i] = tmp;
		}
	}
}

void Model_Solution::RoArrayLines(const varray<Spline>& nl, varray<varray<Spline>>& nls, int num, int mode)
{
	double angle = PI * 2.0 / num;
	varray<Spline> tmp = nl;
	nls.clear();
	nls.resize(num);
	nls[0] = nl;
	if (mode == 3) {
		//��Z������
		for (int i = 1; i < num; ++i) {
			for (auto& l : tmp) {
				RotateLine(l, angle, mode);
			}
			nls[i] = tmp;
		}
	}
}

void Model_Solution::RoArrayLines(const varray<Feature_Line>& fl, varray<varray<Feature_Line>>& fls, int num, int mode)
{
	double angle = PI * 2.0 / num;
	varray<Feature_Line> tmp = fl;
	fls.clear();
	fls.resize(num);
	fls[0] = fl;
	if (mode == 3) {
		//��Z������
		for (int i = 1; i < num; ++i) {
			for (auto& l : tmp) {
				RotateLine(l, angle, mode);
			}
			fls[i] = tmp;
		}
	}
}

void Model_Solution::RoArraySurface(const SplineSurface & sf, varray<SplineSurface>& sfs, int num, int mode)
{
	double angle = PI * 2.0 / num;
	SplineSurface tmp = sf;
	sfs.clear();
	sfs.resize(num);
	sfs[0] = sf;
	if (mode == 3) {
		//��Z������
		for (int i = 1; i < num; ++i) {
			RotateSurface(tmp, angle, mode);
			sfs[i] = tmp;
		}
	}
}

void Model_Solution::RoArraySurface(const Feature_Surface & fs, varray<Feature_Surface>& fss, int num, int mode)
{
	double angle = PI * 2.0 / num;
	Feature_Surface tmp = fs;
	fss.clear();
	fss.resize(num);
	fss[0] = fs;
	if (mode == 3) {
		//��Z������
		for (int i = 1; i < num; ++i) {
			RotateSurface(tmp, angle, mode);
			fss[i] = tmp;
		}
	}
}

void Model_Solution::RoArraySurfaces(const varray<SplineSurface>& sf, varray<varray<SplineSurface>>& sfs, int num, int mode)
{
	double angle = PI * 2.0 / num;
	varray<SplineSurface> tmp = sf;
	sfs.clear();
	sfs.resize(num);
	sfs[0] = sf;
	if (mode == 3) {
		//��Z������
		for (int i = 1; i < num; ++i) {
			for (auto& s : tmp) {
				RotateSurface(s, angle, mode);
			}
			sfs[i] = tmp;
		}
	}
}

void Model_Solution::RoArraySurfaces(const varray<Feature_Surface>& fs, varray<varray<Feature_Surface>>& fss, int num, int mode)
{
	double angle = PI * 2.0 / num;
	varray<Feature_Surface> tmp = fs;
	fss.clear();
	fss.resize(num);
	fss[0] = fs;
	if (mode == 3) {
		//��Z������
		for (int i = 1; i < num; ++i) {
			for (auto& s : tmp) {
				RotateSurface(s, angle, mode);
			}
			fss[i] = tmp;
		}
	}
}

void Model_Solution::RoArrayVol(const Feature_Vol & fv, varray<Feature_Vol>& fvs, int num, int mode)
{
	double angle = PI * 2.0 / num;
	Feature_Vol tmp = fv;
	fvs.clear();
	fvs.resize(num);
	fvs[0] = fv;
	if (mode == 3) {
		//��Z������
		for (int i = 1; i < num; ++i) {
			RotateVol(tmp, angle, mode);
			fvs[i] = tmp;
		}
	}
}

void Model_Solution::RoArrayVols(const varray<Feature_Vol>& fv, varray<varray<Feature_Vol>>& fvs, int num, int mode)
{
	double angle = PI * 2.0 / num;
	varray<Feature_Vol> tmp = fv;
	fvs.clear();
	fvs.resize(num);
	fvs[0] = fv;
	if (mode == 3) {
		//��Z������
		for (int i = 1; i < num; ++i) {
			for (auto& v : tmp) {
				RotateVol(v, angle, mode);
			}
			fvs[i] = tmp;
		}
	}
}

void Model_Solution::RoArrayVols(const varray<Feature_Vol>& fv, varray<Feature_Vol>& fvs, int num, int mode)
{
	varray<varray<Feature_Vol>> tmpfvs;
	RoArrayVols(fv, tmpfvs, num, mode);
	DimReduce(tmpfvs, fvs, true);
}

void Model_Solution::DimReduce(varray<varray<Feature_Vol>>& fvs1, varray<Feature_Vol>& fvs2, bool del)
{
	for (const auto& vs : fvs1) {
		for (const auto& v : vs) {
			fvs2.push_back(v);
		}
	}
	if (del)
		fvs1.clear();
}

varray<Feature_Line> Model_Solution::GetFlWithValue(varray<Feature_Line>& fls, double val, int mode)
{
	varray< Feature_Line> tmpfls;
	assert(mode == 1 || mode == 2 || mode == 3);
	if (mode == 1) {
		//xֵ����val
		for (const auto& l : fls) {
			bool flag = true;
			for (const auto& p : l.m_CtrlPts) {
				if (abs(p.x - val) >= 1e-4) {
					flag = false;
					break;
				}
			}
			if (flag)
				tmpfls.push_back(l);
		}
	}

	else if (mode == 2) {
		//yֵ����val
		for (const auto& l : fls) {
			bool flag = true;
			for (const auto& p : l.m_CtrlPts) {
				if (abs(p.y - val) >= 1e-4) {
					flag = false;
					break;
				}
			}
			if (flag)
				tmpfls.push_back(l);
		}
	}
	else if (mode == 3) {
		//zֵ����val
		for (const auto& l : fls) {
			bool flag = true;
			for (const auto& p : l.m_CtrlPts) {
				if (abs(p.z - val) >= 1e-4) {
					flag = false;
					break;
				}
			}
			if (flag)
				tmpfls.push_back(l);
		}
	}
	return tmpfls;
}

void Model_Solution::SortFeatureLines(varray<Feature_Line>& fls, int mode)
{
	assert(mode == 1 || mode == 2 || mode == 3 || mode == 4 || mode == 5 || mode == 6);
	Vec3 p1, p2;

	if (mode == 1) {
		//x��С����
		for (auto& l : fls) {
			p1 = *(l.m_CtrlPts.begin());
			p2 = *(l.m_CtrlPts.end() - 1);
			if (p1.x > p2.x)
				l.CruveReverse();
		}
		auto cmp = [=](const Feature_Line& l1, const Feature_Line& l2) {
			return l1.m_CtrlPts[0].x < l2.m_CtrlPts[0].x;
		};
		sort(fls.begin(), fls.end(), cmp);
	}
	if (mode == 2) {
		//y��С����
		for (auto& l : fls) {
			p1 = *(l.m_CtrlPts.begin());
			p2 = *(l.m_CtrlPts.end() - 1);
			if (p1.y > p2.y)
				l.CruveReverse();
		}
		auto cmp = [=](const Feature_Line& l1, const Feature_Line& l2) {
			return l1.m_CtrlPts[0].y < l2.m_CtrlPts[0].y;
		};
		sort(fls.begin(), fls.end(), cmp);
	}
	if (mode == 3) {
		//z��С����
		for (auto& l : fls) {
			p1 = *(l.m_CtrlPts.begin());
			p2 = *(l.m_CtrlPts.end() - 1);
			if (p1.z > p2.z)
				l.CruveReverse();
		}
		auto cmp = [=](const Feature_Line& l1, const Feature_Line& l2) {
			return l1.m_CtrlPts[0].z < l2.m_CtrlPts[0].z;
		};
		sort(fls.begin(), fls.end(), cmp);
	}

	if (mode == 4) {
		//x�Ӵ�С
		for (auto& l : fls) {
			p1 = *(l.m_CtrlPts.begin());
			p2 = *(l.m_CtrlPts.end() - 1);
			if (p1.x < p2.x)
				l.CruveReverse();
		}
		auto cmp = [=](const Feature_Line& l1, const Feature_Line& l2) {
			return l1.m_CtrlPts[0].x > l2.m_CtrlPts[0].x;
		};
		sort(fls.begin(), fls.end(), cmp);
	}
	if (mode == 5) {
		//y�Ӵ�С
		for (auto& l : fls) {
			p1 = *(l.m_CtrlPts.begin());
			p2 = *(l.m_CtrlPts.end() - 1);
			if (p1.y < p2.y)
				l.CruveReverse();
		}
		auto cmp = [=](const Feature_Line& l1, const Feature_Line& l2) {
			return l1.m_CtrlPts[0].y > l2.m_CtrlPts[0].y;
		};
		sort(fls.begin(), fls.end(), cmp);
	}
	if (mode == 6) {
		//z�Ӵ�С
		for (auto& l : fls) {
			p1 = *(l.m_CtrlPts.begin());
			p2 = *(l.m_CtrlPts.end() - 1);
			if (p1.z < p2.z)
				l.CruveReverse();
		}
		auto cmp = [=](const Feature_Line& l1, const Feature_Line& l2) {
			return l1.m_CtrlPts[0].z > l2.m_CtrlPts[0].z;
		};
		sort(fls.begin(), fls.end(), cmp);
	}
}

void Model_Solution::OrderCoonsLines(varray<Spline>& coonsline)
{
	varray<Spline> lines;
	Spline tempLine;
	Vec3 tempPoint;
	assert(coonsline.size() == 4);

	bool nextLine = false;
	map<int, bool> founded;
	for (int i = 0; i < 4; ++i) {
		founded[i] = false;
	}
	//��һ��
	tempLine = coonsline[0];
	tempPoint = tempLine.m_CtrlPts[0];
	lines.push_back(tempLine);
	founded[0] = true;
	//�ڶ���
	for (int i = 1; i < 4; ++i) {
		nextLine = false;
		if (founded[i]) continue;
		tempLine = coonsline[i];
		if (JudgeTwoPointsCoincide(tempLine.m_CtrlPts[0], tempPoint))
		{
			nextLine = true;
			founded[i] = true;
		}
		if (nextLine) break;
		if (JudgeTwoPointsCoincide(*(tempLine.m_CtrlPts.end() - 1), tempPoint)) {
			tempLine.CruveReverse();
			nextLine = true;
			founded[i] = true;
		}
		if (nextLine) break;
	}
	lines.push_back(tempLine);
	tempPoint = *(tempLine.m_CtrlPts.end() - 1);
	//������
	for (int i = 1; i < 4; ++i) {
		nextLine = false;
		if (founded[i]) continue;
		tempLine = coonsline[i];
		if (JudgeTwoPointsCoincide(tempLine.m_CtrlPts[0], tempPoint))
		{
			nextLine = true;
			founded[i] = true;
		}
		if (nextLine) break;
		if (JudgeTwoPointsCoincide(*(tempLine.m_CtrlPts.end() - 1), tempPoint)) {
			tempLine.CruveReverse();
			nextLine = true;
			founded[i] = true;
		}
		if (nextLine) break;
	}
	lines.push_back(tempLine);
	tempPoint = *(tempLine.m_CtrlPts.end() - 1);
	//������
	for (int i = 1; i < 4; ++i) {
		nextLine = false;
		if (founded[i]) continue;
		tempLine = coonsline[i];
		if (JudgeTwoPointsCoincide(tempLine.m_CtrlPts[0], tempPoint))
		{
			tempLine.CruveReverse();
			nextLine = true;
			founded[i] = true;
		}
		if (nextLine) break;
		if (JudgeTwoPointsCoincide(*(tempLine.m_CtrlPts.end() - 1), tempPoint)) {
			nextLine = true;
			founded[i] = true;
		}
		if (nextLine) break;
	}
	lines.push_back(tempLine);
	coonsline = lines;
}

varray<SplineSurface> Model_Solution::GetSurfsWithVols(const varray<SplineVolume>& vols, int dir)
{
	varray<SplineSurface> sfs;
	for (auto& v : vols) {
		SplineSurface s = v.GetSingleSurfaces(dir);
		sfs.push_back(s);
	}
	return sfs;
}

void Model_Solution::DecoincideNurbsLine(varray<Spline>& nl)
{
	assert(nl.size());
	for (int i = 0; i < nl.size(); ++i) {
		for (int j = i + 1; j < nl.size();) {
			Vec3 head1, head2, tail1, tail2;
			head1 = nl[i].m_CtrlPts[0];
			head2 = nl[j].m_CtrlPts[0];
			tail1 = *(nl[i].m_CtrlPts.end() - 1);
			tail2 = *(nl[j].m_CtrlPts.end() - 1);
			if (abs(head1.x - head2.x) < 1e-4 &&abs(head1.y - head2.y) < 1e-4&&abs(head1.z - head2.z) < 1e-4) {
				if (abs(tail1.x - tail2.x) < 1e-4 &&abs(tail1.y - tail2.y) < 1e-4&&abs(tail1.z - tail2.z) < 1e-4) {
					nl.erase(&nl[j]);
					continue;
				}
			}
			else if (abs(head1.x - tail2.x) < 1e-4 &&abs(head1.y - tail2.y) < 1e-4&&abs(head1.z - tail2.z) < 1e-4) {
				if (abs(tail1.x - head2.x) < 1e-4 &&abs(tail1.y - head2.y) < 1e-4&&abs(tail1.z - head2.z) < 1e-4) {
					nl.erase(&nl[j]);
					continue;
				}
			}
			++j;
		}
	}
}

void Model_Solution::SetQuadParameter(const varray<Feature_Surface>& fss, varray<varray<Spline>>& edgelines, varray<varray<int>>& seg, varray<bool>& genus, varray<int>& sfnum)
{
	edgelines.clear();
	seg.clear();
	genus.clear();
	sfnum.clear();
	varray<Spline> tmpline;
	varray<int> tmpseg;
	bool flag = false;
	if (fss[0].m_sfNum <= 0)
		flag = true;
	else {
		for (const auto& fs : fss) {
			//������������������(>=1)
			assert(fs.m_sfNum > 0);
		}
		flag = false;
	}

	for (auto& fs : fss) {
		tmpline.clear();
		tmpseg.clear();
		genus.push_back(fs.m_isempty);		//�����ʷֿ����Լ���
		//�����ʷ��������Ǽ���
		if (flag)
			sfnum.push_back(0);
		else
			sfnum.push_back(fs.m_sfNum);
		//�趨�ʷ�������߽����߼��ɷָ���
		for (const auto& l : fs.m_edges) {
			tmpline.push_back(l);
			tmpseg.push_back(l.m_seg);
		}
		edgelines.push_back(tmpline);
		seg.push_back(tmpseg);
	}
}

void Model_Solution::GetQuadPartData(varray<Feature_Surface>& fss, SfCtainTreeNode * root)
{
	assert(root);				//root��Ϊ��
	bool flag = false;
	int pos = -1;				//�������Ӧ��������λ��
	if (fss[0].m_sfNum <= 0)
		flag = true;			//�ʷּ����������������±�һһ��Ӧ(���1)
	else
		flag = false;			//�ʷּ��������������������һһ��Ӧ
	
	queue<SfCtainTreeNode*> q;
	varray<SplineSurface> sfs;
	varray <Spline> edgelines;
	varray<Vec3> edgeponts;
	q.push(root);

	while (!q.empty())
	{
		sfs.clear();
		edgelines.clear();
		edgeponts.clear();
		SfCtainTreeNode*cur = q.front();
		q.pop();

		//for (auto i : cur->quadPolNumber)
		//{	//��ȡ�ü����������������ݣ��ݲ���Ҫ
		//	for (auto j : i) {
		//		allLines.push_back(cur->allLines[j]);
		//	}
		//}

		//��ȡ����NURBS��������
		if (cur->isGenus) {
			//��������(����������)
			continue;
		}
		cur->GetSurfs(sfs);
		//��ȡ���б߽綥������
		cur->GetOutPoints(edgeponts);
		//��ȡ���б߽���������
		cur->GetOutLines(edgelines);

		if (flag) {
			pos = cur->num - 1;
		}
		else {
			for (int i = 0; i < fss.size(); ++i) {
				if (cur->num == fss[i].m_sfNum) {
					pos = i;
					break;
				}
					
			}
		}
		assert(pos >= 0);

		fss[pos].m_edgePoints = edgeponts;		//����߽綥������
		//����߽���������
		fss[pos].m_edges.resize(edgelines.size());
		for (int i = 0; i < edgelines.size(); ++i) {
			fss[pos].m_edges[i] = Feature_Line(edgelines[i], false);//���ɷָ�
		}
		//����NURBS��������
		fss[pos].m_ns = sfs;
		fss[pos].m_iscomplete = true;

		//��ȡ������������
		list<SfCtainTreeNode*>::iterator it = cur->childs.begin();
		for (; it != cur->childs.end(); it++) {
			q.push(*it);
		}
	}
	sfs.clear();
	edgelines.clear();
	edgeponts.clear();
}

void Model_Solution::QuadPartFeatureSurfaces(varray<Feature_Surface>& fss, varray<Spline>& conline)
{
	varray<varray<Spline>> edgeLines;	//������߽���
	varray<Spline> conl;					//������
	varray<varray<int>> seg;				//�߽��߿ɷָ���
	varray<bool> genus;						//������
	varray<int> sfNum;						//��������
	SfCtainTreeNode* root = nullptr;		//����������������

	//�趨�ʷֲ���
	SetQuadParameter(fss, edgeLines,seg, genus, sfNum);
	conl = conline;

	////�������start
	//for (int i = 0; i < edgeLines.size(); ++i) {
	//	cout << "��" << to_string(i+1) << "�������:" << endl;
	//	for (int j = 0; j < edgeLines[i].size(); j++) {
	//		Vec3 p1 = edgeLines[i][j].m_CtrlPts[0];
	//		Vec3 p2 = *(edgeLines[i][j].m_CtrlPts.end() - 1);
	//		cout << "line" << to_string(j + 1) << ":" << " head: (" << p1.x << "," << p1.y << "��" << p1.z << ")" <<
	//			"         tail: (" << p2.x << "," << p2.y << "��" << p2.z << ")" << endl;
	//	}
	//	cout << endl;
	//}
	////�������end

	//�������������
	root = CreateSurfContainTree(edgeLines, conl, seg, genus, sfNum);
	//ִ���ʷ�
	QuadWithContainTree(root);

	//��ȡ���ݴ��������棬���������湹��
	GetQuadPartData(fss, root);

	//�����Ƿ������������Ѿ�����
	for (const auto& fs : fss) {
		if (fs.m_isempty)
			continue;
		if (fs.m_ns.size() != 0) {
			assert(fs.m_iscomplete);
		}
		if (fs.m_iscomplete) {
			assert(fs.m_ns.size());
		}
	}

	delete root;
	root = nullptr;
}

void Model_Solution::QuadPartFeatureSurfaces(Feature_Surface & fs)
{
	varray<Feature_Surface> fss;
	varray<Spline> conLines;
	fss.push_back(fs);
	QuadPartFeatureSurfaces(fss, conLines);
	fs = fss[0];
}

Cir_arc::Cir_arc(double r, double ang, bool seg) :m_ang(ang)
{
	assert(ang < PI);
	Vec3 p(0, r, 0);
	Vec3 p0, p1, p2;
	double ang2 = 0.5*ang;
	double Op1 = r / cos(ang2);
	double w1 = cos(ang2);
	p0 = p.RotateZ(-ang2);
	p2 = p.RotateZ(ang2);
	p1 = Vec3(0, Op1, 0);

	m_arc.m_Degree = 2;
	m_arc.m_Knots.push_back(0);
	m_arc.m_Knots.push_back(0);
	m_arc.m_Knots.push_back(0);
	m_arc.m_Knots.push_back(1);
	m_arc.m_Knots.push_back(1);
	m_arc.m_Knots.push_back(1);
	m_arc.m_CtrlPts.push_back(Vec4(p0.x, p0.y,p0.z, 1));
	m_arc.m_CtrlPts.push_back(Vec4(p1.x,p1.y,p1.z, w1));
	m_arc.m_CtrlPts.push_back(Vec4(p2.x,p2.y,p2.z, 1));
	m_arc.m_seg = seg;
}

Feature_Line Cir_arc::Getarc()
{
	return m_arc;
}

Feature_Line Cir_arc::GetStdarc()
{
	Feature_Line stdArc = m_arc;
	double roAang = (PI - m_ang)*0.5;
	for (auto& p : stdArc.m_CtrlPts) {
		Vec3 tmp = p.RotateZ(-roAang);
		p = Vec4(tmp.x, tmp.y,tmp.z, p.w);
	}
	return stdArc;
}

Feature_Line Cir_arc::GetRoarc(double roang)
{
	Feature_Line roArc = this->GetStdarc();
	for (auto& p : roArc.m_CtrlPts) {
		Vec3 tmp = p.RotateZ(roang);
		p = Vec4(tmp.x,tmp.y,tmp.z, p.w);
	}
	return roArc;
}

Circle_2::Circle_2(double r, bool seg , bool empty , bool complete )
{
	//�����һ��Բ����
	Cir_arc* arc = new Cir_arc(r, PI / 2.0, seg);
	Feature_Line arcLine = arc->Getarc();
	delete arc;
	arc = nullptr;
	
	//���ñ߽綥��
	m_circle.m_edgePoints.resize(4);
	Vec3 p = arcLine.m_CtrlPts[0];
	m_circle.m_edgePoints[0] = p;
	p.y = -(p.y);
	m_circle.m_edgePoints[3] = p;
	p = *(arcLine.m_CtrlPts.end() - 1);
	m_circle.m_edgePoints[1] = p;
	p.y = -(p.y);
	m_circle.m_edgePoints[2] = p;

	//���ñ߽���
	varray<Feature_Line> edgelines;
	ms.RoArrayLine(arcLine, edgelines, 4);
	m_circle.m_edges = edgelines;
	edgelines.clear();

	//���ÿ���������ɶ�
	m_circle.m_isempty = empty;
	m_circle.m_iscomplete = complete;

	if (complete) {
		//ֱ��������湹��
		SplineSurface ns;
		varray<Spline> coonsLine;
		varray<SplineSurface> sfs;

		//�����ڲ�����
		Vec3 p0 = Vec3(0, r / 2.0, 0);
		p0 = p0.RotateZ(-PI / 4.0);
		Vec3 p1, p2, p3;
		ms.MirrorPoint(p0, p1, 1);
		ms.MirrorPoint(p1, p2, 2);
		ms.MirrorPoint(p0, p3, 2);
		Spline nl;
		nl.CreatLineWithTwoPoints(p2, p3);
		coonsLine.push_back(nl);
		nl.CreatLineWithTwoPoints(p2, p1);
		coonsLine.push_back(nl);
		nl.CreatLineWithTwoPoints(p1, p0);
		coonsLine.push_back(nl);
		nl.CreatLineWithTwoPoints(p3, p0);
		coonsLine.push_back(nl);
		ms.OrderCoonsLines(coonsLine);
		ns.CoonsInterpolate(coonsLine);
		m_circle.m_ns.push_back(ns);

		//����ʣ���ı�����
		p1 = m_circle.m_edgePoints[0];
		p2 = p0;
		p2.x = (-p2.x);
		p3= m_circle.m_edgePoints[1];
		coonsLine.clear();
		nl.CreatLineWithTwoPoints(p2, p0);
		coonsLine.push_back(nl);
		nl.CreatLineWithTwoPoints(p2, p3);
		coonsLine.push_back(nl);
		nl = arcLine;
		coonsLine.push_back(nl);
		nl.CreatLineWithTwoPoints(p0, p1);
		coonsLine.push_back(nl);
		ms.OrderCoonsLines(coonsLine);
		ns.CoonsInterpolate(coonsLine);
		ms.RoArraySurface(ns, sfs, 4);

		//�ϲ���������
		for (const auto& s : sfs) {
			m_circle.m_ns.push_back(s);
		}
	}

}

Feature_Surface Circle_2::GetCircle()
{
	return m_circle;
}

Ring::Ring(double r, double R, double angle)
{
	Cir_arc* arc1 = new Cir_arc(R, angle, false);
	Cir_arc* arc2 = new Cir_arc(r, angle, false);
	Feature_Line line1 = arc1->GetStdarc();
	Feature_Line line2 = arc2->GetStdarc();
	
	Vec3 p0, p1, p2, p3;
	p0 = line1.m_CtrlPts[0];
	p1 = *(line1.m_CtrlPts.end() - 1);
	p2 = *(line2.m_CtrlPts.end() - 1);
	p3 = line2.m_CtrlPts[0];

	//���ö���
	m_ring.m_edgePoints.push_back(p0);
	m_ring.m_edgePoints.push_back(p1);
	m_ring.m_edgePoints.push_back(p2);
	m_ring.m_edgePoints.push_back(p3);

	//���ñ߽���
	Feature_Line f1(p1, p2, false);
	Feature_Line f2(p3, p0, false);
	m_ring.m_edges.push_back(line1);
	m_ring.m_edges.push_back(f1);
	m_ring.m_edges.push_back(line2);
	m_ring.m_edges.push_back(f2);

	//������ɶȼ�����
	m_ring.m_iscomplete = true;
	m_ring.m_isempty = false;

	//����NURBS����
	varray<Spline> nls;
	SplineSurface sf;
	for (const auto& l : m_ring.m_edges) {
		nls.push_back(l);
	}
	ms.OrderCoonsLines(nls);
	sf.CoonsInterpolate(nls);
	m_ring.m_ns.push_back(sf);
}

Feature_Surface Ring::GetRing()
{
	return m_ring;
}

Rectangle_2::Rectangle_2(double length, double width, bool empty, bool complete)
{
	m_length = length;
	m_width = width;
	double len = length / 2.0;
	double wid = width / 2.0;
	Vec3 p0, p1, p2, p3;
	p0 = Vec3(len, wid, 0);
	ms.MirrorPoint(p0, p1, 1);
	ms.MirrorPoint(p1, p2, 2);
	ms.MirrorPoint(p0, p3, 2);

	//���ö���
	m_rectangle.m_edgePoints.push_back(p0);
	m_rectangle.m_edgePoints.push_back(p1);
	m_rectangle.m_edgePoints.push_back(p2);
	m_rectangle.m_edgePoints.push_back(p3);

	//���ñ߽�����
	Feature_Line* fl = new Feature_Line(p0, p1);
	m_rectangle.m_edges.push_back(*fl);
	delete fl;
	fl = nullptr;

	fl = new Feature_Line(p1, p2);
	m_rectangle.m_edges.push_back(*fl);
	delete fl;
	fl = nullptr;

	fl = new Feature_Line(p2, p3);
	m_rectangle.m_edges.push_back(*fl);
	delete fl;
	fl = nullptr;

	fl = new Feature_Line(p0, p3);
	m_rectangle.m_edges.push_back(*fl);
	delete fl;
	fl = nullptr;

	//���ÿ���������ɶ�
	m_rectangle.m_isempty = empty;
	m_rectangle.m_iscomplete = complete;

	//����NURBS����
	if (complete) {
		varray<Spline> nls;
		SplineSurface sf;
		for (auto& l : m_rectangle.m_edges) {
			l.m_seg = false;
			nls.push_back(l);
		}
		Spline tmpl = nls[0];
		nls[0] = nls[2];
		nls[2] = tmpl;
		ms.OrderCoonsLines(nls);
		sf.CoonsInterpolate(nls);
		m_rectangle.m_ns.push_back(sf);
	}
}

Feature_Surface Rectangle_2::GetRectangle()
{
	return m_rectangle;
}

Feature_Surface Rectangle_2::GetStdRectangle()
{
	Feature_Surface fs = m_rectangle;
	Vec3 p1(-(m_length)*0.5, -(m_width)*0.5, 0);
	Vec3 p2(0, 0, 0);
	ms.MoveSurface(fs, p1, p2);
	return fs;
}

Feature_Surface::Feature_Surface(const varray<Feature_Line>& fls, bool empty, int sfnum)
	:m_iscomplete(false), m_isempty(empty), m_sfNum(sfnum)
{
	m_edges = fls;
}

Feature_Surface::Feature_Surface(const varray<Feature_Line>& fls, const varray<Vec3>& ps, bool empty, int sfnum)
	:m_iscomplete(false), m_isempty(empty), m_sfNum(sfnum)
{
	m_edges = fls;
	m_edgePoints = ps;
}

Feature_Surface::Feature_Surface(const SplineSurface & s) :m_iscomplete(true), m_isempty(false)
{
	m_ns.push_back(s);
}

Feature_Surface::Feature_Surface(const varray<SplineSurface>& s) : m_iscomplete(true), m_isempty(false)
{
	m_ns = s;
}

varray<SplineSurface> Feature_Surface::GetSurfaces()const 
{
	assert(m_iscomplete);
	return m_ns;
}

Rec_Cir::Rec_Cir(double length, double width, double r)
{
	assert(2 * r < length && 2 * r < width);
	Feature_Surface* fs = new Feature_Surface();
	Rectangle_2* rec = new Rectangle_2(length, width, true, false);	//��������
	Cir_arc* arc = new Cir_arc(r, PI*0.5, false);					//������1/4Բ��
	*fs = rec->GetRectangle();
	m_rec_cir.m_edgePoints = fs->m_edgePoints;						//���ñ߽��
	m_rec_cir.m_edges = fs->m_edges;								//���ñ߽���

	//���ó�Ա����
	m_length = length;
	m_width = width;
	m_rec_cir.m_isempty = false;
	m_rec_cir.m_iscomplete = true;

	//��ȡԲ��������
	Feature_Line* fl = new Feature_Line();
	*fl = arc->Getarc();

	delete rec;
	rec = nullptr;
	delete arc;
	arc = nullptr;

	//����NURBS����
	m_rec_cir.m_ns.resize(4);
	Vec3 p0, p1, p2, p3;
	Spline nl;
	SplineSurface ns, ns2;
	varray<Spline> coonsLine;
	p0 = *(fl->m_CtrlPts.begin());
	p1 = *(fl->m_CtrlPts.end() - 1);
	p2 = fs->m_edgePoints[0];
	p3 = fs->m_edgePoints[1];

	nl = *fl;
	nl.CruveReverse();
	coonsLine.push_back(nl);					//��һ��coons����
	nl.CreatLineWithTwoPoints(p1, p3);
	coonsLine.push_back(nl);					//�ڶ���coons����
	coonsLine.push_back(fs->m_edges[0]);		//������coons����
	nl.CreatLineWithTwoPoints(p0, p2);
	coonsLine.push_back(nl);					//������coons����

	ms.OrderCoonsLines(coonsLine);				//coons��������
	ns.CoonsInterpolate(coonsLine);				//�����湹��
	coonsLine.clear();
	m_rec_cir.m_ns[0] = ns;

	ms.MirrorSuface(ns, ns2, 2);					//�����湹��
	m_rec_cir.m_ns[2] = ns2;

	//������������
	ms.RotateLine(*fl, PI / 2.0);
	p0 = *(fl->m_CtrlPts.begin());
	p1 = *(fl->m_CtrlPts.end() - 1);
	p2 = fs->m_edgePoints[1];
	p3 = fs->m_edgePoints[2];

	nl = *fl;
	nl.CruveReverse();
	coonsLine.push_back(nl);					//��һ��coons����
	nl.CreatLineWithTwoPoints(p1, p3);
	coonsLine.push_back(nl);					//�ڶ���coons����
	coonsLine.push_back(fs->m_edges[1]);		//������coons����
	nl.CreatLineWithTwoPoints(p0, p2);
	coonsLine.push_back(nl);					//������coons����

	ms.OrderCoonsLines(coonsLine);				//coons��������
	ns.CoonsInterpolate(coonsLine);				//�����湹��
	coonsLine.clear();
	m_rec_cir.m_ns[1] = ns;

	ms.MirrorSuface(ns, ns2, 1);				//�����湹��
	m_rec_cir.m_ns[3] = ns2;

	delete fl;
	fl = nullptr;
	delete fs;
	fs = nullptr;
}

Feature_Surface Rec_Cir::GetRecCir()
{
	return m_rec_cir;
}

Feature_Surface Rec_Cir::GetStdRecCir()
{
	Feature_Surface fs = m_rec_cir;
	Vec3 p1(-(m_length)*0.5, -(m_width)*0.5, 0);
	Vec3 p2(0, 0, 0);
	ms.MoveSurface(fs, p1, p2);
	return fs;
}

Mspring::Mspring(double D, double d, double t, double n)
{
	this->InitHighParameter(D, d, t, n);
	Feature_Surface fs;
	Feature_Line fl;
	Feature_Vol fv;

	//��������
	Circle_2* cir = new Circle_2(d*0.5);
	fs = cir->GetCircle();
	delete cir;
	cir = nullptr;
	ms.RotateSurface(fs, PI*0.5, 1);
	ms.MoveSurface(fs, D*0.5, 1);

	//����·������
	Cir_arc* arc = new Cir_arc(D*0.5, PI*0.5, false);
	fl = arc->GetStdarc();
	delete arc;
	arc = nullptr;
	assert(fl.m_CtrlPts.size() == 3);
	fl.m_CtrlPts[1].z += h_t / 8.0;
	fl.m_CtrlPts[2].z += h_t / 4.0;

	//����������
	Feature_Vol* tmpvol = new Feature_Vol(fs, fl, 2);

	int num = 4 * n;//����������
	for (int i = 0; i < num; ++i) {
		//�����幹��NURBS��
		tmpvol->ConstructVold(true);
		m_fvs.push_back(*tmpvol);

		//�����µ�·��
		ms.RotateLine(fl, PI*0.5);
		ms.MoveLine(fl, h_t / 4.0, 3);

		//�����µĽ���
		varray<SplineVolume> vols = tmpvol->GetVOls();
		varray<SplineSurface> sf = ms.GetSurfsWithVols(vols, 6);
		fs = Feature_Surface(sf);

		//�����µ�������
		delete tmpvol;
		tmpvol = nullptr;
		if (i != num - 1) {
			tmpvol = new Feature_Vol(fs, fl, 2);
		}
	}
}

void Mspring::InitHighParameter(double D, double d, double t, int n)
{
	this->h_D = D;
	this->h_d = d;
	this->h_t = t;
	this->h_n = n;
}

varray<SplineVolume> Mspring::GetMspring()
{
	varray<SplineVolume> vols;
	for (const auto& fv : m_fvs) {
		varray<SplineVolume> tmpv = fv.GetVOls();
		for (const auto& v : tmpv) {
			vols.push_back(v);
		}
	}
	return vols;
}

Feature_Vol::Feature_Vol(const Feature_Surface & fs, const Feature_Line & path, int mode, int dir, double dis)
	: m_iscomplete(false), m_mode(mode), m_dir(dir), m_dis(dis)
{
	m_surf.push_back(fs);
	m_path.push_back(path);
}

Feature_Vol::Feature_Vol(const varray<Feature_Surface>& fs, const Feature_Line & path, int mode, int dir, double dis)
	: m_iscomplete(false), m_mode(mode), m_dir(dir), m_dis(dis)
{
	m_surf = fs;
	m_path.push_back(path);
}

Feature_Vol::Feature_Vol(const Feature_Surface & fs, const varray<Feature_Line>& path, int mode, int dir, double dis)
	: m_iscomplete(false), m_mode(mode), m_dir(dir), m_dis(dis)
{
	m_surf.push_back(fs);
	m_path = path;
}

Feature_Vol::Feature_Vol(const varray<Feature_Surface>& fs, const varray<Feature_Line>& path, int mode, int dir, double dis)
	: m_iscomplete(false), m_mode(mode), m_dir(dir), m_dis(dis)
{
	m_surf = fs;
	m_path = path;
}

Feature_Vol::Feature_Vol(const Feature_Surface & fs, int mode, int dir, double dis)
	: m_iscomplete(false), m_mode(mode), m_dir(dir), m_dis(dis)
{
	m_surf.push_back(fs);
}

Feature_Vol::Feature_Vol(const varray<Feature_Surface>& fs, int mode, int dir, double dis)
	: m_iscomplete(false), m_mode(mode), m_dir(dir), m_dis(dis)
{
	m_surf = fs;
}

void Feature_Vol::ConstructVold(bool del)
{
	if (m_iscomplete) {
		//���Ѿ��������ģ�����
		if (del)
		{
			m_surf.clear();
			m_path.clear();
		}
		return;
	}

	assert(m_surf.size());//��������������Ϊ0
	if (m_mode == 2)
		assert(m_path.size());
	m_fv.clear();

	varray<SplineSurface> sfs;

	//����
	if (this->m_mode == 1) {
		assert(this->m_dir);
		assert(abs(this->m_dis - 0.0) > 1e-4);

		for (const auto& fs : m_surf)
		{
			if (!fs.m_isempty) {
				varray<SplineSurface> sf = fs.GetSurfaces();
				for (const auto& s : sf) {
					sfs.push_back(s);
				}
			}
		}
		
		//����·��
		if (this->m_dir == 1) {
			//X��������
			double x0 = sfs[0].m_CtrlPts[0].x;
			Feature_Line path(Vec3(x0, 0, 0), Vec3(x0 + m_dis, 0, 0), false);
			m_path.push_back(path);
		}
		else if (this->m_dir == 2) {
			//Y��������
			double y0 = sfs[0].m_CtrlPts[0].y;
			Feature_Line path(Vec3(0, y0, 0), Vec3(0, y0 + m_dis, 0), false);
			m_path.push_back(path);
		}
		else if (this->m_dir == 3) {
			//Z��������
			double z0 = sfs[0].m_CtrlPts[0].z;
			Feature_Line path(Vec3(0, 0, z0), Vec3(0, 0, z0 + m_dis), false);
			m_path.push_back(path);
		}

		//RWGeometric rw;
		//rw.WriteNurbsSurface("D:\\quadTest\\FeatureNetwork\\gear\\gear_surf.txt", sfs);

		int num_sf = sfs.size();//���������
		for (const auto& s : sfs) {
			SplineVolume v;
			v.CreateTransSweepSplineVolume(m_path[0], s);//���治��ֱ��·��
			m_fv.push_back(v);
		}


	}

	//��һ��ɨ��
	if (m_mode == 2) {
		//��ȡ���н���
		for (const auto& fs : m_surf)
		{
			if (!fs.m_isempty) {
				varray<SplineSurface> sf = fs.GetSurfaces();
				for (const auto& s : sf) {
					sfs.push_back(s);
				}
			}
		}

		int num_sf = sfs.size();//���������
		//��������NURBS��
		for (int i = 0; i < m_path.size(); ++i) {
			if (i != 0) {
				//��ȡ֮ǰɨ����UV����,W=1�Ľ���
				sfs.clear();
				for (int j = (i-1) * num_sf; j < i*num_sf; ++j) {
					SplineSurface s = m_fv[j].GetSingleSurfaces(6);
					sfs.push_back(s);
				}
			}
			for (const auto& s : sfs) {
				SplineVolume v;
				v.CreateTransSweepSplineVolume (m_path[i], s);
				m_fv.push_back(v);
			}
		}
	}

	//�ڶ���ɨ��
	if (m_mode == 3) {
		//��ȡ���н���
		for (const auto& fs : m_surf)
		{
			if (!fs.m_isempty) {
				varray<SplineSurface> sf = fs.GetSurfaces();
				for (const auto& s : sf) {
					sfs.push_back(s);
				}
			}
		}

		int num_sf = sfs.size();//���������
		//��������NURBS��
		for (int i = 0; i < m_path.size(); ++i) {
			if (i != 0) {
				//��ȡ֮ǰɨ����UV����,W=1�Ľ���
				sfs.clear();
				for (int j = (i-1) * num_sf; j < i*num_sf; ++j) {
					SplineSurface s = m_fv[j].GetSingleSurfaces(6);
					sfs.push_back(s);
				}
			}
			for (const auto& s : sfs) {
				SplineVolume v;
				v.CreateSweepSplineVolume(m_path[i], s, m_path[i].m_CtrlPts.size() - 1);
				m_fv.push_back(v);
			}
		}
	}

	if (m_mode == 4) {
		//����
		for (const auto& fs : m_surf) {
			for (const auto& s : fs.m_ns) {
				assert(int(s.m_CtrlPts.size() % 2));//�˴���ʱ�趨���Ƶ��������Ϊ����
			}
		}
		if (m_surf.size() == 2) {
			//ֻ����������
			Vec3 p1, p2;
			SplineSurface s1, s2;
			varray< SplineSurface> sfs;
			Spline path;
			SplineVolume v;
			s1 = m_surf[0].m_ns[0];
			s2 = m_surf[1].m_ns[0];
			sfs.push_back(s1);
			sfs.push_back(s2);
			p1 = s1.m_CtrlPts[s1.m_CtrlPts.size() / 2];
			p2 = s2.m_CtrlPts[s2.m_CtrlPts.size() / 2];
			//if (m_dir == 1) {
			//	//x�������
			//	p2.y = p1.y;
			//	p2.z = p1.z;
			//}
			//else if (m_dir == 2) {
			//	//y�������
			//	p2.x = p1.x;
			//	p2.z = p1.z;
			//}
			//else if (m_dir == 3) {
			//	//z�������
			//	p2.x = p1.x;
			//	p2.y = p1.y;
			//}
			path.CreatLineWithTwoPoints(p1, p2);
			v.LoftingSplineVolume(path, s1, s2);
			m_fv.push_back(v);
		}
		if (m_surf.size() > 2) {
			//������������
		}
	}

	m_iscomplete = true;//������NURBS�幹��
	if (del)
	{
		m_surf.clear();
		m_path.clear();
	}
}

bool Feature_Vol::Iscomplete()const
{
	return m_iscomplete;
}

varray<SplineVolume> Feature_Vol::GetVOls()const
{
	return m_fv;
}

Gear_Wheel::Gear_Wheel(double m, int z, double alph, double hax, double cx, double B, double x, double Dk)
	:m(m), z(z), alph(alph), hax(hax), cx(cx), B(B), x(x), Dk(Dk)
{
	RWGeometric rw;
	Gear_Straight* gs = new Gear_Straight(m, z, alph, hax, cx, B, x, Dk);
	varray<SplineVolume> spvol = gs->getGear();
	//rw.WriteSplineVolume("D:\\quadTest\\FeatureNetwork\\gear\\gear_all_��ʼ.txt", spvol);
	varray<Spline> edgeLines;			//��ݵı߽�����
	for (int i = 1; i < 6; i ++) {
		if (i == 1 || i == 2 || i == 4) {
			SplineVolume v = spvol[i];
			SplineSurface s = v.GetSingleSurfaces(1);
			varray<Spline> ls;
			s.GetEdgeLines(ls);
			for (auto& l : ls) {
				edgeLines.push_back(l);
			}
		}

	}
	spvol.clear();
	delete gs;
	gs = nullptr;

	//for (auto& l : edgeLines) {
	//	Vec3 p1 = l.m_CtrlPts[0];
	//	Vec3 p2 = *(l.m_CtrlPts.end() - 1);
	//	cout << "p1= (" << p1.x << "," << p1.y << "," << p1.z << ")" << "        "
	//		<< "p2= (" << p2.x << "," << p2.y << "," << p2.z << ")" << endl;
	//}
	//cout << endl << endl;

	//����ȥ��
	ms.DecoincideNurbsLine(edgeLines);

	//for (auto& l : edgeLines) {
	//	Vec3 p1 = l.m_CtrlPts[0];
	//	Vec3 p2 = *(l.m_CtrlPts.end() - 1);
	//	cout << "p1= (" << p1.x << "," << p1.y << "," << p1.z << ")" << "        "
	//		<< "p2= (" << p2.x << "," << p2.y << "," << p2.z << ")" << endl;
	//}
	//cout << endl << endl;

	//��ȡ�߽�����
	varray<int> flag(edgeLines.size(), 0);		//���߱��
	vector<Vec3> points1;					//����y�������ߵĶ���
	vector<Vec3> points2;					//ֻ��һ������y�������ߵĶ��������Сֵ����
	for (int i = 0; i < edgeLines.size();++i) {
		Vec3 p1, p2;
		p1 = *(edgeLines[i].m_CtrlPts.begin());
		p2 = *(edgeLines[i].m_CtrlPts.end() - 1);
		if (abs(p1.x - 0.0) < 1e-4 && abs(p2.x - 0.0) < 1e-4) {
			flag[i] = 1;
			bool flag2 = true;
			for (const auto& p : points1) {
				if (JudgeTwoPointsCoincide(p1, p)) {
					flag2 = false;
					break;
				}
			}
			if (flag2) points1.push_back(p1);
			flag2 = true;
			for (const auto& p : points1) {
				if (JudgeTwoPointsCoincide(p2, p)) {
					flag2 = false;
					break;
				}
			}
			if (flag2) points1.push_back(p2);
			continue;
		}
		else if (abs(p1.x - 0.0) < 1e-4) {
			flag[i] = 2;
			points2.push_back(p2);
			continue;
		}
		else if (abs(p2.x - 0.0) < 1e-4) {
			points2.push_back(p1);
			edgeLines[i].CruveReverse();
			flag[i] = 2;
			continue;
		}
	}

	assert(points1.size() == 3);
	assert(points2.size() == 4);
	//��y����С->������
	auto cmp = [=](const Vec3& p1, const Vec3& p2) {
		return p1.y < p2.y;
	};
	sort(points1.begin(), points1.end(), cmp);
	sort(points2.begin(), points2.end(), cmp);

	varray<Spline> tmpLines;
	for (int i = 0; i < edgeLines.size(); ++i) {
		if (flag[i] == 0) {
			tmpLines.push_back(edgeLines[i]);
		}
		else if (flag[i] == 2) {
			Vec3 p = *(edgeLines[i].m_CtrlPts.end() - 1);
			if (JudgeTwoPointsCoincide(p, points2[0]) || JudgeTwoPointsCoincide(p, *(points2.end() - 1)))
				tmpLines.push_back(edgeLines[i]);
		}
	}
	Spline l;
	l.CreatLineWithTwoPoints(points1[0], *(points1.end() - 1));
	tmpLines.push_back(l);
	edgeLines = tmpLines;
	tmpLines.clear();
	points1.clear();
	points2.clear();

	//�������������
	varray< Feature_Line> edge_fls;		//ĳ��������ı߽������߼���
	SfCtainTreeNode* root = nullptr;	//����������������
	varray<Vec3> edgepoints = OrderLinesAntioclock(edgeLines);
	//rw.WriteNurbsLine("D:\\quadTest\\FeatureNetwork\\gear\\��ݱ߽���.txt", edgeLines);
	
	for (const auto& l : edgeLines) {
		Feature_Line fl(l);
		edge_fls.push_back(fl);
	}
	edgeLines.clear();
	Feature_Surface fs(edge_fls, edgepoints);
	edgepoints.clear();

	//����ʷֲ���
	ms.QuadPartFeatureSurfaces(fs);

	//��������幹��
	Feature_Vol fv1(fs, 1, 3, B);
	//���������徵��
	Feature_Vol fv2;
	varray< Feature_Vol> fv_onegear;
	ms.MirrorVol(fv1, fv2, 1);
	fv_onegear.push_back(fv1);
	fv_onegear.push_back(fv2);

	//���л�ȡȫ��
	ms.RoArrayVols(fv_onegear, this->m_fvs, z);

	//����nurbs��
	for (auto& v : m_fvs) {
		v.ConstructVold();
	}
}

varray<SplineVolume> Gear_Wheel::GetGearVols()
{
	varray<SplineVolume> tmpvols;
	for (const auto& fv : m_fvs) {
		assert(fv.Iscomplete());
		varray<SplineVolume> vs = fv.GetVOls();
		for (const auto& v : vs) {
			tmpvols.push_back(v);
		}
	}
	return tmpvols;
}

Box_Quad::Box_Quad(double L, double L1, double L2, double L3, double W, double H, double H1, double SL1, double SL2, double r1, double r2, double r3, double r4, double t1, double t2, double tk1, double tk2, int n)
	:h_L(L), h_L1(L1), h_L2(L2), h_L3(L3), h_W(W), h_H(H), h_H1(H1), h_SL1(SL1), h_SL2(SL2), h_r1(r1), h_r2(r2), h_r3(r3), h_r4(r4), h_t1(t1), h_t2(t2), h_tk1(tk1), h_tk2(tk2), h_n(n)
{
	//ִ��ӳ��
	Map1();
	Map2();

	Vec3 p1, p2;
	Feature_Line* tmpFl=nullptr;				//��ʱ������
	Feature_Surface tmpFs, tmpFs2;				//��ʱ������
	Feature_Vol* tmpFv = nullptr;				//��ʱ������
	varray<Feature_Line> fls,tmpfls;			//�����߼���
	varray<Feature_Surface> fss, tmpfss;		//�����漯��
	varray<Feature_Vol> fvs,tmpfvs;				//�����弯��
	Spline tmpL;						
	varray<Spline> conls;					//�����߼���
	RWGeometric rw;								//���������

	//ǰ��3�������湹��
	//��������������
	Rectangle_2* rec = new Rectangle_2(m_L1, m_H1, false, false);
	tmpFs = rec->GetStdRectangle();
	fss.push_back(tmpFs);
	delete rec;
	rec = nullptr;

	//����С����
	rec= new Rectangle_2(m_SL1, m_SL1, true, false);
	tmpFs = rec->GetRectangle();
	ms.MoveSurface(tmpFs, m_L2, 1);
	ms.MoveSurface(tmpFs, m_H2, 2);
	fss.push_back(tmpFs);
	delete rec;
	rec = nullptr;

	//������Բ
	Circle_2* cir = new Circle_2(m_r1, false, true, false);
	tmpFs = cir->GetCircle();
	ms.MoveSurface(tmpFs, m_L4, 1);
	ms.MoveSurface(tmpFs, m_H2, 2);
	fss.push_back(tmpFs);
	delete cir;
	cir = nullptr;

	//����������
	p1 = fss[0].m_edgePoints[1];
	p2 = fss[1].m_edgePoints[1];
	tmpL.CreatLineWithTwoPoints(p1, p2);
	conls.push_back(tmpL);


	p1 = fss[0].m_edgePoints[2];
	p2 = fss[1].m_edgePoints[2];
	tmpL.CreatLineWithTwoPoints(p1, p2);
	conls.push_back(tmpL);


	p1 = fss[0].m_edgePoints[0];
	p2 = fss[2].m_edgePoints[0];
	tmpL.CreatLineWithTwoPoints(p1, p2);
	conls.push_back(tmpL);


	p1 = fss[0].m_edgePoints[3];
	p2 = fss[2].m_edgePoints[3];
	tmpL.CreatLineWithTwoPoints(p1, p2);
	conls.push_back(tmpL);


	//ִ���ʷ�
	ms.QuadPartFeatureSurfaces(fss, conls);
	conls.clear();

	//����Բ��
	Cir_Ring* m_ciring = new Cir_Ring(m_r3, m_r1);
	tmpFs = m_ciring->GetCirRing();
	ms.MoveSurface(tmpFs, m_L4, 1);
	ms.MoveSurface(tmpFs, m_H2, 2);
	fss.pop_back();
	fss.push_back(tmpFs);
	delete m_ciring;
	m_ciring = nullptr;

	////�������Start
	//for (int i = 0; i < fss.size(); ++i) {
	//	varray<Spline> ls;
	//	for (auto& l : fss[i].m_edges) {
	//		ls.push_back(l);
	//	}
	//	rw.WriteNurbsLine("D:\\quadTest\\FeatureNetwork\\Box_Quad\\l" + to_string(i + 1) + ".txt", ls);
	//	if (fss[i].m_ns.size()) {
	//		rw.WriteNurbsSurface("D:\\quadTest\\FeatureNetwork\\Box_Quad\\s" + to_string(i + 1) + ".txt", fss[i].m_ns);
	//	}
	//}

	//��תǰ��������90��
	ms.RotateSurfaces(fss, PI*0.5, 1);
	//����ǰ��������
	tmpFv = new Feature_Vol(fss, 1, 2, -m_tk1);
	m_fvs.push_back(*tmpFv);


	//����ǰСԲͲ������
	fss.clear();
	tmpFs = *(tmpFv->m_surf.end() - 1);
	delete tmpFv;
	tmpFv = nullptr;

	ms.MoveSurface(tmpFs, -m_tk1, 2);
	fss.push_back(tmpFs);
	for (int i = 1; i < h_n; ++i) {	
		ms.MoveSurface(tmpFs, -(m_t1 + m_t2), 2);
		fss.push_back(tmpFs);
	}
	tmpFv = new Feature_Vol(fss, 1, 2, -m_t2);
	m_fvs.push_back(*tmpFv);


	//����ǰ��ԲͲ������
	fss.clear();
	tmpFs = *(tmpFv->m_surf.begin());
	delete tmpFv;
	tmpFv = nullptr;

	ms.MoveSurface(tmpFs, -m_t2, 2);
	tmpfss.push_back(tmpFs);

	m_ciring = new Cir_Ring(m_r1, m_r2);
	tmpFs = m_ciring->GetCirRing();
	ms.MoveSurface(tmpFs, m_L4, 1);
	ms.MoveSurface(tmpFs, m_H2, 2);
	ms.RotateSurface(tmpFs, PI*0.5, 1);
	delete m_ciring;
	m_ciring = nullptr;

	ms.MoveSurface(tmpFs, -m_tk1 - m_t2, 2);
	tmpfss.push_back(tmpFs);
	for (auto& s : tmpfss) {
		fss.push_back(s);
	}
	for (int i = 1; i < h_n; ++i) {
		ms.MoveSurfaces(tmpfss, -(m_t1 + m_t2), 2);
		for (auto& s : tmpfss) {
			fss.push_back(s);
		}
	}
	tmpFv = new Feature_Vol(fss, 1, 2, -m_t1);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;
	tmpfss.clear();
	fss.clear();

	//������弰��ԲͲ������
	fvs = m_fvs;
	ms.MoveVols(fvs, -0.5*m_L1, 1);
	ms.MirrorVols(fvs, tmpfvs, 1);
	fvs = tmpfvs;
	tmpfvs.clear();
	ms.MoveVols(fvs, 0.5*m_L1, 1);
	ms.MoveVols(fvs, -0.5*m_W1, 2);
	ms.MirrorVols(fvs, tmpfvs, 2);
	ms.MoveVols(tmpfvs, 0.5*m_W1, 2);
	fvs = tmpfvs;
	tmpfvs.clear();
	for (auto&v : fvs) {
		//�޸����췽��
		v.m_dis = -v.m_dis;
		m_fvs.push_back(v);
	}
	fvs.clear();

	//�������Ұ�������
	fss.clear();
	rec = new Rectangle_2(m_W1, m_H1, false, true);
	tmpFs = rec->GetStdRectangle();
	delete rec;
	rec = nullptr;
	ms.RotateSurface(tmpFs, PI*0.5, 1);
	ms.RotateSurface(tmpFs, PI*0.5, 3);
	fss.push_back(tmpFs);
	ms.MoveSurface(tmpFs, m_L1 + m_tk1, 1);
	fss.push_back(tmpFs);
	tmpFv = new Feature_Vol(fss, 1, 1, -m_tk1);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;
	fss.clear();

	//�����װ�������
	tmpfls = ms.GetFlWithValue(m_fvs[0].m_surf[0].m_edges, 0.0, 3);
	ms.MoveLines(tmpfls, -m_L1 * 0.5, 1);
	ms.MirrorLines(tmpfls, fls, 1);
	ms.MoveLines(tmpfls, m_L1 * 0.5, 1); //�ָ���ԭλ
	ms.MoveLines(fls, m_L1 * 0.5, 1);

	ms.MoveLines(fls, -m_W1*0.5, 2);
	varray<Feature_Line> tmpbaselines;
	ms.MirrorLines(fls, tmpbaselines, 2);
	fls = tmpbaselines;
	tmpbaselines.clear();
	ms.MoveLines(fls, m_W1 * 0.5, 2);

	//��ӡ����
	for (const auto& l : tmpfls) {
		Vec3 p1, p2;
		p1 = *(l.m_CtrlPts.begin());
		p2 = *(l.m_CtrlPts.end()-1);
		cout << "head: (" << p1.x << "," << p1.y << "," << p1.z << ")          " << "tail: (" << p2.x << "," << p2.y << "," << p2.z << ")" << endl;
	}
	
	cout << endl << endl;
	for (const auto& l : fls) {
		Vec3 p1, p2;
		p1 = *(l.m_CtrlPts.begin());
		p2 = *(l.m_CtrlPts.end() - 1);
		cout << "head: (" << p1.x << "," << p1.y << "," << p1.z << ")          " << "tail: (" << p2.x << "," << p2.y << "," << p2.z << ")" << endl;
	}
	//


	for (const auto& l : tmpfls) {
		fls.push_back(l);
	}

	

	tmpfls.clear();
	tmpFl = new Feature_Line(l_p9, l_p10, false);
	fls.push_back(*tmpFl);
	ms.MoveLine(*tmpFl, m_L1, 1);
	fls.push_back(*tmpFl);
	delete tmpFl;
	tmpFl = nullptr;
	tmpFs = Feature_Surface(fls);
	fls.clear();
	ms.QuadPartFeatureSurfaces(tmpFs);
	tmpFv = new Feature_Vol(tmpFs, 1, 3, -m_tk2);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;

	//����x����ǰͻ����
	rec = new Rectangle_2(m_tk1, m_tk2, false, true);
	tmpFs = rec->GetStdRectangle();
	delete rec;
	rec = nullptr;
	ms.RotateSurface(tmpFs, PI*0.5, 1);
	ms.RotateSurface(tmpFs, -PI * 0.5, 3);
	ms.MoveSurface(tmpFs, -m_tk2, 3);
	tmpfls = ms.GetFlWithValue(m_fvs[0].m_surf[0].m_edges, 0.0, 3);
	ms.SortFeatureLines(tmpfls, 1);
	tmpFv = new Feature_Vol(tmpFs, tmpfls, 2);		//ģʽΪ3ע��
	m_fvs.push_back(*tmpFv);
	ms.MoveVol(*tmpFv, -m_tk1, 2);
	m_fvs.push_back(*tmpFv);
	ms.MoveVol(*tmpFv, m_tk1, 2);
	ms.MoveVol(*tmpFv, m_H1+m_tk2, 3);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;

	//����x�����ͻ����
	ms.MirrorLines(tmpfls, fls, 1);
	tmpfls.clear();
	ms.MoveLines(fls, m_L1, 1);
	ms.SortFeatureLines(fls, 1);
	ms.MoveSurface(tmpFs, m_W1 + m_tk1, 2);
	ms.MoveLines(fls, m_W1, 2);

	//
	for (const auto& l : fls) {
		Vec3 p1, p2;
		p1 = *(l.m_CtrlPts.begin());
		p2 = *(l.m_CtrlPts.end() - 1);
		cout << "head: (" << p1.x << "," << p1.y << "," << p1.z << ")          " << "tail: (" << p2.x << "," << p2.y << "," << p2.z << ")" << endl;
	}
	//

	tmpFv = new Feature_Vol(tmpFs, fls, 2);		//ģʽΪ3ע��
	fls.clear();
	m_fvs.push_back(*tmpFv);
	ms.MoveVol(*tmpFv, m_tk1, 2);
	m_fvs.push_back(*tmpFv);
	ms.MoveVol(*tmpFv, -m_tk1, 2);
	ms.MoveVol(*tmpFv, m_H1 + m_tk2, 3);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;

	//����y����ͻ����
	rec = new Rectangle_2(m_tk1, m_W1, false, true);
	tmpFs = rec->GetStdRectangle();
	delete rec;
	rec = nullptr;
	ms.MoveSurface(tmpFs, -m_tk1, 1);
	fss.push_back(tmpFs);
	ms.MoveSurface(tmpFs, -m_tk1, 1);
	fss.push_back(tmpFs);
	for (auto& s : fss) {
		tmpfss.push_back(s);
	}
	ms.MoveSurfaces(tmpfss, m_H1 + m_tk2, 3);
	for (auto& s : tmpfss) {
		fss.push_back(s);
	}
	ms.MoveSurfaces(fss, -m_L1*0.5, 1);
	ms.MirrorSufaces(fss, tmpfss, 1);
	for (auto&s : tmpfss) {
		fss.push_back(s);
	}
	ms.MoveSurfaces(fss, m_L1*0.5, 1);
	tmpfss.clear();
	tmpFv = new Feature_Vol(fss,1,3,-m_tk2);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;
	fss.clear();

	//����z����ͻ����
	rec = new Rectangle_2(m_tk1, m_tk1, false, true);
	tmpFs = rec->GetStdRectangle();
	delete rec;
	rec = nullptr;
	ms.MoveSurface(tmpFs, -m_tk1, 1);
	ms.MoveSurface(tmpFs, -m_tk1, 2);
	tmpfss.push_back(tmpFs);
	ms.MoveSurface(tmpFs, m_tk1 + m_W1, 2);
	tmpfss.push_back(tmpFs);

	fss = tmpfss;
	for (auto& s: tmpfss) {
		ms.MoveSurface(s, m_tk1 + m_L1, 1);
		fss.push_back(s);
	}
	tmpfss.clear();

	tmpFv = new Feature_Vol(fss, 1, 3, m_H1);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;
	fss.clear();

	//�����Ľ�start
	rec = new Rectangle_2(m_tk1, m_tk1, false, true);
	tmpFs = rec->GetStdRectangle();
	delete rec;
	rec = nullptr;
	ms.MoveSurface(tmpFs, -m_tk1, 1);
	ms.MoveSurface(tmpFs, -m_tk1, 2);
	fss.push_back(tmpFs);
	tmpFs2 = tmpFs;
	ms.MoveSurface(tmpFs2, m_H1 + m_tk2, 3);
	fss.push_back(tmpFs2);
	tmpFs2 = tmpFs;
	ms.MoveSurface(tmpFs2, -m_SL2, 1);
	fss.push_back(tmpFs2);
	tmpFs2 = tmpFs;
	ms.MoveSurface(tmpFs2, -m_SL2, 2);
	fss.push_back(tmpFs2);

	Rec_Cir* rc = new Rec_Cir(m_SL2, m_SL2, m_r4);
	tmpFs = rc->GetStdRecCir();
	delete rc;
	rc = nullptr;
	ms.MoveSurface(tmpFs, -m_tk1 - m_SL2, 1);
	ms.MoveSurface(tmpFs, -m_tk1 - m_SL2, 2);
	tmpFs2 = tmpFs;
	ms.MoveSurface(tmpFs2, m_SL2, 2);
	ms.MoveSurface(tmpFs2, m_H1 + m_tk2, 3);
	fss.push_back(tmpFs);
	fss.push_back(tmpFs2);

	//�������ĸ���
	ms.MoveSurfaces(fss, -m_W1 * 0.5, 2);
	ms.MirrorSufaces(fss, tmpfss, 2);
	for (auto& s : tmpfss) {
		fss.push_back(s);
	}
	tmpfss.clear();
	ms.MoveSurfaces(fss, m_W1 * 0.5, 2);

	ms.MoveSurfaces(fss, -m_L1 * 0.5, 1);
	ms.MirrorSufaces(fss, tmpfss, 1);
	for (auto& s : tmpfss) {
		fss.push_back(s);
	}
	tmpfss.clear();
	ms.MoveSurfaces(fss, m_L1 * 0.5, 1);

	tmpFv = new Feature_Vol(fss, 1, 3, -m_tk2);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;
	fss.clear();
	//�����Ľ�end

	//����nurbs��
	for (auto& v : m_fvs) {
		v.ConstructVold();
	}

}

void Box_Quad::Map1()
{
	m_L = h_L;
	m_L1 = h_L1;
	m_L2 = h_L2;
	m_L3 = h_L3;
	m_W = h_W;
	m_H = h_H;
	m_H1 = h_H1;
	m_SL1 = h_SL1;
	m_SL2 = h_SL2;
	m_r1 = h_r1;
	m_r2 = h_r2;
	m_r3 = h_r3;
	m_r4 = h_r4;
	m_t1 = h_t1;
	m_t2 = h_t2;
	m_tk1 = h_tk1;
	m_tk2 = h_tk2;

	m_L4 = m_L1 - m_L3;
	m_H2 = 0.5*m_H1;
	m_W1 = m_W - 2 * m_SL2 - 2 * m_tk1;
}

void Box_Quad::Map2()
{
	Vec3 p(0, 0, 0);
	varray<Vec3> ps1, ps2;
	l_p9 = p;
	l_p10 = l_p9;
	ms.MovePoint(l_p10, m_W1, 2);

	l_p1 = p;
	ms.MovePoint(l_p1, -m_tk1 - 0.5*m_SL2, 1);
	ms.MovePoint(l_p1, -0.5*m_SL2, 2);
	l_p2 = l_p1;
	ms.MovePoint(l_p2, m_tk1, 2);
	ms.MovePoint(l_p2, m_tk2 + m_H1, 3);
	ps1.clear();
	ps1.push_back(l_p1);
	ps1.push_back(l_p2);
	ms.MovePoints(ps1, -m_W1 * 0.5, 2);
	ms.MirrorPoints(ps1, ps2, 2);
	ms.MovePoints(ps2, m_W1 * 0.5, 2);
	l_p3 = ps2[0];
	l_p4 = ps2[1];
	ps1.clear();
	ps2.clear();
	ps1.resize(4);
	ps1[0] = l_p1;
	ps1[1] = l_p2;
	ps1[2] = l_p3;
	ps1[3] = l_p4;
	ms.MovePoints(ps1, -m_L1 * 0.5, 1);
	ms.MirrorPoints(ps1, ps2, 2);
	l_p5 = ps2[0];
	l_p6 = ps2[1];
	l_p7 = ps2[2];
	l_p8 = ps2[3];
	ps1.clear();
	ps2.clear();
}

varray<SplineVolume> Box_Quad::GetBox()
{
	varray<SplineVolume> tmpvols;
	for (const auto& fv : m_fvs) {
		assert(fv.Iscomplete());
		varray<SplineVolume> vs = fv.GetVOls();
		for (const auto& v : vs) {
			tmpvols.push_back(v);
		}
	}
	return tmpvols;
}

Cir_Ring::Cir_Ring(double r, double R)
{
	varray<Feature_Surface> sfs;//�����Ķ�Բ��
	Ring* m_ring = new Ring(r, R);
	Feature_Surface surf_ring = m_ring->GetRing();
	delete m_ring;
	m_ring = nullptr;
	ms.RotateSurface(surf_ring, PI*0.25);
	ms.RoArraySurface(surf_ring, sfs, 4);
	for (const auto& s : sfs) {
		m_cir_ring.m_edgePoints.push_back(s.m_edgePoints[0]);
		m_cir_ring.m_edges.push_back(s.m_edges[0]);
		for (const auto& ns : s.m_ns) {
			m_cir_ring.m_ns.push_back(ns);
		}
	}
	m_cir_ring.m_iscomplete = true;
	m_cir_ring.m_isempty = false;
}

Feature_Surface Cir_Ring::GetCirRing()
{
	return m_cir_ring;
}

Bearing_Chock::Bearing_Chock(double L, double L1, double L2, double L3, double W, double W1, double W2, double W3, double W4, double W5, double H1, double H2, double d1, double r1, double r2, double sita)
	:h_L(L), h_L1(L1), h_L2(L2), h_L3(L3), h_W(W), h_W1(W1), h_W2(W2), h_W3(W3), h_W4(W4), h_W5(W5), h_H1(H1), h_H2(H2), h_d1(d1), h_r1(r1), h_r2(r2), h_sita(sita)
{
	//ִ��ӳ��
	Map1();
	Map2();

	Vec3 p1, p2;
	Feature_Line* tmpFl = nullptr;				//��ʱ������
	Feature_Surface tmpFs, tmpFs2;				//��ʱ������
	Feature_Vol* tmpFv = nullptr;				//��ʱ������
	varray<Feature_Line> fls, tmpfls;			//�����߼���
	varray<Feature_Surface> fss, tmpfss;		//�����漯��
	varray<Feature_Vol> fvs, tmpfvs;			//�����弯��
	Spline tmpL;
	varray<Spline> conls, tmpnls;			//�����߼���
	SplineSurface tmpns;
	RWGeometric rw;								//���������

	//�����װ�����������(0)
	Rectangle_2* rec = new Rectangle_2(m_L, m_W, false, false);
	tmpFs = rec->GetRectangle();
	ms.MoveSurface(tmpFs, 0.5*m_W, 2);
	fss.push_back(tmpFs);
	delete rec;
	rec = nullptr;

	//�����м����(��)(1)
	rec = new Rectangle_2(m_L3, m_W2, false, false);
	tmpFs = rec->GetRectangle();
	for (auto& l : tmpFs.m_edges) {
		l.m_seg = false;
	}
	ms.MoveSurface(tmpFs,Vec3(0,0,0),l_p1);
	fss.push_back(tmpFs);
	delete rec;
	rec = nullptr;

	//�����м����(��)(2)
	rec = new Rectangle_2(m_L3,m_W3, false, false);
	tmpFs = rec->GetRectangle();
	for (auto& l : tmpFs.m_edges) {
		l.m_seg = false;
	}
	ms.MoveSurface(tmpFs, Vec3(0, 0, 0), l_p2);
	fss.push_back(tmpFs);
	delete rec;
	rec = nullptr;

	//���������(3)
	rec = new Rectangle_2(m_L5, m_W3, false, false);
	tmpFs = rec->GetRectangle();
	for (auto& l : tmpFs.m_edges) {
		l.m_seg = false;
	}
	ms.MoveSurface(tmpFs, Vec3(0, 0, 0), l_p3);
	fss.push_back(tmpFs);
	delete rec;
	rec = nullptr;

	//�����Ҿ���(����)(4)
	ms.MirrorSuface(tmpFs, tmpFs2, 1);
	fss.push_back(tmpFs2);

	//��������Բ(5,6)
	Circle_2* cir = new Circle_2(0.5*m_d1, false, true, false);
	tmpFs = cir->GetCircle();
	ms.RotateSurface(tmpFs, h_sita);
	ms.MoveSurface(tmpFs, Vec3(0, 0, 0), l_p4);
	ms.MirrorSuface(tmpFs, tmpFs2, 1);
	fss.push_back(tmpFs);
	fss.push_back(tmpFs2);
	delete cir;
	cir = nullptr;

	//����������

	p1 = fss[0].m_edgePoints[1];
	p2 = fss[5].m_edgePoints[1];
	tmpL.CreatLineWithTwoPoints(p1, p2);
	conls.push_back(tmpL);

	p1 = fss[3].m_edgePoints[2];
	p2 = fss[5].m_edgePoints[0];
	tmpL.CreatLineWithTwoPoints(p1, p2);
	conls.push_back(tmpL);

	//p1 = fss[0].m_edgePoints[2];
	//p2 = fss[5].m_edgePoints[2];
	//tmpL.CreatLineWithTwoPoints(p1, p2);
	//conls.push_back(tmpL);

	ms.MirrorLines(conls, tmpnls, 1);
	for (auto& l : tmpnls) {
		conls.push_back(l);
	}
	tmpnls.clear();

	//ִ���ʷ�
	ms.QuadPartFeatureSurfaces(fss, conls);
	conls.clear();

	//����װ�������
	tmpFv = new Feature_Vol(fss, 1, 3, -m_H2);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;

	//����ײ�Բ����(0)
	fss.clear();
	Ring* ring = new Ring(m_r1, m_r2, m_sita_1);
	tmpFs = ring->GetRing();
	ms.RotateSurface(tmpFs, 0.5*(PI - m_sita_1));
	ms.RotateSurface(tmpFs, PI);
	fss.push_back(tmpFs);
	delete ring;
	ring = nullptr;

	//��������Բ����(1,2)
	ring = new Ring(m_r1, m_r2, m_sita_2);
	tmpFs = ring->GetRing();
	ms.RotateSurface(tmpFs, 0.5*(PI - m_sita_2));
	ms.RotateSurface(tmpFs, PI - 0.5*(m_sita_1 + m_sita_2));
	fss.push_back(tmpFs);
	delete ring;
	ring = nullptr;

	ms.MirrorSuface(tmpFs, tmpFs2, 1);
	fss.push_back(tmpFs);
	fss.push_back(tmpFs2);

	//���춥��Բ����(3)
	ring = new Ring(m_r1, m_r2, m_sita_3);
	tmpFs = ring->GetRing();
	ms.RotateSurface(tmpFs, 0.5*(PI - m_sita_3));
	fss.push_back(tmpFs);
	delete ring;
	ring = nullptr;

	ms.RotateSurfaces(fss, 0.5*PI, 1);
	ms.MoveSurfaces(fss, Vec3(0, 0, 0), l_p5);

	//����ԲͲ������
	tmpFv= tmpFv = new Feature_Vol(fss, 1, 2, m_W6);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;

	tmpFv = tmpFv = new Feature_Vol(fss, 1, 2, -m_W3);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;

	ms.MoveSurfaces(fss, -m_W3, 2);
	tmpFv = tmpFv = new Feature_Vol(fss, 1, 2, -m_W5);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;
	ms.MoveSurfaces(fss, m_W3, 2);

	//�����м�֧�Ű�
	tmpfss.clear();
	tmpnls.clear();
	tmpnls.push_back(fss[0].m_edges[0]);
	p1 = fss[0].m_edgePoints[0];
	p2 = p1;
	p2.z = 0;
	tmpL.CreatLineWithTwoPoints(p1, p2);
	tmpnls.push_back(tmpL);

	p1 = p2;
	p1.x += m_L3;
	tmpL.CreatLineWithTwoPoints(p1, p2);
	tmpnls.push_back(tmpL);

	p2 = fss[0].m_edgePoints[1];
	tmpL.CreatLineWithTwoPoints(p1, p2);
	tmpnls.push_back(tmpL);

	ms.OrderCoonsLines(tmpnls);
	tmpns.CoonsInterpolate(tmpnls);

	tmpFs = Feature_Surface(tmpns);
	tmpfss.push_back(tmpFs);

	//����֧�Ű�
	tmpnls.clear();
	tmpnls.push_back(fss[1].m_edges[0]);
	p1 = fss[1].m_edgePoints[0];
	p2 = Vec3(0, 0, 0);
	ms.MovePoint(p2, m_W, 2);
	ms.MovePoint(p2, -0.5*m_L2, 1);
	tmpL.CreatLineWithTwoPoints(p1, p2);
	tmpnls.push_back(tmpL);

	p1 = p2;
	ms.MovePoint(p1, m_L5, 1);
	tmpL.CreatLineWithTwoPoints(p1, p2);
	tmpnls.push_back(tmpL);

	p2= fss[1].m_edgePoints[1];
	tmpL.CreatLineWithTwoPoints(p1, p2);
	tmpnls.push_back(tmpL);

	ms.OrderCoonsLines(tmpnls);
	tmpns.CoonsInterpolate(tmpnls);

	tmpFs = Feature_Surface(tmpns);
	ms.MirrorSuface(tmpFs, tmpFs2, 1);
	tmpfss.push_back(tmpFs);
	tmpfss.push_back(tmpFs2);
	tmpnls.clear();

	//����֧�Ű�������
	tmpFv = tmpFv = new Feature_Vol(tmpfss, 1, 2, -m_W3);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;
	tmpfss.clear();
	fss.clear();

	//�����߰�
	varray<SplineSurface> nss;//���������
	fss.clear();
	tmpnls.clear();
	tmpFs = (m_fvs.end()-2)->m_surf[0];
	tmpFv = new Feature_Vol(tmpFs, 1, 2, -m_W5);
	tmpFv->ConstructVold();
	varray<SplineVolume> vol = tmpFv->GetVOls();
	//for (int i = 0; i < vol[0].m_CtrlPts.size(); i += 3) {
	//	
	//	for (int j = i; j < i + 3; ++j) {
	//		Vec4 p = vol[0].m_CtrlPts[j];
	//		cout << "��" << to_string(j) << ": (" << p.x << "," << p.y << "," << p.z << "," << p.w<< ")       ";
	//	}
	//	cout << endl;
	//}

	tmpns = vol[0].GetSingleSurfaces(5);
	nss.push_back(tmpns);//�������
	tmpFs = Feature_Surface(tmpns);
	fss.push_back(tmpFs);
	delete tmpFv;
	tmpFv = nullptr;

	rec = new Rectangle_2(m_L3, m_W2, false, false);
	tmpFs = rec->GetRectangle();
	for (auto& l : tmpFs.m_edges) {
		l.m_seg = false;
	}
	ms.MoveSurface(tmpFs, Vec3(0, 0, 0), l_p1);
	delete rec;
	rec = nullptr;
	for (auto& l : tmpFs.m_edges) {
		tmpnls.push_back(l);
	}
	//tmpnls[0].CruveReverse();
	ms.OrderCoonsLines(tmpnls);
	tmpns.CoonsInterpolate(tmpnls);
	tmpns.OrderCtrlPts();
	nss.push_back(tmpns);//�������
	tmpFs = Feature_Surface(tmpns);
	fss.push_back(tmpFs);
	//rw.WriteSplineSurface("D:\\quadTest\\FeatureNetwork\\bearing\\fangyang.txt", nss);

	tmpFv = new Feature_Vol(fss, 4, 3);
	m_fvs.push_back(*tmpFv);
	delete tmpFv;
	tmpFv = nullptr;

	//����nurbs��
	for (auto& v : m_fvs) {
		v.ConstructVold();
	}
}

void Bearing_Chock::Map1()
{
	m_L = h_L;
	m_L1 = h_L1;
	m_L2 = h_L2;
	m_L3 = h_L3;
	m_W = h_W;
	m_W1 = h_W1;
	m_W2 = h_W2;
	m_W3 = h_W3;
	m_W4 = h_W4;
	m_W5 = h_W5;
	m_H1 = h_H1;
	m_H2 = h_H2;
	m_d1 = h_d1;
	m_r1 = h_r1;
	m_r2 = h_r2;

	m_L4 = 0.5*m_L - m_L1;
	m_L5 = 0.5*(m_L2 - m_L3);
	m_W6 = m_W4 - m_W5 - m_W3;
	m_H3 = m_H1 - m_H2;
	m_sita_1 = 2 * asin(0.5*m_L3 / m_r2);

	double angle = atan(0.5*m_L2 / m_H3);
	double lth = sqrt(pow(0.5*m_L2, 2) + pow(m_H3, 2));
	angle = angle + acos(m_r2 / lth);

	m_sita_3 = 2 * PI - 2 * angle;
	m_sita_2 = angle - 0.5*m_sita_1;
}

void Bearing_Chock::Map2()
{
	Vec3 p(0, 0, 0);

	l_p1 = p;
	ms.MovePoint(l_p1, 0.5*m_W2, 2);

	l_p2 = p;
	ms.MovePoint(l_p2, m_W - 0.5*m_W3, 2);

	l_p3 = l_p2;
	ms.MovePoint(l_p3, -(0.5*m_L2 - 0.5*0.5*(m_L2 - m_L3)), 1);

	l_p4 = p;
	ms.MovePoint(l_p4, -m_L4, 1);
	ms.MovePoint(l_p4, m_W1, 2);

	l_p5 = p;
	ms.MovePoint(l_p5, m_W, 2);
	ms.MovePoint(l_p5, m_H3, 3);

}

varray<SplineVolume> Bearing_Chock::GetBearing()
{
	varray<SplineVolume> tmpvols;
	for (const auto& fv : m_fvs) {
		assert(fv.Iscomplete());
		varray<SplineVolume> vs = fv.GetVOls();
		for (const auto& v : vs) {
			tmpvols.push_back(v);
		}
	}
	return tmpvols;
}
