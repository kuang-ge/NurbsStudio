#pragma once

#include "stdafx.h"
#include "quadPart.h"
#include <math.h>


double sharpAngle = cos(10 * PI / 180); //������ֵ(��ʱȡ20��)
double Sharp = 10 * PI / 180;



QuadPart::QuadPart(varray<Spline> outLines, varray<Spline> inLines, varray<Spline> addLines, varray<int> outQuality, const varray<varray<Vec3>> removePol)
{
	this->removepoints_Pol = removePol;
	this->CreateAllPointAndLine(outLines, inLines, addLines, outQuality);
}

QuadPart::QuadPart(varray<Spline> outLines, varray<Spline> inLines, varray<Spline> addLines, varray<int> outQuality)
{
	this->CreateAllPointAndLine(outLines, inLines, addLines, outQuality);
}

void QuadPart::Init()
{
	this->in_Points.clear();
	this->in_NurbsLines.clear();
	this->out_Points.clear();
	this->out_NurbsLines.clear();
	this->all_Points.clear();
	this->all_NurbsLines.clear();
	this->quad_Polygon.clear();
	this->large_Polygon.clear();
	this->g_Adjlist.adjlist.clear();
	this->g_Adjlist.numEdges = 0;
	this->g_Adjlist.numNodes = 0;
}

//�ʷֲ���
void QuadPart::ExecuteQuadPart()
{
	this->OrderOutPointsAntioclock();
	this->CreateALGraph();//��������ͼ���ڽӱ�ṹ
	this->CreatePolygon();//���ɶ����
	this->OrderPolygens();//���������
	this->ResetPolygens();//��������
	this->RemovePolygon();//������Ƴ�

	varray<MyPolygon> endPol;
	varray<PartNurbsLine> endL;


	bool res;
	//ȥ�����
	if (sharpMode == 1) {
		res = this->QuadOperation(quad_Polygon, large_Polygon, all_NurbsLines, endPol, endL, sharpMode);
		if (!res)
		{
			endPol.clear();
			endL.clear();
			res = this->QuadOperation(quad_Polygon, large_Polygon, all_NurbsLines, endPol, endL, 0);
		}
	}
	else {
		res = this->QuadOperation(quad_Polygon, large_Polygon, all_NurbsLines, endPol, endL, sharpMode);
	}


	assert(res);
	//�����趨�ı��μ���������
	this->all_NurbsLines = endL;
	this->quad_Polygon = endPol;
	//��������
	this->ResetData();
}

void QuadPart::RemovePolygon()
{
	if (quad_Polygon.size() > 0) {

		for (int i = 0; i < removepoints_Pol.size();i++) {
			//��ȥ���ıߵ�
			varray<MyPolygon>::iterator it;
			for (it = quad_Polygon.begin(); it != quad_Polygon.end();) {
				if (removepoints_Pol[i].size() != (*it).p_Points.size()) {
					//����ζ����������ͬ
					it++;
					continue;
				}
				bool flag = false;
				for (int j = 0; j < removepoints_Pol[i].size(); j++) {
					flag = false;
					for (int k = 0; k < (*it).p_Points.size(); k++) {
						if (JudgeTwoPointsCoincide(removepoints_Pol[i][j], (*it).p_Points[k].vetx)) {
							flag = true;
							break;
						}
					}
					if (!flag) break;
				}
				if (flag) {
					remove_Polygon.push_back(*it);
					if (remove_genus[i]) {
						//�˶�����ǿ����
						int pos = 0;
						for (auto num : (*it).objectInAllLines) {
							//�趨��Ӧ�����߷ָ���
							//all_NurbsLines[num].ifSegLine = 1;
							all_NurbsLines[num].ifSegLine = remove_seg[i][pos];
							++pos;	
						}
					}
					quad_Polygon.erase(it);
					break;
				}
				it++;
			}
		}


		//for (auto p : removepoints_Pol) {
		//	//��ȥ���ıߵ�
		//	varray<MyPolygon>::iterator it;
		//	for (it = quad_Polygon.begin(); it != quad_Polygon.end();) {
		//		if (p.size() != (*it).p_Points.size()) {
		//			//����ζ����������ͬ
		//			it++;
		//			continue;
		//		}
		//		bool flag = false;
		//		for (int j = 0; j < p.size(); j++) {
		//			flag = false;
		//			for (int k = 0; k < (*it).p_Points.size(); k++) {
		//				if (JudgeTwoPointsCoincide(p[j], (*it).p_Points[k].vetx)) {
		//					flag = true;
		//					break;
		//				}
		//			}
		//			if (!flag) break;
		//		}
		//		if (flag) {
		//			remove_Polygon.push_back(*it);
		//			quad_Polygon.erase(it);
		//			break;
		//		}
		//		it++;
		//	}
		//}
	}

	if (large_Polygon.size() > 0) {

		for (int i = 0;i< removepoints_Pol.size(); i++) {
			//ȥ����ߵ�
			varray<MyPolygon>::iterator it;
			for (it = large_Polygon.begin(); it != large_Polygon.end();) {
				if (removepoints_Pol[i].size() != (*it).p_Points.size()) {
					//����ζ����������ͬ
					it++;
					continue;
				}
				bool flag = false;
				for (int j = 0; j < removepoints_Pol[i].size(); j++) {
					flag = false;
					for (int k = 0; k < (*it).p_Points.size(); k++) {
						if (JudgeTwoPointsCoincide(removepoints_Pol[i][j], (*it).p_Points[k].vetx)) {
							flag = true;
							break;
						}
					}
					if (!flag) break;
				}
				if (flag) {
					remove_Polygon.push_back(*it);
					if (remove_genus[i]) {
						//�˶�����ǿ����
						int pos = 0;
						for (auto num : (*it).objectInAllLines) {
							//�趨��Ӧ�����߿ɷָ���
							//all_NurbsLines[num].ifSegLine = 1;
							all_NurbsLines[num].ifSegLine = remove_seg[i][pos];
							++pos;
						}
					}
					large_Polygon.erase(it);
					break;
				}
				it++;
			}
		}



		//for (auto p : removepoints_Pol) {
		//	//ȥ����ߵ�
		//	varray<MyPolygon>::iterator it;
		//	for (it = large_Polygon.begin(); it != large_Polygon.end();) {
		//		if (p.size() != (*it).p_Points.size()) {
		//			//����ζ����������ͬ
		//			it++;
		//			continue;
		//		}
		//		bool flag = false;
		//		for (int j = 0; j < p.size(); j++) {
		//			flag = false;
		//			for (int k = 0; k < (*it).p_Points.size(); k++) {
		//				if (JudgeTwoPointsCoincide(p[j], (*it).p_Points[k].vetx)) {
		//					flag = true;
		//					break;
		//				}
		//			}
		//			if (!flag) break;
		//		}
		//		if (flag) {
		//			remove_Polygon.push_back(*it);
		//			large_Polygon.erase(it);
		//			break;
		//		}
		//		it++;
		//	}
		//}
	}

	assert(remove_Polygon.size() == removepoints_Pol.size());
}

void QuadPart::ResetData()
{
	map<int, int> mp;
	varray<PartNurbsLine> tmpAllLines;
	for (const auto &i : quad_Polygon) {
		for (auto j : i.objectInAllLines) {
			if (mp.find(j) == mp.end()) {
				//û�ҵ�key
				mp[j] = 1;
			}
		}		
	}

	//for (int i = 0; i < remove_Polygon.size();i++) {
	//	//if (remove_genus[i]) {
	//	//	//���������
	//	//	continue;
	//	//}
	//	for (auto j : remove_Polygon[i].objectInAllLines) {
	//		if (mp.find(j) == mp.end()) {
	//			mp[j] = 1;
	//		}
	//	}
	//}

	for (const auto &i : remove_Polygon) {
		for (auto j : i.objectInAllLines) {
			if (mp.find(j) == mp.end()) {
				mp[j] = 1;
			}
		}		
	}

	int length = tmpAllLines.size();
	map<int, int>::iterator it = mp.begin();
	for (; it != mp.end(); it++) {
		tmpAllLines.push_back(all_NurbsLines[it->first]);
		it->second = length++;
	}
	//�����趨���ж�������������
	for (auto &i : quad_Polygon) {
		for (int j = 0; j < i.objectInAllLines.size(); j++) {
			i.objectInAllLines[j] = mp[i.objectInAllLines[j]];
		}
	}
	for (auto &i : remove_Polygon) {
		for (int j = 0; j < i.objectInAllLines.size(); j++) {
			i.objectInAllLines[j] = mp[i.objectInAllLines[j]];
		}
	}
	mp.clear();
	this->all_NurbsLines = tmpAllLines;

	//�����趨����������
	varray<SubPoint> tmpOutP;
	varray<int> chooseLines;
	tmpOutP.push_back(out_Points[0]);//��һ����yֵ�϶���С

	//��ȡ��һ������
	for (int i = 0; i < tmpAllLines.size(); i++) {
		if (JudgeTwoPointsCoincide(*(tmpAllLines[i].m_CtrlPts.begin()), tmpOutP[0].vetx) ||
			JudgeTwoPointsCoincide(*(tmpAllLines[i].m_CtrlPts.end() - 1), tmpOutP[0].vetx))
			chooseLines.push_back(i);
	}
	
	if (chooseLines.size() == 0) assert(0);
	int tmpNum = 0;
	for (int i = 1; i < chooseLines.size(); i++) {
		Spline tmpLine = tmpAllLines[chooseLines[tmpNum]];
		Vec3 p = tmpOutP[0].vetx;
		Vec3 p1 = JudgeTwoPointsCoincide(p, tmpLine.m_CtrlPts[0]) ? *(tmpLine.m_CtrlPts.end() - 1) : *(tmpLine.m_CtrlPts.begin());
		Vec3 p2= JudgeTwoPointsCoincide(p, tmpAllLines[chooseLines[i]].m_CtrlPts[0]) ? *(tmpAllLines[chooseLines[i]].m_CtrlPts.end() - 1) : *(tmpAllLines[chooseLines[i]].m_CtrlPts.begin());
		Vec3 vec1 = p1 - p;
		Vec3 vec2 = p2 - p;
		double angle = GetAngle(vec1, vec2);
		if (angle > PI) {
			tmpNum = i;
		}
	}
	out_Number.push_back(chooseLines[tmpNum]);
	Vec3 p = JudgeTwoPointsCoincide((*(tmpOutP.end() - 1)).vetx, *(tmpAllLines[*(out_Number.end() - 1)].m_CtrlPts.begin())) ?
		*(tmpAllLines[*(out_Number.end() - 1)].m_CtrlPts.end() - 1) : *(tmpAllLines[*(out_Number.end() - 1)].m_CtrlPts.begin());
	tmpOutP.push_back(SubPoint(1, p));
	//ѭ����ȡ�µ�����������
	while (true)
	{
		chooseLines.clear();
		p = (tmpOutP.end() - 1)->vetx;
		Vec3 vec = (p - (tmpOutP.end() - 2)->vetx).Normalize();//��ǰ��͸�ǰһ�������
		//��ȡ�õ��������
		for (int i = 0; i < tmpAllLines.size(); i++) {
			if (JudgeTwoPointsCoincide(*(tmpAllLines[i].m_CtrlPts.begin()), p) ||
				JudgeTwoPointsCoincide(*(tmpAllLines[i].m_CtrlPts.end() - 1), p))
			{
				if (*(out_Number.end() - 1) == i)
					continue;
				chooseLines.push_back(i);
			}		
		}
		if (chooseLines.size() == 0) assert(0);
		tmpNum = 0;
		for (int i = 1; i < chooseLines.size(); i++) {
			Spline tmpLine = tmpAllLines[chooseLines[tmpNum]];
			Vec3 p1 = JudgeTwoPointsCoincide(p, tmpLine.m_CtrlPts[0]) ? *(tmpLine.m_CtrlPts.end() - 1) : *(tmpLine.m_CtrlPts.begin());
			Vec3 p2 = JudgeTwoPointsCoincide(p, tmpAllLines[chooseLines[i]].m_CtrlPts[0]) ? *(tmpAllLines[chooseLines[i]].m_CtrlPts.end() - 1) : *(tmpAllLines[chooseLines[i]].m_CtrlPts.begin());
			Vec3 vec1 = p1 - p;
			Vec3 vec2 = p2 - p;

			//2021.2.5�޸�
			double angle1 = GetAngle(vec, vec1);
			double angle2 = GetAngle(vec, vec2);
			if (angle1 <= PI &&angle2 <= PI) {
				//������С��180��
				if (angle2 < angle1)
					tmpNum = i;
			}
			else if (angle1 >= PI &&angle2 >= PI) {
				//����������180��
				if (angle2 < angle1)
					tmpNum = i;
			}
			else {
				if (angle2 > PI)
					tmpNum = i;
			}
			//double angle = GetAngle(vec1, vec2);
			//if (angle > PI) {
			//	tmpNum = i;
			//}
		}
		out_Number.push_back(chooseLines[tmpNum]);
		p = JudgeTwoPointsCoincide((*(tmpOutP.end() - 1)).vetx, *(tmpAllLines[*(out_Number.end() - 1)].m_CtrlPts.begin())) ?
			*(tmpAllLines[*(out_Number.end() - 1)].m_CtrlPts.end() - 1) : *(tmpAllLines[*(out_Number.end() - 1)].m_CtrlPts.begin());
		if (JudgeTwoPointsCoincide(p, tmpOutP[0].vetx)) break;
		tmpOutP.push_back(SubPoint(1, p));
	}
	//�滻��Ӧ����
	out_Points = tmpOutP;
}

void QuadPart::CreateAllPointAndLine(varray<Spline> outLines, varray<Spline> inLines, varray<Spline> addLines, varray<int> &outQuality)
{
	this->Init();
	//in_Points.clear();
	//out_Points.clear();
	//all_Points.clear();
	//in_NurbsLines.clear();
	//out_NurbsLines.clear();
	//all_NurbsLines.clear();
	//�ȶ�������������߷ָȥ��
	OperateOutandInLines(inLines, outLines, outQuality);

	SubPoint tempPoint;
	PartNurbsLine tempLine;
	int headFlag, tailFlag;

	//��Ϊ�����
	if (outQuality.size() == 0)
	{
		//������������Ϊ���������ߴ�С������ʼ��Ϊ1�������ɷָ�
		outQuality.resize(outLines.size(), 1);
	}

	//�����������㼯�ͱ߼�
	for (int i = 0; i < outLines.size(); i++)
	{
		headFlag = 1;
		tailFlag = 1;
		for (int j = 0; j < out_Points.size(); j++)
		{
			int ctpLength = outLines[i].m_CtrlPts.size();
			if (headFlag != 0)
			{

				if (JudgeTwoPointsCoincide(outLines[i].m_CtrlPts[0], out_Points[j].vetx))
				{
					headFlag = 0;
				}
			}
			
			if (tailFlag != 0)
			{
				if (JudgeTwoPointsCoincide(outLines[i].m_CtrlPts[ctpLength - 1], out_Points[j].vetx))
				{
					tailFlag = 0;
				}
			}
			
			if (headFlag == 0 && tailFlag == 0)
			{
				break;
			}

		}

		if (headFlag == 1)
		{
			tempPoint.ifOutPoint = 1;
			tempPoint.vetx = Vec3(outLines[i].m_CtrlPts[0].x, outLines[i].m_CtrlPts[0].y, outLines[i].m_CtrlPts[0].z);
			out_Points.push_back(tempPoint);
		}
		if (tailFlag == 1)
		{
			tempPoint.ifOutPoint = 1;
			int length = outLines[i].m_CtrlPts.size();
			tempPoint.vetx = Vec3(outLines[i].m_CtrlPts[length-1].x, outLines[i].m_CtrlPts[length - 1].y, outLines[i].m_CtrlPts[length - 1].z);
			out_Points.push_back(tempPoint);
		}

		tempLine.CreatPartNurbsLine(outLines[i], outQuality[i]); //�������ı߼��Ѿ�������Ϊ�ɷָ���(����outQuality������)
		out_NurbsLines.push_back(tempLine);
	}

	all_Points = out_Points;
	all_NurbsLines = out_NurbsLines;

	//�����������㼯�ͱ߼�
	for (int i = 0; i < inLines.size(); i++) 
	{
		//����0��ʾ�غ�
		headFlag = 1;
		tailFlag = 1;
		bool lineCon = false;
		Vec3 pHead = inLines[i].m_CtrlPts[0];
		Vec3 pTail = *(inLines[i].m_CtrlPts.end() - 1);

		//��ѯ�Ƿ���������������غ�
		if (in_Points.size() > 0)
		{
			for (int j = 0; j < in_Points.size(); j++)
			{
				int ctpLength = inLines[i].m_CtrlPts.size();
				if (headFlag != 0)
				{
					if (JudgeTwoPointsCoincide(inLines[i].m_CtrlPts[0], in_Points[j].vetx))
					{
						headFlag = 0;
					}
				}
				if (tailFlag != 0)
				{
					if (JudgeTwoPointsCoincide(inLines[i].m_CtrlPts[ctpLength - 1], in_Points[j].vetx))
					{
						tailFlag = 0;
					}
				}

				if (headFlag == 0 && tailFlag == 0)
				{
					break;
				}
			}
		}

		//��ѯ�Ƿ���������������غ�
		if (headFlag != 0) {
			for (auto p : out_Points) {
				if (JudgeTwoPointsCoincide(pHead, p.vetx)) {
					headFlag = 0;
					break;
				}
			}
		}

		if (tailFlag != 0) {
			for (auto p : out_Points) {
				if (JudgeTwoPointsCoincide(pTail, p.vetx)) {
					tailFlag = 0;
					break;
				}
			}
		}

		
		//��ѯ�Ƿ���������������غ�
		for (int m = 0; m < out_NurbsLines.size(); m++)
		{
			Vec3 p1, p2, p3, p4;
			p1 = out_NurbsLines[m].m_CtrlPts[0];
			p2 = *(out_NurbsLines[m].m_CtrlPts.end() - 1);
			p3 = inLines[i].m_CtrlPts[0];
			p4 = *(inLines[i].m_CtrlPts.end() - 1);
			bool compare1 = JudgeTwoPointsCoincide(p1, p3);
			bool compare2 = JudgeTwoPointsCoincide(p1, p4);
			bool compare3 = JudgeTwoPointsCoincide(p2, p3);
			bool compare4 = JudgeTwoPointsCoincide(p2, p4);
			if ((compare1 && compare4) || (compare2 && compare3)) //��ֱ���غ�
			{
				out_NurbsLines[m].ifSegLine = 0;  //���غϵ����������߶���Ϊ���ɷָ���
				lineCon = true;
				break;
			}
		}

		//��ѯ�Ƿ����������������غ�
		if ((!lineCon) && (in_NurbsLines.size()!=0))
		{
			for (int m = 0; m < in_NurbsLines.size(); m++)
			{
				Vec3 p1, p2, p3, p4;
				p1 = in_NurbsLines[m].m_CtrlPts[0];
				p2 = *(in_NurbsLines[m].m_CtrlPts.end() - 1);
				p3 = inLines[i].m_CtrlPts[0];
				p4 = *(inLines[i].m_CtrlPts.end() - 1);
				bool compare1 = JudgeTwoPointsCoincide(p1, p3);
				bool compare2 = JudgeTwoPointsCoincide(p1, p4);
				bool compare3 = JudgeTwoPointsCoincide(p2, p3);
				bool compare4 = JudgeTwoPointsCoincide(p2, p4);
				if ((compare1 && compare4) || (compare2 && compare3)) //��ֱ���غ�
				{
					lineCon = true;
					break;
				}
			}
		}
		



		//if (!lineCon)
		//{
		//	//��ѯ�Ƿ���������������غ�(��Ҫ�޸�)
		//	if (!(headFlag == 0 && tailFlag == 0))
		//	{
		//		for (int k = 0; k < out_Points.size(); k++)
		//		{
		//			int ctpLength = inLines[i].m_CtrlPts.size();
		//			if (headFlag != 0)
		//			{
		//				if (JudgeTwoPointsCoincide(inLines[i].m_CtrlPts[0], out_Points[k].vetx)
		//					/*inLines[i].m_CtrlPts[0].x == out_Points[k].vetx.x && inLines[i].m_CtrlPts[0].y == out_Points[k].vetx.y
		//					&& inLines[i].m_CtrlPts[0].z == out_Points[k].vetx.z*/)
		//				{
		//					headFlag = 0;
		//				}
		//			}
		//			if (tailFlag != 0)
		//			{
		//				if (JudgeTwoPointsCoincide(inLines[i].m_CtrlPts[ctpLength - 1], out_Points[k].vetx)
		//					/*inLines[i].m_CtrlPts[ctpLength - 1].x == out_Points[k].vetx.x && inLines[i].m_CtrlPts[ctpLength - 1].y == out_Points[k].vetx.y
		//					&& inLines[i].m_CtrlPts[ctpLength - 1].z == out_Points[k].vetx.z*/)
		//				{
		//					tailFlag = 0;
		//				}
		//			}

		//			if (headFlag == 0 && tailFlag == 0)
		//			{
		//				break;
		//			}
		//		}
		//	}

		//	if (headFlag == 1)
		//	{
		//		tempPoint.ifOutPoint = 0;
		//		tempPoint.vetx = Vec3(inLines[i].m_CtrlPts[0].x, inLines[i].m_CtrlPts[0].y, inLines[i].m_CtrlPts[0].z);
		//		in_Points.push_back(tempPoint);
		//	}
		//	if (tailFlag == 1)
		//	{
		//		tempPoint.ifOutPoint = 0;
		//		int length = inLines[i].m_CtrlPts.size();
		//		tempPoint.vetx = Vec3(inLines[i].m_CtrlPts[length - 1].x, inLines[i].m_CtrlPts[length - 1].y, inLines[i].m_CtrlPts[length - 1].z);
		//		in_Points.push_back(tempPoint);
		//	}

		//	tempLine.CreatPartNurbsLine(inLines[i], 0);   //������������Ϊ���ɷָ���
		//	in_NurbsLines.push_back(tempLine);
		//}

		if (headFlag == 1)
		{
			tempPoint.ifOutPoint = 0;
			tempPoint.vetx = Vec3(inLines[i].m_CtrlPts[0].x, inLines[i].m_CtrlPts[0].y, inLines[i].m_CtrlPts[0].z);
			in_Points.push_back(tempPoint);
		}
		if (tailFlag == 1)
		{
			tempPoint.ifOutPoint = 0;
			int length = inLines[i].m_CtrlPts.size();
			tempPoint.vetx = Vec3(inLines[i].m_CtrlPts[length - 1].x, inLines[i].m_CtrlPts[length - 1].y, inLines[i].m_CtrlPts[length - 1].z);
			in_Points.push_back(tempPoint);
		}

		if (!lineCon) {
			tempLine.CreatPartNurbsLine(inLines[i], 0);   //������������Ϊ���ɷָ���
			in_NurbsLines.push_back(tempLine);
		}
	}

	//�����Լ���ӵ����ߵ����������߼�����
	if (addLines.size() != 0)
	{
		for (int i = 0; i < addLines.size(); i++)
		{
			tempLine.CreatPartNurbsLine(addLines[i], 1);
			in_NurbsLines.push_back(tempLine);
		}
		
	}

	for (int i = 0; i < in_Points.size(); i++)
	{
		all_Points.push_back(in_Points[i]);
	}
	for (int i = 0; i < in_NurbsLines.size(); i++)
	{
		all_NurbsLines.push_back(in_NurbsLines[i]);
	}
}


void QuadPart::OperateOutandInLines(varray<Spline>& IL, varray<Spline>& OL, varray<int>& outQuality)
{
	SISLCurve *Icurve = nullptr;
	SISLCurve *Ocurve = nullptr;
	double *intpar1 = nullptr;
	double *intpar2 = nullptr;

	//�ȶ��ڲ�����ȥ��
	varray<Spline>::iterator it1, it2, tmpit;
	for (it1 = IL.begin(); it1 != IL.end();) {
		it2 = it1 + 1;
		for (; it2 != IL.end();) {
			if (JudgeTwoLinesCoincide(*it1, *it2)) {
				//�������غ�
				//tmpit = it2;
				//tmpit++;
				IL.erase(it2);
				//it2 = tmpit;
			}
			else
			{
				it2++;
			}
		}
		it1++;
	}


	//for (it1 = IL.begin(); it1 != IL.end(); it1++) {
	//	if (Icurve) freeCurve(Icurve);
	//	Icurve = NurbsLineToSislLine(*it1,2); //ת��������ΪSISL��ʽ
	//	for (it2 = OL.begin(); it2 != OL.end();) {
	//		if (Ocurve) freeCurve(Ocurve);
	//		Ocurve = NurbsLineToSislLine(*it2,2); //ת��������ΪSISL��ʽ
	//		if (intpar1)
	//		{
	//			delete[] intpar1;
	//			intpar1 = nullptr;
	//		}
	//		if (intpar2)
	//		{
	//			delete[] intpar2;
	//			intpar2 = nullptr;
	//		}
	//		int flag = TwoNurbsLineIntersect(Icurve, Ocurve, intpar1, intpar2); //����������,������״ֵ̬flag
	//		if (flag == 0) {
	//			//�޽��㣬����ѭ��
	//			it2++;
	//			continue;
	//		}
	//		else if (flag == 1)//���ڽ���
	//		{
	//			int Dvalue1 = CalU(*intpar1);
	//			int Dvalue2 = CalU(*intpar2);
	//			if ((Dvalue1 == 1 || Dvalue1 == 0) && (Dvalue2 == 1 || Dvalue2 == 0))//ֻ��һ���˵�Ϊ�ཻ�㣬�������κβ���
	//			{
	//				it2++;
	//				continue;
	//			}
	//			else //�������м䱻�ض�
	//			{

	//				varray<Spline> tmpLines;
	//				Vec3 p1 = (Dvalue1 == 0) ? (it1->m_CtrlPts[0]) : *(it1->m_CtrlPts.end() - 1);
	//				Vec3 p2 = it2->GetLinePoint(*intpar2);
	//				assert(JudgeTwoPointsCoincide(p1, p2));

	//				it2->Segmentation(*intpar2, tmpLines);
	//				if (tmpLines.size() == 0) assert(0);
	//				//����ԭ����ɾ���������߲������
	//				for (int i = 0; i < tmpLines.size(); i++) {
	//					OL.insert(it2, tmpLines[0]);
	//					it2++;
	//					if (it2 == NULL) assert(0);
	//				}
	//				//����it2 eraseӦ��ָ��Ļ�����һ���Ĳ���
	//				OL.erase(it2);		
	//				if (it2 == NULL) assert(0);
	//			}
	//		}
	//		else if (flag == 2) {
	//			//�����غ��߶�
	//			if (JudgeTwoLinesCoincide(*it1, *it2)) //�������غϣ�����Ҫ�����޸Ĳ���
	//			{
	//				it2++;
	//				continue;
	//			}
	//			//ȡ�㣬���Ե��Ƿ�����������
	//			double* p1 = Point3dToSislPoint(it1->m_CtrlPts[0], 2);
	//			double* p2 = Point3dToSislPoint(*(it1->m_CtrlPts.end() - 1), 2);
	//			double u1 = PointIntersectNurbsLine(p1, Ocurve, 2);
	//			double u2 = PointIntersectNurbsLine(p2, Ocurve, 2);
	//			int flag2 = CalU(u1);
	//			if (flag2 == 0 || flag2 == 1) {
	//				//��һ�����Ƕ˵�
	//				varray<Spline> tmpLines;
	//				it2->Segmentation(u2, tmpLines);
	//				assert(tmpLines.size());
	//				//����ԭ����ɾ���������߲������
	//				for (int i = 0; i < tmpLines.size(); i++) {
	//					OL.insert(it2, tmpLines[0]);
	//					it2++;
	//					if (it2 == NULL) assert(0);
	//				}
	//				OL.erase(it2);
	//				if (it2 == NULL) assert(0);

	//				delete[] p1;
	//				delete[] p2;
	//				p1 = NULL;
	//				p2 = NULL;
	//				continue;
	//			}
	//			flag2 = CalU(u2);
	//			if (flag2 == 0 || flag2 == 1) {
	//				//��һ�����Ƕ˵�
	//				varray<Spline> tmpLines;
	//				it2->Segmentation(u1, tmpLines);
	//				//����ԭ����ɾ���������߲������
	//				for (int i = 0; i < tmpLines.size(); i++) {
	//					OL.insert(it2, tmpLines[0]);
	//					it2++;
	//					if (it2 == NULL) assert(0);
	//				}
	//				OL.erase(it2);
	//				if (it2 == NULL) assert(0);

	//				delete[] p1;
	//				delete[] p2;
	//				p1 = NULL;
	//				p2 = NULL;
	//				continue;
	//			}
	//			//if()//�������߾�����һ�����غ�,��������Ȳ�����

	//			if (u1 > 0 && u1 < 1 && u2>0 && u2 < 1)//�����߰���������
	//			{
	//				varray<Spline> tmpLines;
	//				if (u1 < u2)
	//				{
	//					it2->Segmentation(u1, tmpLines);
	//					assert(tmpLines.size());
	//					Spline tmp = tmpLines[1];
	//					freeCurve(Ocurve);
	//					Ocurve = NurbsLineToSislLine(tmp,2);
	//					u2 = PointIntersectNurbsLine(p2, Ocurve, 2);
	//					varray<Spline> tmpLines2;
	//					tmp.Segmentation(u2, tmpLines2);
	//					assert(tmpLines2.size());
	//					tmpLines.pop_back();
	//					tmpLines.push_back(tmpLines2[0]);
	//					tmpLines.push_back(tmpLines2[1]);
	//					tmpLines2.clear();
	//				}
	//				else
	//				{
	//					it2->Segmentation(u2, tmpLines);
	//					assert(tmpLines.size());
	//					Spline tmp = tmpLines[1];
	//					freeCurve(Ocurve);
	//					Ocurve = NurbsLineToSislLine(tmp,2);
	//					u1 = PointIntersectNurbsLine(p1, Ocurve, 2);
	//					varray<Spline> tmpLines2;
	//					tmp.Segmentation(u1, tmpLines2);
	//					assert(tmpLines2.size());
	//					tmpLines.pop_back();
	//					tmpLines.push_back(tmpLines2[0]);
	//					tmpLines.push_back(tmpLines2[1]);
	//					tmpLines2.clear();
	//				}

	//				//����ԭ����ɾ���������߲������
	//				for (int i = 0; i < tmpLines.size(); i++) {
	//					OL.insert(it2, tmpLines[0]);
	//					it2++;
	//					if (it2 == NULL) assert(0);
	//				}
	//				OL.erase(it2);
	//				if (it2 == NULL) assert(0);

	//				delete[] p1;
	//				delete[] p2;
	//				p1 = NULL;
	//				p2 = NULL;
	//				continue;
	//			}
	//		}
	//	}
	//}
	

	for (int i = 0; i < IL.size(); i++)
	{
		if (Icurve)
		{
			freeCurve(Icurve);
		}
		Icurve = NurbsLineToSislLine(IL[i]); //ת��������ΪSISL��ʽ
		for (int j = 0; j < OL.size(); j++)
		{
			if (Ocurve)
			{
				freeCurve(Ocurve);
			}
			Ocurve = NurbsLineToSislLine(OL[j]); //ת��������ΪSISL��ʽ

			if(intpar1)
			{
				delete[] intpar1;
				intpar1 = nullptr;
			}
			if (intpar2)
			{
				delete[] intpar2;
				intpar2 = nullptr;
			}
			int flag = TwoNurbsLineIntersect(Icurve, Ocurve, intpar1, intpar2); //����������,������״ֵ̬flag

			if (flag == 0) //�޽��㣬����ѭ��
			{
				continue;
			}
			else if (flag == 1)//���ڽ���
			{
				int Dvalue1 = CalU(*intpar1);
				int Dvalue2 = CalU(*intpar2);
				if ((Dvalue1 == 1 || Dvalue1 == 0) && (Dvalue2 == 1 || Dvalue2 == 0))//ֻ��һ���˵�Ϊ�ཻ�㣬�������κβ���
				{
					continue;
				}
				else //�������м䱻�ض�
				{
					varray<Spline> tmpLines;
					Vec3 p1 = (Dvalue1 == 0) ? (IL[i].m_CtrlPts[0]) : *(IL[i].m_CtrlPts.end() - 1);
					Vec3 p2 = OL[j].GetLinePoint(*intpar2);
					assert(JudgeTwoPointsCoincide(p1, p2));

					OL[j].Segmentation(*intpar2, tmpLines);
					//����ԭ����ɾ���������߲������
					if (j == OL.size() - 1)
					{
						OL.pop_back();
						OL.push_back(tmpLines[0]);
						OL.push_back(tmpLines[1]);
						//�����ɷָ���
						outQuality.pop_back();
						outQuality.push_back(1);
						outQuality.push_back(1);
						j += 1;
					}
					else
					{
						OL.erase(&OL[j]);
						outQuality.erase(&outQuality[j]);
						for (int m = 0; m < tmpLines.size(); m++)
						{
							OL.insert(&OL[j], tmpLines[m]);
							//
							outQuality.insert(&outQuality[j], 1);
							j++;
						}
						j--;
					}
				}
			}
			else if (flag == 2)//�����غ��߶�
			{
				if (JudgeTwoLinesCoincide(IL[i], OL[j])) //�������غϣ�����Ҫ�����޸Ĳ���
				{
					outQuality[j] = 0;
					continue;
				}

				//ȡ�㣬���Ե��Ƿ�����������
				double* p1 = Point3dToSislPoint(IL[i].m_CtrlPts[0], 3);
				double* p2 = Point3dToSislPoint(*(IL[i].m_CtrlPts.end() - 1), 3);
				double u1 = PointIntersectNurbsLine(p1, Ocurve, 3);
				double u2 = PointIntersectNurbsLine(p2, Ocurve, 3);
				int flag2 = CalU(u1);
				if (flag2 == 0 || flag2 == 1)//��һ�����Ƕ˵�
				{
					varray<Spline> tmpLines;
					OL[j].Segmentation(u2, tmpLines);
					if (j == OL.size() - 1) //ԭ���߲����ɾ��
					{
						OL.pop_back();
						OL.push_back(tmpLines[0]);
						OL.push_back(tmpLines[1]);
						outQuality.pop_back();
						if (JudgeTwoLinesCoincide(IL[i], tmpLines[0])) {
							outQuality.push_back(0);
							outQuality.push_back(1);
						}
						else if (JudgeTwoLinesCoincide(IL[i], tmpLines[1])) {
							outQuality.push_back(1);
							outQuality.push_back(0);
						}
						else
						{
							assert(0);
						}
						j += 1;
					}
					else
					{
						OL.erase(&OL[j]);
						outQuality.erase(&outQuality[j]);
						bool flagU1 = false;
						for (int m = 0; m < tmpLines.size(); m++)
						{
							OL.insert(&OL[j], tmpLines[m]);
							if (JudgeTwoLinesCoincide(IL[i], tmpLines[m])) {
								flagU1 = true;
								outQuality.insert(&outQuality[j], 0);
							}
							else {
								outQuality.insert(&outQuality[j], 1);
							}
							j++;
						}
						assert(flagU1);
						j--;
					}
					delete[] p1;
					delete[] p2;
					p1 = NULL;
					p2 = NULL;
					continue;
				}
				flag2 = CalU(u2);
				if (flag2 == 0 || flag2 == 1)//��һ�����Ƕ˵�
				{
					varray<Spline> tmpLines;
					OL[j].Segmentation(u1, tmpLines);
					if (j == OL.size() - 1) //ԭ���߲����ɾ��
					{
						OL.pop_back();
						OL.push_back(tmpLines[0]);
						OL.push_back(tmpLines[1]);

						outQuality.pop_back();
						if (JudgeTwoLinesCoincide(IL[i], tmpLines[0])) {
							outQuality.push_back(0);
							outQuality.push_back(1);
						}
						else if (JudgeTwoLinesCoincide(IL[i], tmpLines[1])) {
							outQuality.push_back(1);
							outQuality.push_back(0);
						}
						else
						{
							assert(0);
						}

						j += 1;
					}
					else
					{
						OL.erase(&OL[j]);
						outQuality.erase(&outQuality[j]);
						bool flagU1 = false;
						for (int m = 0; m < tmpLines.size(); m++)
						{
							OL.insert(&OL[j], tmpLines[m]);
							if (JudgeTwoLinesCoincide(IL[i], tmpLines[m])) {
								flagU1 = true;
								outQuality.insert(&outQuality[j], 0);
							}
							else {
								outQuality.insert(&outQuality[j], 1);
							}
							j++;
						}
						assert(flagU1);
						j--;
					}
					delete[] p1;
					delete[] p2;
					p1 = NULL;
					p2 = NULL;
					continue;
				}

				//if()//�������߾�����һ�����غ�,��������Ȳ�����

				if (u1 > 0 && u1 < 1 && u2>0 && u2 < 1)//�����߰���������
				{
					varray<Spline> tmpLines;
					if (u1 < u2)
					{
						OL[j].Segmentation(u1, tmpLines);
						Spline tmp = tmpLines[1];
						freeCurve(Ocurve);
						Ocurve = NurbsLineToSislLine(tmp);
						u2 = PointIntersectNurbsLine(p2, Ocurve, 3);
						varray<Spline> tmpLines2;
						tmp.Segmentation(u2, tmpLines2);
						tmpLines.pop_back();
						tmpLines.push_back(tmpLines2[0]);
						tmpLines.push_back(tmpLines2[1]);
						tmpLines2.clear();
					}
					else
					{
						OL[j].Segmentation(u2, tmpLines);
						Spline tmp = tmpLines[1];
						freeCurve(Ocurve);
						Ocurve = NurbsLineToSislLine(tmp);
						u1 = PointIntersectNurbsLine(p1, Ocurve, 3);
						varray<Spline> tmpLines2;
						tmp.Segmentation(u1, tmpLines2);
						tmpLines.pop_back();
						tmpLines.push_back(tmpLines2[0]);
						tmpLines.push_back(tmpLines2[1]);
						tmpLines2.clear();
					}
					//���в���ɾ������
					if (j == OL.size() - 1) //ԭ���߲����ɾ��
					{
						OL.pop_back();
						outQuality.pop_back();
						for (int m = 0; m < tmpLines.size(); m++)
						{
							OL.push_back(tmpLines[m]);
							if (m == 1) {
								outQuality.push_back(0);
							}
							else {
								outQuality.push_back(1);
							}
							j++;
						}
						j--;
					}
					else
					{
						OL.erase(&OL[j]);
						outQuality.erase(&outQuality[j]);
						for (int m = 0; m < tmpLines.size(); m++)
						{
							OL.insert(&OL[j], tmpLines[m]);
							if (m == 1) {
								outQuality.insert(&outQuality[j], 0);
							}
							else {
								outQuality.insert(&outQuality[j], 1);
							}
							j++;
						}
						j--;
					}
					delete[] p1;
					delete[] p2;
					p1 = NULL;
					p2 = NULL;
					continue;
				}
			}
		}
	}

	if (Icurve)
	{
		freeCurve(Icurve);
		Icurve = nullptr;
	}
	if (Ocurve)
	{
		freeCurve(Ocurve);
		Ocurve = nullptr;
	}
	if (intpar1)
	{
		delete[] intpar1;
		intpar1 = nullptr;
	}
	if (intpar2)
	{
		delete[] intpar2;
		intpar2 = nullptr;
	}	
}

//�����������㰴��ʱ���������
void QuadPart::OrderOutPointsAntioclock()
{
	int firstPoint = 0;
	int line1 = -1;
	int line2 = -1;
	Vec3 secondPoint1, secondPoint2;
	varray<SubPoint> tempPoints;
	SubPoint tempPoint;
	//�ҵ�y����ֵ��С�ĵ�
	for (int i = 0; i < out_Points.size(); i++)
	{
		if (out_Points[i].vetx.y <= out_Points[firstPoint].vetx.y)
		{
			firstPoint = i;
		}
	}
	tempPoints.push_back(out_Points[firstPoint]);

	//�ҵ��õ�����ϵ��������
	for (int i = 0; i < out_NurbsLines.size(); i++)
	{
		int cptLength = out_NurbsLines[i].m_CtrlPts.size();
		if (JudgeTwoPointsCoincide(out_NurbsLines[i].m_CtrlPts[0], out_Points[firstPoint].vetx))
		{
			if (line1 == -1 && line2 == -1)
			{
				secondPoint1 = Vec3(out_NurbsLines[i].m_CtrlPts[cptLength - 1].x, out_NurbsLines[i].m_CtrlPts[cptLength - 1].y,
					out_NurbsLines[i].m_CtrlPts[cptLength - 1].z);
				line1 = i;
			}
			else
			{
				secondPoint2= Vec3(out_NurbsLines[i].m_CtrlPts[cptLength - 1].x, out_NurbsLines[i].m_CtrlPts[cptLength - 1].y,
					out_NurbsLines[i].m_CtrlPts[cptLength - 1].z);
				line2 = i;
			}
		}

		else if (JudgeTwoPointsCoincide(out_NurbsLines[i].m_CtrlPts[cptLength - 1], out_Points[firstPoint].vetx))
		{
			if (line1 == -1 && line2 == -1)
			{
				secondPoint1 = Vec3(out_NurbsLines[i].m_CtrlPts[0].x, out_NurbsLines[i].m_CtrlPts[0].y,
					out_NurbsLines[i].m_CtrlPts[0].z);
				line1 = i;
			}
			else
			{
				secondPoint2 = Vec3(out_NurbsLines[i].m_CtrlPts[0].x, out_NurbsLines[i].m_CtrlPts[0].y,
					out_NurbsLines[i].m_CtrlPts[0].z);
				line2 = i;
			}
		}
		
		if (line1 != -1 && line2 != -1)
		{
			break;
		}
	}

	Vec3 vec1 = secondPoint1 - out_Points[firstPoint].vetx;
	Vec3 vec2 = secondPoint2 - out_Points[firstPoint].vetx;
	vec1.Normalize(); //��λ��
	vec2.Normalize();
	Vec3 cross = vec1.Cross(vec2);//���������
	if (cross.z > 0)
	{
		tempPoint.ifOutPoint = 1;
		tempPoint.vetx = secondPoint1;
		tempPoints.push_back(tempPoint);
	}
	else if (cross.z < 0)
	{
		tempPoint.ifOutPoint = 1;
		tempPoint.vetx = secondPoint2;
		tempPoints.push_back(tempPoint);
	}
	else if (fabs(cross.z - 0) < 1e-5)
	{
		if (vec1.x > vec2.x)
		{
			tempPoint.ifOutPoint = 1;
			tempPoint.vetx = secondPoint1;
			tempPoints.push_back(tempPoint);
		}
		else
		{
			tempPoint.ifOutPoint = 1;
			tempPoint.vetx = secondPoint2;
			tempPoints.push_back(tempPoint);
		}
	}

	//ѭ�������ʱ��˳��ĵ�
	for (int i = 0; i < out_NurbsLines.size(); i++)
	{
		int tempPointsLength = tempPoints.size();
		int cptLength = out_NurbsLines[i].m_CtrlPts.size();
		Vec3 headPoint = out_NurbsLines[i].m_CtrlPts[0];
		Vec3 tailPoint = out_NurbsLines[i].m_CtrlPts[cptLength - 1];
		Vec3 beforeTempPoint = tempPoints[tempPointsLength - 2].vetx;
		tempPoint = tempPoints[tempPointsLength - 1];

		if (JudgeTwoPointsCoincide(headPoint, tempPoint.vetx))
		{
			if (JudgeTwoPointsCoincide(tailPoint, tempPoints[0].vetx))
			{
				if (tempPointsLength == out_Points.size())
				{
					break;
				}
				else
				{
					continue;
				}
				
			}

			if (!JudgeTwoPointsCoincide(tailPoint, beforeTempPoint))
			{
				tempPoint.ifOutPoint = 1;
				tempPoint.vetx = tailPoint;
				tempPoints.push_back(tempPoint);
				i = -1;
			}
		}

		else if (JudgeTwoPointsCoincide(tailPoint, tempPoint.vetx))
		{
			if (JudgeTwoPointsCoincide(headPoint, tempPoints[0].vetx))
			{
				if (tempPointsLength == out_Points.size())
				{
					break;
				}
				else
				{
					continue;
				}
			}

			if (!JudgeTwoPointsCoincide(headPoint, beforeTempPoint))
			{
				tempPoint.ifOutPoint = 1;
				tempPoint.vetx = headPoint;
				tempPoints.push_back(tempPoint);
				i = -1;
			}
		}
	}

	//�滻�������㼯��
	if (tempPoints.size() == out_Points.size())
	{
		out_Points.clear();
		for (int i = 0; i < tempPoints.size(); i++)
		{
			out_Points.push_back(tempPoints[i]);
		}
	}

	//�˴��ɼ���������ȵľ����������
}


//����ͼ���ڽӱ�ṹ
void QuadPart::CreateALGraph()
{
	int i, j, k, m, CptLength;
	//EdgeNode *e,*p;
	VertexNode tempNode;
	Vec3 headPoint, tailPoint;
	g_Adjlist.numNodes = out_Points.size();
	g_Adjlist.numEdges = out_NurbsLines.size();

	for (int i = 0; i < g_Adjlist.numNodes; i++)   //����������������Ϣ(��ʱ������õ�)�������������Ķ����ͱ߱�
	{
		if (i < g_Adjlist.numNodes - 1)
		{
			tempNode.vexdata = out_Points[i];
			EdgeNode* e = new EdgeNode();
			e->adjvex = i + 1;
			e->next = NULL;
			for (j = 0; j < g_Adjlist.numEdges; j++)
			{
				k = out_NurbsLines[j].m_CtrlPts.size();
				if (JudgeTwoPointsCoincide(out_NurbsLines[j].m_CtrlPts[0], out_Points[i].vetx)
					&& JudgeTwoPointsCoincide(out_NurbsLines[j].m_CtrlPts[k - 1], out_Points[i + 1].vetx))
				{
					e->EdgeLine = out_NurbsLines[j];
				}
				else if (JudgeTwoPointsCoincide(out_NurbsLines[j].m_CtrlPts[0], out_Points[i + 1].vetx)
					&& JudgeTwoPointsCoincide(out_NurbsLines[j].m_CtrlPts[k - 1], out_Points[i].vetx))
				{
					e->EdgeLine = out_NurbsLines[j];
				}
				//�ɼ������߶�����������ж�
			}
			tempNode.firstedge = e;
			e = NULL;
			delete e;
		}

		else
		{
			tempNode.vexdata = out_Points[i];
			EdgeNode* e = new EdgeNode();
			e->adjvex = 0;
			e->next = NULL;
			for (j = 0; j < g_Adjlist.numEdges; j++)
			{
				k = out_NurbsLines[j].m_CtrlPts.size();
				if (JudgeTwoPointsCoincide(out_NurbsLines[j].m_CtrlPts[0], out_Points[i].vetx)
					&& JudgeTwoPointsCoincide(out_NurbsLines[j].m_CtrlPts[k - 1], out_Points[0].vetx))
				{
					e->EdgeLine = out_NurbsLines[j];
				}
				else if (JudgeTwoPointsCoincide(out_NurbsLines[j].m_CtrlPts[0], out_Points[0].vetx)
					&& JudgeTwoPointsCoincide(out_NurbsLines[j].m_CtrlPts[k - 1], out_Points[i].vetx))
				{
					e->EdgeLine = out_NurbsLines[j];
				}
				//�ɼ������߶�����������ж�
			}
			tempNode.firstedge = e;
			e = NULL;
			delete e;
		}
		g_Adjlist.adjlist.push_back(tempNode);
		
	}

	g_Adjlist.numNodes = in_Points.size();
	for (int i = 0; i < g_Adjlist.numNodes; i++)  //����������������Ϣ�����������������
	{
		tempNode.vexdata = in_Points[i];  //���붥����Ϣ
		tempNode.firstedge = NULL;   //���߱���Ϊ�ձ�
		g_Adjlist.adjlist.push_back(tempNode);
	}

	g_Adjlist.numNodes = g_Adjlist.adjlist.size();
	g_Adjlist.numEdges = in_NurbsLines.size();

	for (k = 0; k < g_Adjlist.numEdges; k++)  //�����ڱ߱߱�
	{
		CptLength = in_NurbsLines[k].m_CtrlPts.size();
		headPoint = Vec3(in_NurbsLines[k].m_CtrlPts[0].x, in_NurbsLines[k].m_CtrlPts[0].y, in_NurbsLines[k].m_CtrlPts[0].z);
		tailPoint = Vec3(in_NurbsLines[k].m_CtrlPts[CptLength-1].x, in_NurbsLines[k].m_CtrlPts[CptLength-1].y, in_NurbsLines[k].m_CtrlPts[CptLength-1].z);
		
		i = -1;
		j = -1;
		for (m = 0; m < g_Adjlist.numNodes; m++)
		{
			if (JudgeTwoPointsCoincide(g_Adjlist.adjlist[m].vexdata.vetx, headPoint))
			{
				i = m;
			}
			else if(JudgeTwoPointsCoincide(g_Adjlist.adjlist[m].vexdata.vetx, tailPoint))
			{
				j = m;
			}
			if (i != -1 && j != -1)
			{
				break;
			}
		}
		assert(!(i == -1 || j == -1));

		//�趨��i�ı߱�ṹ
		EdgeNode* e = new EdgeNode();  
		e->adjvex = j;
		e->next = NULL;
		e->EdgeLine = in_NurbsLines[k];
		EdgeNode* p;
		p = g_Adjlist.adjlist[i].firstedge;
		if (p == NULL)
		{
			g_Adjlist.adjlist[i].firstedge = e;
			
		}
		else if (p != NULL) //���в������
		{
			e->next = p;
			g_Adjlist.adjlist[i].firstedge = e;
		}
		p = NULL;
		e = NULL;
		delete e;

		//�趨��j�ı߱�ṹ
		EdgeNode* r = new EdgeNode();
		r->adjvex = i;
		r->next = NULL;
		r->EdgeLine = in_NurbsLines[k];
		p = g_Adjlist.adjlist[j].firstedge;
		if (p == NULL)
		{
			g_Adjlist.adjlist[j].firstedge = r;
			
		}
		else if (p != NULL) //���в������
		{
			r->next = p;
			g_Adjlist.adjlist[j].firstedge = r;
		}
		p = NULL;
		r = NULL;
		delete r;
	}
	g_Adjlist.numEdges = all_NurbsLines.size();
}
//��������
void QuadPart::CreatePolygon()
{
	int numAreas = g_Adjlist.numEdges - g_Adjlist.numNodes + 1;
	while (large_Polygon.size() < numAreas)
	{
		for (int i = 0; i < g_Adjlist.numNodes; i++)
		{
			if (g_Adjlist.adjlist[i].firstedge != NULL)
			{
				MyPolygon tempPolygon;
				int head = i;
				Vec3 beginVetx, endVetx;
				VertexNode* p ;
				VertexNode* q ;
				EdgeNode* tempEdge ;
				EdgeNode* tempEdge2 ;
				p = &(g_Adjlist.adjlist[i]);
				tempPolygon.p_Points.push_back(p->vexdata);
				//�ڶ������ѡȡ
				q = p;
				p = &(g_Adjlist.adjlist[q->firstedge->adjvex]);
				tempPolygon.p_Points.push_back(p->vexdata);
				tempPolygon.p_NurbsLines.push_back(q->firstedge->EdgeLine);
				//ɾ����һ����ָ��ڶ�����ı߱�
				if (q->firstedge->next == NULL)
				{
					tempEdge = q->firstedge;
					q->firstedge = NULL;
					//free(tempEdge);
					delete tempEdge;   
					tempEdge = NULL;
				}
				else if (q->firstedge->next != NULL)
				{
					tempEdge = q->firstedge;
					q->firstedge = q->firstedge->next;
					//free(tempEdge);
					delete tempEdge;   
					tempEdge = NULL;
				}
				beginVetx = p->vexdata.vetx - q->vexdata.vetx;
				beginVetx=beginVetx.Normalize();

				//ѭ����ȡ����κ���ĵ㣬ֱ����һ��������ʼ��
				while (!(p->vexdata.vetx == g_Adjlist.adjlist[head].vexdata.vetx))  //��������е����⣬��Ҫ�޸���
				{
					int tag = 0;
					int caseChoice = -1;
					int nextAdj=0;
					//if (p->firstedge->adjvex != head)
					
						varray<Vec3> crossCollect;
						varray<double> dotCollect;
						VertexNode tempNode;
						tempEdge = p->firstedge;
						while (tempEdge)        //���������ڽӱߵĵ�˺Ͳ��
						{
							if (JudgeTwoPointsCoincide(g_Adjlist.adjlist[tempEdge->adjvex].vexdata.vetx, q->vexdata.vetx))
							{
								dotCollect.push_back(-1.0);
								crossCollect.push_back(Vec3(0, 0, -2.0));
								tempEdge = tempEdge->next;
								continue;
							}

							//������Ҫ��һ�����ߵ��ж�(������ͬ)

							else
							{
								tempNode = g_Adjlist.adjlist[tempEdge->adjvex];
								endVetx = tempNode.vexdata.vetx - p->vexdata.vetx;
								endVetx=endVetx.Normalize();
								double dot = beginVetx.Dot(endVetx);
								Vec3 cross = beginVetx.Cross(endVetx);
								dotCollect.push_back(dot);
								crossCollect.push_back(cross);
								
							}	
							tempEdge = tempEdge->next;
						}
						for (int i = 0; i < crossCollect.size(); i++)
						{
							if (crossCollect[i].z >= 0)
							{
								tag++;
							}
						}

						if (tag == 0)   
						{
							caseChoice = 0;
						}
						else if (tag==1)
						{
							caseChoice = 1;
						}
						else if (tag > 1)
						{
							caseChoice = 2;
						}

						nextAdj = 0;  //��ʼΪ0
						switch (caseChoice)
						{
						case 0:
							for (int i = 0; i < dotCollect.size(); i++)
							{
								if (dotCollect[i] >= dotCollect[nextAdj])
								{
									nextAdj = i;
								}
							}
							break;

						case 1:
							for (int i = 0; i < dotCollect.size(); i++)
							{
								if (crossCollect[i].z >= 0)
								{
									nextAdj = i;
								}
							}
							break;

						case 2:

							for (int i = 0; i < dotCollect.size(); i++)
							{
								if (crossCollect[i].z < 0)
								{
									dotCollect[i] = 2.0;
								}
								else if (crossCollect[i].z >= 0)
								{
									if (dotCollect[i] <= dotCollect[nextAdj] && dotCollect[i] != -1.0)
									{
										nextAdj = i;
									}
								}
							}
							break;
						}

						q = p;   
						tempEdge2 = p->firstedge;
						if (nextAdj > 0)
						{
							for (int i = 0; i < nextAdj; i++)
							{
								tempEdge2 = tempEdge2->next;
							}
						}
						p = &(g_Adjlist.adjlist[tempEdge2->adjvex]);
						if (!JudgeTwoPointsCoincide(p->vexdata.vetx, g_Adjlist.adjlist[head].vexdata.vetx))
						{
							tempPolygon.p_Points.push_back(p->vexdata);
						}
						tempPolygon.p_NurbsLines.push_back(tempEdge2->EdgeLine);

						//ɾ���õ�ָ����һ����ı߱�
						if (nextAdj == 0)
						{
							tempEdge2 = q->firstedge;
							if (tempEdge2->next == NULL)
							{
								q->firstedge = NULL;
								delete tempEdge2;   
								tempEdge2 = NULL;
							}
							else if (tempEdge2->next != NULL)
							{
								q->firstedge = tempEdge2->next;
								delete tempEdge2;   
								tempEdge2 = NULL;
							}
						}

						else if (nextAdj == 1)
						{
							tempEdge = q->firstedge;
							tempEdge2 = tempEdge->next;
							if (tempEdge2->next == NULL)
							{
								tempEdge->next = NULL;
								delete tempEdge2;   
								tempEdge2 = NULL;
								tempEdge = NULL;
							}
							else if (tempEdge2->next != NULL)
							{
								tempEdge->next = tempEdge2->next;
								delete tempEdge2;  
								tempEdge2 = NULL;
								tempEdge = NULL;
							}
						}

						else if (nextAdj > 1)
						{
							tempEdge = q->firstedge;
							tempEdge2 = tempEdge->next;
							for (int i = 1; i < nextAdj; i++)
							{
								tempEdge = tempEdge->next;
								tempEdge2 = tempEdge2->next;
							}
							if (tempEdge2->next == NULL)
							{
								tempEdge->next = NULL;
								delete tempEdge2;   
								tempEdge2 = NULL;
								tempEdge = NULL;
							}
							else if (tempEdge2->next != NULL)
							{
								tempEdge->next = tempEdge2->next;
								delete tempEdge2;   
								tempEdge2 = NULL;
								tempEdge = NULL;
							}
						}
						beginVetx = p->vexdata.vetx - q->vexdata.vetx;
						beginVetx=beginVetx.Normalize();
				}

				large_Polygon.push_back(tempPolygon);
				
				i = -1;  //���´ӵ�һ���㿪ʼ�����Ƿ������һ�����ӵ�
			}
		}
	}
}
//��������
Vec3 QuadPart::CalGravity(varray<SubPoint> polygenPoints)
{
	double martix, xaddx, yaddy, gx, gy, gz, denominator, x_numerator, y_numerator;
	varray<double> martixCollect, xCollect, yCollect;
	Vec3 gcoordinate;
	int pointsNumber = polygenPoints.size();
	for (int i = 0; i < pointsNumber; i++)
	{
		if (i == (pointsNumber-1))
		{
			martix = (polygenPoints[i].vetx.x*polygenPoints[0].vetx.y) - (polygenPoints[i].vetx.y*polygenPoints[0].vetx.x);
			xaddx = polygenPoints[i].vetx.x + polygenPoints[0].vetx.x;
			yaddy= polygenPoints[i].vetx.y + polygenPoints[0].vetx.y;
			martixCollect.push_back(martix);
			xCollect.push_back(xaddx);
			yCollect.push_back(yaddy);
		}

		else
		{
			martix = (polygenPoints[i].vetx.x*polygenPoints[i+1].vetx.y) - (polygenPoints[i].vetx.y*polygenPoints[i+1].vetx.x);
			xaddx = polygenPoints[i].vetx.x + polygenPoints[i+1].vetx.x;
			yaddy = polygenPoints[i].vetx.y + polygenPoints[i+1].vetx.y;
			martixCollect.push_back(martix);
			xCollect.push_back(xaddx);
			yCollect.push_back(yaddy);
		}
	}

	x_numerator = 0;
	y_numerator = 0;
	denominator = 0;

	for (int i = 0; i < martixCollect.size(); i++)
	{
		x_numerator += (xCollect[i] * martixCollect[i]);
		y_numerator += (yCollect[i] * martixCollect[i]);
		denominator += martixCollect[i];
	}
	denominator = denominator * 3;
	gx = x_numerator / denominator;
	gy = y_numerator / denominator;
	gz = polygenPoints[0].vetx.z;
	gcoordinate = Vec3(gx, gy, gz);
	return gcoordinate;
}
//���������
void QuadPart::OrderPolygens()
{
	int flag;
	double tempDistance;
	Vec3 outGravityCoord, tempGravityCoord, tempVec;  //����������
	varray<Vec3> polygensGravityCoord;
	varray<double> gravityDistance;
	varray<int> order;
	varray<MyPolygon> tempPolygens;

	outGravityCoord = CalGravity(out_Points);
	int polygensLength = large_Polygon.size();

	for (int i = 0; i < polygensLength; i++)
	{
		tempGravityCoord = CalGravity(large_Polygon[i].p_Points);
		polygensGravityCoord.push_back(tempGravityCoord);
	}

	for (int i = 0; i < polygensGravityCoord.size(); i++)
	{
		if (preProcessQuadPol) {
			if (large_Polygon[i].p_Points.size() != 4)
			{
				assert(large_Polygon[i].p_Points.size() > 2);
				tempVec = outGravityCoord - polygensGravityCoord[i];
				tempDistance = tempVec.Magnitude();
				gravityDistance.push_back(tempDistance);
			}
			else if (large_Polygon[i].p_Points.size() == 4)
			{
				gravityDistance.push_back(-1);
				quad_Polygon.push_back(large_Polygon[i]); //���ı������ݵ����ı��μ���
			}
		}
		else {
			assert(large_Polygon[i].p_Points.size() > 2);
			tempVec = outGravityCoord - polygensGravityCoord[i];
			tempDistance = tempVec.Magnitude();
			gravityDistance.push_back(tempDistance);
		}
	}

	while (order.size() < gravityDistance.size()) //������θ������Ľ�������
	{
		flag = 0;
		tempDistance = gravityDistance[0];
		for (int i = 0; i < gravityDistance.size(); i++)
		{
			if ((gravityDistance[i] > -1) && (gravityDistance[i] >= tempDistance))
			{
				tempDistance = gravityDistance[i];
				flag = i;
			}

			else if ((gravityDistance[i] == -1) && (gravityDistance[i] >= tempDistance))
			{
				tempDistance = gravityDistance[i];
				flag = i;

			}

			if (i == gravityDistance.size() - 1)
			{
				order.push_back(flag);
				gravityDistance[flag] = -2;
			}
		}
	}

	for (int i = 0; i < order.size(); i++)
	{
		if (preProcessQuadPol) {
			if (large_Polygon[order[i]].p_Points.size() != 4)
			{
				assert(large_Polygon[i].p_Points.size() > 2);
				tempPolygens.push_back(large_Polygon[order[i]]);
			}
			else if (large_Polygon[order[i]].p_Points.size() == 4)
			{
				continue;
			}
		}
		else {
			assert(large_Polygon[i].p_Points.size() > 2);
			tempPolygens.push_back(large_Polygon[order[i]]);
		}
	}

	large_Polygon.clear();
	large_Polygon = tempPolygens; //���ô����ıߵĶ���μ���
	if (priorityQuadPol) {
		//���ı��η��������ʷ֣�������ǰ
		tempPolygens.clear();
		tempPolygens.resize(large_Polygon.size());
		int i = 0;
		for (const auto& pol : large_Polygon) {
			if (pol.p_Points.size() == 4) {
				tempPolygens[i] = pol;
				++i;
			}
		}
		for (const auto& pol : large_Polygon) {
			if (pol.p_Points.size() == 4) {
				continue;
			}
			tempPolygens[i] = pol;
			++i;
		}
		assert(i == tempPolygens.size());
		large_Polygon.clear();
		large_Polygon = tempPolygens; //���ö���μ���
	}
}

void QuadPart::ResetPolygens()
{
	int largeLength = large_Polygon.size();
	int quadLength = quad_Polygon.size();
	for (int i = 0; i < largeLength; i++)
	{
		int linesNumber = large_Polygon[i].p_NurbsLines.size();
		for (int j = 0; j < linesNumber; j++)
		{
			Vec4 firstPoint = large_Polygon[i].p_NurbsLines[j].m_CtrlPts[0];
			Vec4 lastPoint = large_Polygon[i].p_NurbsLines[j].m_CtrlPts[large_Polygon[i].p_NurbsLines[j].m_CtrlPts.size() - 1];
			for (int k = 0; k < all_NurbsLines.size(); k++)
			{
				Vec4 tempFirstPoint = all_NurbsLines[k].m_CtrlPts[0];
				Vec4 tempLastPoint = all_NurbsLines[k].m_CtrlPts[all_NurbsLines[k].m_CtrlPts.size() - 1];
				if (JudgeTwoPointsCoincide(firstPoint, tempFirstPoint)
					&& JudgeTwoPointsCoincide(lastPoint, tempLastPoint))
				{
					large_Polygon[i].objectInAllLines.push_back(k);
				}
			}
		}
		large_Polygon[i].p_NurbsLines.clear(); //������߾������ݵ�����
	}

	for (int i = 0; i < quadLength; i++)
	{
		int linesNumber = quad_Polygon[i].p_NurbsLines.size();
		for (int j = 0; j < linesNumber; j++)
		{
			Vec4 firstPoint = quad_Polygon[i].p_NurbsLines[j].m_CtrlPts[0];
			Vec4 lastPoint = quad_Polygon[i].p_NurbsLines[j].m_CtrlPts[quad_Polygon[i].p_NurbsLines[j].m_CtrlPts.size() - 1];
			for (int k = 0; k < all_NurbsLines.size(); k++)
			{
				Vec4 tempFirstPoint = all_NurbsLines[k].m_CtrlPts[0];
				Vec4 tempLastPoint = all_NurbsLines[k].m_CtrlPts[all_NurbsLines[k].m_CtrlPts.size() - 1];
				if (JudgeTwoLinesCoincide(quad_Polygon[i].p_NurbsLines[j], all_NurbsLines[k]))
				{
					quad_Polygon[i].objectInAllLines.push_back(k);
					all_NurbsLines[k].ifSegLine = 0;
				}
			}
		}
		quad_Polygon[i].p_NurbsLines.clear(); //������߾������ݵ�����
	}
}



void QuadPart::FindVisiblePoint(MyPolygon &pendingPolygen, varray<VisiblePoints> &visiblePointPair)
{
	int pointsNumber = pendingPolygen.p_Points.size();
	VisiblePoints tempPointsPair;
	varray<VisiblePoints> pendingPointsPair; //���п��ܵĿ��ӵ�ԣ�������

	//��Ϊ����Σ�����(����ηŵ�����Ĺ����д���)
	if (pointsNumber == 5) 
	{
		return;
	}

	else if (pointsNumber > 5)
	{
		int flag1, flag2;
		for (int i = 0; i < pointsNumber - 3; i++)
		{
			//�������۵�һ���㣬�ڶ����㣬�Լ�֮��ĵ����ȡ��������
			if (i == 0)
			{
				flag1 = 3;
				flag2 = pointsNumber - 3;
				for (int j = flag1; j <= flag2; j++)
				{
					tempPointsPair.firstPoint = i;
					tempPointsPair.secondPoint = j;
					pendingPointsPair.push_back(tempPointsPair);
				}
			}
			else if (i == 1)
			{
				flag1 = 4;
				flag2 = pointsNumber - 2;
				for (int j = flag1; j <= flag2; j++)
				{
					tempPointsPair.firstPoint = i;
					tempPointsPair.secondPoint = j;
					pendingPointsPair.push_back(tempPointsPair);
				}
			}
			else
			{
				flag1 = i + 3;
				flag2 = pointsNumber - 1;
				for (int j = flag1; j <= flag2; j++)
				{
					tempPointsPair.firstPoint = i;
					tempPointsPair.secondPoint = j;
					pendingPointsPair.push_back(tempPointsPair);
				}
			}

			//��iȡ���������ĸ��㣬����ѭ��
			//if (i == pointsNumber - 4)
			//{
			//	break;
			//}
		}
	}

	//�����еĿ��ܿ��ӵ�Խ����ж�
	//�����ж��ǲ���ȫ͹
	bool allConcave = true;
	for (int i = 0; i < pointsNumber; i++)
	{
		if (pendingPolygen.convexHull[i] == 1)
		{
			allConcave = false;
			break;
		}
	}
	if (allConcave)
	{
		visiblePointPair = pendingPointsPair;
		return;
	}

	//ֱ��������ֱ�߽���
	int pairNumber = pendingPointsPair.size(); //����������Կ��ӵ����Ŀ

	if (cgalJudge) {
		Segment_2 *segments = new Segment_2[pointsNumber + 1];
		for (int i = 0; i < pointsNumber; i++)
		{
			if (i == (pointsNumber - 1))
			{
				segments[i] = Segment_2(Point_2(pendingPolygen.p_Points[i].vetx.x, pendingPolygen.p_Points[i].vetx.y),
					Point_2(pendingPolygen.p_Points[0].vetx.x, pendingPolygen.p_Points[0].vetx.y));
			}
			else
			{
				segments[i] = Segment_2(Point_2(pendingPolygen.p_Points[i].vetx.x, pendingPolygen.p_Points[i].vetx.y),
					Point_2(pendingPolygen.p_Points[i + 1].vetx.x, pendingPolygen.p_Points[i + 1].vetx.y));
			}
		}

		for (int i = 0; i < pairNumber; i++)
		{
			//�޶��������ߵĽǶ�����
			double angle1, angle2;
			Vec3 fvec, bvec, abvec;
			int flag = 1;
			abvec = pendingPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx
				- pendingPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx;
			GetVecs(pendingPolygen, pendingPointsPair[i].firstPoint, fvec, bvec);
			angle1 = GetAngle(bvec, fvec);
			angle2 = GetAngle(bvec, abvec);
			if (!((angle2 > 0) && (angle2 < angle1)))
			{
				flag = 0;
			}
			if (flag == 0)
			{
				continue;
			}

			abvec = -abvec;
			GetVecs(pendingPolygen, pendingPointsPair[i].secondPoint, fvec, bvec);
			angle1 = GetAngle(bvec, fvec);
			angle2 = GetAngle(bvec, abvec);
			if (!((angle2 > 0) && (angle2 < angle1)))
			{
				flag = 0;
			}
			if (flag == 0)
			{
				continue;
			}

			//�߶���
			segments[pointsNumber] = Segment_2(Point_2(pendingPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx.x, pendingPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx.y),
				Point_2(pendingPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx.x, pendingPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx.y));

			std::list<Point_2> pts;
			CGAL::compute_intersection_points(segments, segments + pointsNumber + 1,
				std::back_inserter(pts));

			if (pts.size() == 0)
			{
				tempPointsPair = pendingPointsPair[i];
				visiblePointPair.push_back(tempPointsPair);
			}
			else
			{
				pts.clear();
				continue;
			}
			pts.clear();
		}
		delete[] segments;
		segments = nullptr;
	}

	//visiblePointPair.clear();
	//CGAL��һ���ж�����
	else {
		//���ж�
		vector< Segment_22>segs;
		stringstream ss;
		for (int i = 0; i < pointsNumber; i++)
		{
			if (i == (pointsNumber - 1))
			{
				segs.push_back(Segment_22(Point_2(pendingPolygen.p_Points[i].vetx.x, pendingPolygen.p_Points[i].vetx.y),
					Point_2(pendingPolygen.p_Points[0].vetx.x, pendingPolygen.p_Points[0].vetx.y)));
			}
			else
			{
				segs.push_back(Segment_22(Point_2(pendingPolygen.p_Points[i].vetx.x, pendingPolygen.p_Points[i].vetx.y),
					Point_2(pendingPolygen.p_Points[i + 1].vetx.x, pendingPolygen.p_Points[i + 1].vetx.y)));
			}
		}

		for (int i = 0; i < pairNumber; i++)
		{
			//�޶��������ߵĽǶ�����
			double angle1, angle2;
			Vec3 fvec, bvec, abvec;
			int flag = 1;
			abvec = pendingPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx
				- pendingPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx;
			GetVecs(pendingPolygen, pendingPointsPair[i].firstPoint, fvec, bvec);
			angle1 = GetAngle(bvec, fvec);
			angle2 = GetAngle(bvec, abvec);
			if (!((angle2 > 0) && (angle2 < angle1)))
			{
				flag = 0;
			}
			if (flag == 0)
			{
				continue;
			}

			abvec = -abvec;
			GetVecs(pendingPolygen, pendingPointsPair[i].secondPoint, fvec, bvec);
			angle1 = GetAngle(bvec, fvec);
			angle2 = GetAngle(bvec, abvec);
			if (!((angle2 > 0) && (angle2 < angle1)))
			{
				flag = 0;
			}
			if (flag == 0)
			{
				continue;
			}


			Segment_22 line(Point_2(pendingPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx.x, pendingPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx.y),
				Point_2(pendingPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx.x, pendingPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx.y));
			bool isinter = false;//Ĭ�ϲ��ཻ
			for (int j = 0; j < pointsNumber; j++) {
				CGAL::cpp11::result_of<Intersect_2(Segment_22, Segment_22)>::type
					result = intersection(segs[j], line);
				if (!result) {
					continue;
				}
				else {
					if (const Segment_22* s = boost::get<Segment_22>(&*result)) {
						isinter = true;
						break;
					}
					else {
						const Point_2* p = boost::get<Point_2 >(&*result);
						ss.clear();
						double x, y;
						ss << p->x();
						ss >> x;
						ss.clear();
						ss << p->y();
						ss >> y;
						Vec3 point(x, y, 0);
						if (JudgeTwoPointsCoincide(point, pendingPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx)
							|| JudgeTwoPointsCoincide(point, pendingPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx)) {
							continue;
						}
						isinter = true;
						break;
					}
				}
			}
			if (!isinter) {
				//���ཻ
				tempPointsPair = pendingPointsPair[i];
				visiblePointPair.push_back(tempPointsPair);
			}
		}
	}	
	return;
}

void QuadPart::FindAnotherVp(const MyPolygon & inputPolygen, varray<VisiblePoints>& anotherVp)
{
	int pointsNumber = inputPolygen.p_Points.size();
	VisiblePoints tempPointsPair;
	varray<VisiblePoints> pendingPointsPair; //���п��ܵĿ��ӵ�ԣ�������

		int flag1, flag2;
		for (int i = 0; i < pointsNumber - 2; i++)
		{
			//�������۵�һ���㣬�ڶ����㣬�Լ�֮��ĵ����ȡ��������
			if (i == 0)
			{
				flag1 = 2;
				flag2 = pointsNumber - 2;
				for (int j = flag1; j <= flag2; j++)
				{
					tempPointsPair.firstPoint = i;
					tempPointsPair.secondPoint = j;
					pendingPointsPair.push_back(tempPointsPair);
				}
			}
			else if (i == 1)
			{
				flag1 = 3;
				flag2 = pointsNumber - 1;
				for (int j = flag1; j <= flag2; j++)
				{
					tempPointsPair.firstPoint = i;
					tempPointsPair.secondPoint = j;
					pendingPointsPair.push_back(tempPointsPair);
				}
			}
			else
			{
				flag1 = i + 2;
				flag2 = pointsNumber - 1;
				for (int j = flag1; j <= flag2; j++)
				{
					tempPointsPair.firstPoint = i;
					tempPointsPair.secondPoint = j;
					pendingPointsPair.push_back(tempPointsPair);
				}
			}

			//��iȡ���������ĸ��㣬����ѭ��
			//if (i == pointsNumber - 4)
			//{
			//	break;
			//}
		}

	//�����еĿ��ܿ��ӵ�Խ����ж�
		bool allConcave = true;
		for (int i = 0; i < pointsNumber; i++)
		{
			if (inputPolygen.convexHull[i] == 1)
			{
				allConcave = false;
				break;
			}
		}
		if (allConcave)
		{
			anotherVp = pendingPointsPair;
			return;
		}


		//ֱ��������ֱ�߽���
		int pairNumber = pendingPointsPair.size(); //����������Կ��ӵ����Ŀ

		if (cgalJudge) {
			Segment_2 *segments = new Segment_2[pointsNumber + 1];
			for (int i = 0; i < pointsNumber; i++)
			{
				if (i == (pointsNumber - 1))
				{
					segments[i] = Segment_2(Point_2(inputPolygen.p_Points[i].vetx.x, inputPolygen.p_Points[i].vetx.y),
						Point_2(inputPolygen.p_Points[0].vetx.x, inputPolygen.p_Points[0].vetx.y));
				}
				else
				{
					segments[i] = Segment_2(Point_2(inputPolygen.p_Points[i].vetx.x, inputPolygen.p_Points[i].vetx.y),
						Point_2(inputPolygen.p_Points[i + 1].vetx.x, inputPolygen.p_Points[i + 1].vetx.y));
				}
			}

			for (int i = 0; i < pairNumber; i++)
			{
				//�޶��������ߵĽǶ�����
				double angle1, angle2;
				Vec3 fvec, bvec, abvec;
				int flag = 1;
				abvec = inputPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx
					- inputPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx;
				GetVecs(inputPolygen, pendingPointsPair[i].firstPoint, fvec, bvec);
				angle1 = GetAngle(bvec, fvec);
				angle2 = GetAngle(bvec, abvec);
				if (!((angle2 > 0) && (angle2 < angle1)))
				{
					flag = 0;
				}
				if (flag == 0)
				{
					continue;
				}

				abvec = -abvec;
				GetVecs(inputPolygen, pendingPointsPair[i].secondPoint, fvec, bvec);
				angle1 = GetAngle(bvec, fvec);
				angle2 = GetAngle(bvec, abvec);
				if (!((angle2 > 0) && (angle2 < angle1)))
				{
					flag = 0;
				}
				if (flag == 0)
				{
					continue;
				}


				segments[pointsNumber] = Segment_2(Point_2(inputPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx.x, inputPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx.y),
					Point_2(inputPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx.x, inputPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx.y));

				std::list<Point_2> pts;
				CGAL::compute_intersection_points(segments, segments + pointsNumber + 1,
					std::back_inserter(pts));

				if (pts.size() == 0)
				{
					tempPointsPair = pendingPointsPair[i];
					anotherVp.push_back(tempPointsPair);
				}
				else
				{
					pts.clear();
					continue;
				}
				pts.clear();
			}
			delete[] segments;
		}
		
		//CGAL��һ���ж�����
		else {
			//���ж�
			vector< Segment_22>segs;
			stringstream ss;
			for (int i = 0; i < pointsNumber; i++)
			{
				if (i == (pointsNumber - 1))
				{
					segs.push_back(Segment_22(Point_2(inputPolygen.p_Points[i].vetx.x, inputPolygen.p_Points[i].vetx.y),
						Point_2(inputPolygen.p_Points[0].vetx.x, inputPolygen.p_Points[0].vetx.y)));
				}
				else
				{
					segs.push_back(Segment_22(Point_2(inputPolygen.p_Points[i].vetx.x, inputPolygen.p_Points[i].vetx.y),
						Point_2(inputPolygen.p_Points[i + 1].vetx.x, inputPolygen.p_Points[i + 1].vetx.y)));
				}
			}

			for (int i = 0; i < pairNumber; i++)
			{
				//�޶��������ߵĽǶ�����
				double angle1, angle2;
				Vec3 fvec, bvec, abvec;
				int flag = 1;
				abvec = inputPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx
					- inputPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx;
				GetVecs(inputPolygen, pendingPointsPair[i].firstPoint, fvec, bvec);
				angle1 = GetAngle(bvec, fvec);
				angle2 = GetAngle(bvec, abvec);
				if (!((angle2 > 0) && (angle2 < angle1)))
				{
					flag = 0;
				}
				if (flag == 0)
				{
					continue;
				}

				abvec = -abvec;
				GetVecs(inputPolygen, pendingPointsPair[i].secondPoint, fvec, bvec);
				angle1 = GetAngle(bvec, fvec);
				angle2 = GetAngle(bvec, abvec);
				if (!((angle2 > 0) && (angle2 < angle1)))
				{
					flag = 0;
				}
				if (flag == 0)
				{
					continue;
				}


				Segment_22 line(Point_2(inputPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx.x, inputPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx.y),
					Point_2(inputPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx.x, inputPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx.y));
				bool isinter = false;//Ĭ�ϲ��ཻ
				for (int j = 0; j < pointsNumber; j++) {
					CGAL::cpp11::result_of<Intersect_2(Segment_22, Segment_22)>::type
						result = intersection(segs[j], line);
					if (!result) {
						continue;
					}
					else {
						if (const Segment_22* s = boost::get<Segment_22>(&*result)) {
							isinter = true;
							break;
						}
						else {
							const Point_2* p = boost::get<Point_2 >(&*result);
							ss.clear();
							double x, y;
							ss << p->x();
							ss >> x;
							ss.clear();
							ss << p->y();
							ss >> y;
							Vec3 point(x, y, 0);
							if (JudgeTwoPointsCoincide(point, inputPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx)
								|| JudgeTwoPointsCoincide(point, inputPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx)) {
								continue;
							}
							isinter = true;
							break;
						}
					}
				}
				if (!isinter) {
					//���ཻ
					tempPointsPair = pendingPointsPair[i];
					anotherVp.push_back(tempPointsPair);
				}
			}
		}

		return;

		//������
	//int pairNumber = pendingPointsPair.size(); //����������Կ��ӵ����Ŀ
	//double A, B, C;
	//for (int i = 0; i < pairNumber; i++)
	//{
	//	int flag = 1;
	//	int leftLength, rightlengt; //���ҵ������
	//	varray<int> leftCollect, rightCollect;

	//	//��ȡλ�������������ߵĵ�
	//	if (pendingPointsPair[i].firstPoint == 0)
	//	{
	//		for (int j = pendingPointsPair[i].secondPoint + 1; j < pointsNumber; j++)
	//		{
	//			leftCollect.push_back(j);
	//		}
	//		for (int j = 1; j < pendingPointsPair[i].secondPoint; j++)
	//		{
	//			rightCollect.push_back(j);
	//		}
	//	}
	//	else
	//	{
	//		for (int j = pendingPointsPair[i].secondPoint + 1; j < pointsNumber; j++)
	//		{
	//			leftCollect.push_back(j);
	//		}
	//		for (int j = 0; j < pendingPointsPair[i].firstPoint; j++)
	//		{
	//			leftCollect.push_back(j);
	//		}
	//		for (int j = pendingPointsPair[i].firstPoint + 1; j < pendingPointsPair[i].secondPoint; j++)
	//		{
	//			rightCollect.push_back(j);
	//		}
	//	}

	//	//�趨����ֱ�߷��̺�������ABC
	//	A = inputPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx.y - inputPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx.y;
	//	B = inputPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx.x - inputPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx.x;
	//	C = inputPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx.y*inputPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx.x
	//		- inputPolygen.p_Points[pendingPointsPair[i].firstPoint].vetx.x*inputPolygen.p_Points[pendingPointsPair[i].secondPoint].vetx.y;

	//	//�����������������к������㣬�����о����������򽫵�Դ����赼���Ŀ��ӵ�Լ�����
	//	leftLength = leftCollect.size();
	//	rightlengt = rightCollect.size();

	//	for (int k = 0; k < leftLength; k++)
	//	{
	//		double Fxy = inputPolygen.p_Points[leftCollect[k]].vetx.x*A + inputPolygen.p_Points[leftCollect[k]].vetx.y*B + C;
	//		if (Fxy >= 0)
	//		{
	//			flag = 0;
	//			break;
	//		}
	//	}
	//	if (flag == 0)
	//	{
	//		continue;
	//	}

	//	for (int k = 0; k < rightlengt; k++)
	//	{
	//		double Fxy = inputPolygen.p_Points[rightCollect[k]].vetx.x*A + inputPolygen.p_Points[rightCollect[k]].vetx.y*B + C;
	//		if (Fxy <= 0)
	//		{
	//			flag = 0;
	//			break;
	//		}
	//	}
	//	if (flag == 0)
	//	{
	//		continue;
	//	}
	//	else
	//	{
	//		tempPointsPair = pendingPointsPair[i];
	//		anotherVp.push_back(tempPointsPair);
	//	}
	//}
}

//�ж��Ƿ���ڰ���
void QuadPart::JudgeConcavePoint(MyPolygon &pendingPolygen)
{
		int pointsNumber = pendingPolygen.p_Points.size();
		for (int i = 0; i < pointsNumber; i++)
		{
			Vec3 fvec, bvec;
			double angle;
			GetVecs(pendingPolygen, i, fvec, bvec);
			angle = GetAngle(bvec, fvec);
			if (angle > 0 && abs(PI-angle)>1e-4 && (angle-PI)<1e-4)
			{
				pendingPolygen.convexHull.push_back(0);
			}
			else
			{
				pendingPolygen.convexHull.push_back(1);
			}
		}

		//for (int j = 0; j < pointsNumber; j++)
		//{
		//	if (j == 0)  //�׵��ж�
		//	{
		//		Vec3 beginVex = pendingPolygen.p_Points[0].vetx - pendingPolygen.p_Points[pointsNumber - 1].vetx;
		//		Vec3 endVex = pendingPolygen.p_Points[1].vetx - pendingPolygen.p_Points[0].vetx;
		//		beginVex = beginVex.Normalize();
		//		endVex = endVex.Normalize();
		//		Vec3 cross = beginVex.Cross(endVex);
		//		
		//		if (cross.z > 0)
		//		{
		//			pendingPolygen.convexHull.push_back(0);
		//		}
		//		else
		//		{
		//			pendingPolygen.convexHull.push_back(1);
		//		}
		//	}

		//	else if (j == pointsNumber - 1)  //β���ж�
		//	{
		//		Vec3 beginVex = pendingPolygen.p_Points[pointsNumber - 1].vetx - pendingPolygen.p_Points[pointsNumber - 2].vetx;
		//		Vec3 endVex = pendingPolygen.p_Points[0].vetx - pendingPolygen.p_Points[pointsNumber - 1].vetx;
		//		beginVex = beginVex.Normalize();
		//		endVex = endVex.Normalize();
		//		Vec3 cross = beginVex.Cross(endVex);
		//		if (cross.z > 0)
		//		{
		//			pendingPolygen.convexHull.push_back(0);
		//		}
		//		else
		//		{
		//			pendingPolygen.convexHull.push_back(1);
		//		}
		//	}

		//	else  //������ж�
		//	{
		//		Vec3 beginVex = pendingPolygen.p_Points[j].vetx - pendingPolygen.p_Points[j-1].vetx;
		//		Vec3 endVex = pendingPolygen.p_Points[j+1].vetx - pendingPolygen.p_Points[j].vetx;
		//		beginVex = beginVex.Normalize();
		//		endVex = endVex.Normalize();
		//		Vec3 cross = beginVex.Cross(endVex);
		//		if (cross.z > 0)
		//		{
		//			pendingPolygen.convexHull.push_back(0);
		//		}
		//		else
		//		{
		//			pendingPolygen.convexHull.push_back(1);
		//		}
		//	}
		//}
}

void QuadPart::VisiblePointConcaveJudge(MyPolygon &inputPolygen, varray<VisiblePoints>& visiblePointPair)
{
	if(visiblePointPair.size()!=0)
	{
		for (int i = 0; i < visiblePointPair.size(); i++)
		{
			int flag1 = inputPolygen.convexHull[visiblePointPair[i].firstPoint];
			int flag2 = inputPolygen.convexHull[visiblePointPair[i].secondPoint];
			if (flag1 == 1 || flag2 == 1)
			{
				visiblePointPair[i].ifConvex = 1;
			}
			else
			{
				visiblePointPair[i].ifConvex = 0;
			}
		}
	}
	else
	{
		return;
	}	
}

void QuadPart::VisiblePointSharpJudeg(MyPolygon &inputPolygen, varray<VisiblePoints>& visiblePointPair, const varray<PartNurbsLine> &allLines)
{
	if (visiblePointPair.size() != 0)
	{
		//��AB����֮���������(����ֱ���Լ����ߵ�������)
		for (int i = 0; i < visiblePointPair.size(); i++)
		{
			int firstP = visiblePointPair[i].firstPoint;
			int lastP = visiblePointPair[i].secondPoint;
			int pNumber = inputPolygen.p_Points.size();
			Vec3 vecAB = inputPolygen.p_Points[lastP].vetx - inputPolygen.p_Points[firstP].vetx;
			Vec3 vecBA = inputPolygen.p_Points[firstP].vetx - inputPolygen.p_Points[lastP].vetx;
			Vec3 vecA0, vecA00, vecA1, vecA11, vecB0, vecB00, vecB1, vecB11;
			if (firstP == 0) //�׵�Ϊ0�����(����֮ǰ�Ĳ�����A�㲻�����Ƕ���ε��յ�)
			{
				//A��
				vecA0 = inputPolygen.p_Points[pNumber - 1].vetx - inputPolygen.p_Points[firstP].vetx;
				vecA1 = inputPolygen.p_Points[firstP + 1].vetx - inputPolygen.p_Points[firstP].vetx;
				
				int objectLine = inputPolygen.objectInAllLines[pNumber - 1];
				Vec3 temp = Vec3(allLines[objectLine].m_CtrlPts[0].x, allLines[objectLine].m_CtrlPts[0].y, allLines[objectLine].m_CtrlPts[0].z);
				if (inputPolygen.p_Points[firstP].vetx == temp)
				{
					vecA00 = allLines[objectLine].CalBeginDirecVec();	
				}
				else
				{
					vecA00 = allLines[objectLine].CalEndDirecVec();		
				}

				objectLine = inputPolygen.objectInAllLines[firstP];
				temp = Vec3(allLines[objectLine].m_CtrlPts[0].x, allLines[objectLine].m_CtrlPts[0].y, allLines[objectLine].m_CtrlPts[0].z);
				if (inputPolygen.p_Points[firstP].vetx == temp)
				{
					vecA11 = allLines[objectLine].CalBeginDirecVec();
				}
				else
				{
					vecA11 = allLines[objectLine].CalEndDirecVec();
				}

				//B��
				vecB0 = inputPolygen.p_Points[lastP - 1].vetx - inputPolygen.p_Points[lastP].vetx;
				vecB1 = inputPolygen.p_Points[lastP + 1].vetx - inputPolygen.p_Points[lastP].vetx;
				objectLine = inputPolygen.objectInAllLines[lastP - 1];
				temp = Vec3(allLines[objectLine].m_CtrlPts[0].x, allLines[objectLine].m_CtrlPts[0].y, allLines[objectLine].m_CtrlPts[0].z);
				if (inputPolygen.p_Points[lastP].vetx == temp)
				{
					vecB00 = allLines[objectLine].CalBeginDirecVec();
				}
				else
				{
					vecB00 = allLines[objectLine].CalEndDirecVec();
				}

				objectLine = inputPolygen.objectInAllLines[lastP];
				temp = Vec3(allLines[objectLine].m_CtrlPts[0].x, allLines[objectLine].m_CtrlPts[0].y, allLines[objectLine].m_CtrlPts[0].z);
				if (inputPolygen.p_Points[lastP].vetx == temp)
				{
					vecB11 = allLines[objectLine].CalBeginDirecVec();
				}
				else
				{
					vecB11 = allLines[objectLine].CalEndDirecVec();
				}
			}

			else //ͬ����B�㲻�����Ƕ���ε���ʼ��,���ǿ���Ϊ��ֹ��
			{
				//A��
				vecA0 = inputPolygen.p_Points[firstP - 1].vetx - inputPolygen.p_Points[firstP].vetx;
				vecA1 = inputPolygen.p_Points[firstP + 1].vetx - inputPolygen.p_Points[firstP].vetx;

				int objectLine = inputPolygen.objectInAllLines[firstP-1];
				Vec3 temp = Vec3(allLines[objectLine].m_CtrlPts[0].x, allLines[objectLine].m_CtrlPts[0].y, allLines[objectLine].m_CtrlPts[0].z);
				if (inputPolygen.p_Points[firstP].vetx == temp)
				{
					vecA00 = allLines[objectLine].CalBeginDirecVec();
				}
				else
				{
					vecA00 = allLines[objectLine].CalEndDirecVec();
				}

				objectLine = inputPolygen.objectInAllLines[firstP];
				temp = Vec3(allLines[objectLine].m_CtrlPts[0].x, allLines[objectLine].m_CtrlPts[0].y, allLines[objectLine].m_CtrlPts[0].z);
				if (inputPolygen.p_Points[firstP].vetx == temp)
				{
					vecA11 = allLines[objectLine].CalBeginDirecVec();
				}
				else
				{
					vecA11 = allLines[objectLine].CalEndDirecVec();
				}

				//B��,�����
				if (lastP == (pNumber - 1))
				{
					vecB0 = inputPolygen.p_Points[lastP - 1].vetx - inputPolygen.p_Points[lastP].vetx;
					vecB1 = inputPolygen.p_Points[0].vetx - inputPolygen.p_Points[lastP].vetx;
					objectLine = inputPolygen.objectInAllLines[lastP - 1];
					temp = Vec3(allLines[objectLine].m_CtrlPts[0].x, allLines[objectLine].m_CtrlPts[0].y, allLines[objectLine].m_CtrlPts[0].z);
					if (inputPolygen.p_Points[lastP].vetx == temp)
					{
						vecB00 = allLines[objectLine].CalBeginDirecVec();
					}
					else
					{
						vecB00 = allLines[objectLine].CalEndDirecVec();
					}

					objectLine = inputPolygen.objectInAllLines[lastP];
					temp = Vec3(allLines[objectLine].m_CtrlPts[0].x, allLines[objectLine].m_CtrlPts[0].y, allLines[objectLine].m_CtrlPts[0].z);
					if (inputPolygen.p_Points[lastP].vetx == temp)
					{
						vecB11 = allLines[objectLine].CalBeginDirecVec();
					}
					else
					{
						vecB11 = allLines[objectLine].CalEndDirecVec();
					}
				}
				else
				{
					vecB0 = inputPolygen.p_Points[lastP - 1].vetx - inputPolygen.p_Points[lastP].vetx;
					vecB1 = inputPolygen.p_Points[lastP + 1].vetx - inputPolygen.p_Points[lastP].vetx;
					objectLine = inputPolygen.objectInAllLines[lastP - 1];
					temp = Vec3(allLines[objectLine].m_CtrlPts[0].x, allLines[objectLine].m_CtrlPts[0].y, allLines[objectLine].m_CtrlPts[0].z);
					if (inputPolygen.p_Points[lastP].vetx == temp)
					{
						vecB00 = allLines[objectLine].CalBeginDirecVec();
					}
					else
					{
						vecB00 = allLines[objectLine].CalEndDirecVec();
					}

					objectLine = inputPolygen.objectInAllLines[lastP];
					temp = Vec3(allLines[objectLine].m_CtrlPts[0].x, allLines[objectLine].m_CtrlPts[0].y, allLines[objectLine].m_CtrlPts[0].z);
					if (inputPolygen.p_Points[lastP].vetx == temp)
					{
						vecB11 = allLines[objectLine].CalBeginDirecVec();
					}
					else
					{
						vecB11 = allLines[objectLine].CalEndDirecVec();
					}
				}
			}	

			//�����������AB��cosֵ,��A��B
			visiblePointPair[i].ifSharp = 0;
			double cosA0 = (vecAB.Dot(vecA0)) / (vecAB.Magnitude() * vecA0.Magnitude());
			double cosA1 = (vecAB.Dot(vecA1)) / (vecAB.Magnitude() * vecA1.Magnitude());
			double cosA00 = (vecAB.Dot(vecA00)) / (vecAB.Magnitude() * vecA00.Magnitude());
			double cosA11 = (vecAB.Dot(vecA11)) / (vecAB.Magnitude() * vecA11.Magnitude());
			if ((cosA0 < 1 && cosA0 > sharpAngle) && (cosA00 < 1 && cosA00 > sharpAngle))
			{
				visiblePointPair[i].ifSharp = 1;
				continue;
			}
			if ((cosA1 < 1 && cosA1 > sharpAngle) && (cosA11 < 1 && cosA11 > sharpAngle))
			{
				visiblePointPair[i].ifSharp = 1;
				continue;
			}
			double cosB0 = (vecBA.Dot(vecB0)) / (vecBA.Magnitude() * vecB0.Magnitude());
			double cosB1 = (vecBA.Dot(vecB1)) / (vecBA.Magnitude() * vecB1.Magnitude());
			double cosB00 = (vecBA.Dot(vecB00)) / (vecBA.Magnitude() * vecB00.Magnitude());
			double cosB11 = (vecBA.Dot(vecB11)) / (vecBA.Magnitude() * vecB11.Magnitude());
			if ((cosB0 < 1 && cosB0 > sharpAngle) && (cosB00 < 1 && cosB00 > sharpAngle))
			{
				visiblePointPair[i].ifSharp = 1;
				continue;
			}
			if ((cosB1 < 1 && cosB1 > sharpAngle) && (cosB11 < 1 && cosB11 > sharpAngle))
			{
				visiblePointPair[i].ifSharp = 1;
				continue;
			}
		}
	}

	else
	{
		return;
	}
}

//�Ƿ���ڼ���ж�
bool QuadPart::LinkLineSharpJudge(const Vec3 &frontVec, const Vec3 &lastVec, const Vec3 &pt, const Vec3 &pu)
{
	Vec3 vecLeft = frontVec;
	Vec3 vecRight = lastVec;
	Vec3 vecU = pu - pt;
	vecLeft = vecLeft.Normalize();
	vecRight = vecRight.Normalize();
	vecU = vecU.Normalize();

	double angle1 = GetAngle(vecRight, vecU);
	double angle2 = GetAngle(vecU, vecLeft);
	if (angle1 > 0 && angle1 < Sharp)
	{
		return false;
	}
	if(angle2 > 0 && angle2 < Sharp)
	{
		return false;
	}
	return true;
}

void QuadPart::GetFBpNumber(const MyPolygon & inputPolygen, int pointNumber, int & frontP, int & behindP)
{
	int pLength = inputPolygen.p_Points.size();
	if (pointNumber == 0)
	{
		frontP = pLength - 1;
		behindP = pointNumber + 1;
	}
	else if (pointNumber == pLength - 1)
	{
		frontP = pointNumber - 1;
		behindP =0;
	}
	else
	{
		frontP = pointNumber - 1;
		behindP = pointNumber + 1;
	}
}

void QuadPart::GetFBcoord(const MyPolygon & inputPolygen, int pointNumber, Vec3 & frontCoord, Vec3 & behindCoord)
{
	int pLength = inputPolygen.p_Points.size();
	if (pointNumber == 0)
	{
		frontCoord = inputPolygen.p_Points[pLength - 1].vetx;
		behindCoord = inputPolygen.p_Points[pointNumber + 1].vetx ;
	}
	else if (pointNumber == pLength - 1)
	{
		frontCoord = inputPolygen.p_Points[pointNumber - 1].vetx ;
		behindCoord = inputPolygen.p_Points[0].vetx;
	}
	else
	{
		frontCoord = inputPolygen.p_Points[pointNumber - 1].vetx;
		behindCoord = inputPolygen.p_Points[pointNumber + 1].vetx;
	}
}

void DecoincideNurbsLine(varray<Spline>& nl)
{
	assert(nl.size());
	for (int i = 0; i < nl.size();++i) {
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

void ReSetTreeNode(SfCtainTreeNode * father)
{
	if (!father->childs.empty())
	{
		list<SfCtainTreeNode*>::iterator left = father->childs.begin();
		list<SfCtainTreeNode*>::iterator right;
		//�����жϸ��ڵ��������ӽڵ�Ļ��������ϵ
		for (; left != (father->childs.end());)
		{
			right = left;
			bool flag1 = false, flag2 = false;
			for (++right; right != father->childs.end();)
			{
				//���ж�right�Ƿ���left�ڲ�
				flag1 = PolygonsRelationship((*right)->outLines, (*left)->outLines);
				if (flag1) {
					//��right��left�ڲ�
					list<SfCtainTreeNode*>::iterator tmpit = right;
					tmpit++;
					//��right����left���ӽڵ�list��
					(*left)->childs.push_back(*right);
					//ɾ��right
					father->childs.erase(right);
					right = tmpit;
					continue;
				}

				//��right����left�ڲ�,�ж�left�Ƿ���right�ڲ�
				flag2 = PolygonsRelationship((*left)->outLines, (*right)->outLines);
				if (flag2) {
					//��left��right�ڲ�
					list<SfCtainTreeNode*>::iterator tmpit = left;
					tmpit++;
					//��left����right�ӽڵ�list��
					(*right)->childs.push_back(*left);
					//ɾ��left
					father->childs.erase(left);
					left = tmpit;
					break;
				}
				right++;
			}
			if (!flag2)
				left++;
		}

		//��������������
		for (left = father->childs.begin(); left != father->childs.end(); left++)
		{
			if (father->conLines.size() == 0) break;//��û�и����ߵ����
			if ((*left)->childs.empty()) continue;//û���ӽڵ㣬���踨����
			varray<Spline> line;
			bool flag;
			varray<Spline>::iterator it = father->conLines.begin();
			for (; it != father->conLines.end();) {
				line.clear();
				line.push_back(*it);
				flag = PolygonsRelationship(line, (*left)->outLines);//�жϸ�������������ϵ
				if (flag) {
					//���������������ߴ����ӽڵ���
					(*left)->conLines.push_back(*it);
					//����erase���ܴ������⣡����
					father->conLines.erase(it);
				}
				else {
					it++;
				}
			}
		}
	}
}

SfCtainTreeNode * CreateSurfContainTree(const varray<varray<Spline>>& surf, const varray<Spline> conlines, varray<varray<int>> seg, varray<bool> genus)
{
	if (surf.size() == 0) return nullptr;
	//��genusΪ��,Ĭ�����ж��Ƿǿ����(�������涼Ҫִ���ʷ�)
	if (genus.size() == 0) genus.resize(surf.size(), false);
	//��segΪ��,Ĭ�����пɷָ�
	if (seg.size() == 0) {
		seg.resize(surf.size());
		for (int i = 0; i < seg.size(); ++i) {
			//ÿ�������е����߶���Ϊ�ɷָ�
			seg[i].resize(surf[i].size(), 1);
		}
	}

	assert(surf.size() == genus.size());//���� ���ʽΪ�٣���ӡ������Ϣ������abort��ֹ����

	queue< SfCtainTreeNode *> nodes;
	//��һ���������������Χ����,������������������Ϊroot��������������ߣ������ٽ��д���
	SfCtainTreeNode* root = new SfCtainTreeNode(surf[0], conlines);
	root->isSeg = seg[0];	//���Ͽɷָ�����־
	root->isGenus = genus[0];

	//���
	root->num = 1;

	SfCtainTreeNode* cur = nullptr;
	for (int i = 1; i < surf.size(); ++i)
	{
		cur = new SfCtainTreeNode(surf[i]);
		cur->isSeg = seg[i];
		cur->isGenus = genus[i];
		cur->num = i+1;//�趨���
		//������������������Ϊ���ڵ���ӽڵ�
		root->childs.push_back(cur);
		//����
		//cur->vertex = OrderLinesAntioclock(cur->outLines);
		//cur->isSort = true;
	}

	//��������������нڵ�
	nodes.push(root);
	while (!nodes.empty())
	{
		cur = nodes.front();
		nodes.pop();
		if (!cur->childs.empty()) {
			ReSetTreeNode(cur);
			//����cur�ڵ�������ӽڵ�
			list< SfCtainTreeNode *>::iterator it = cur->childs.begin();
			for (; it != cur->childs.end();) {
				nodes.push(*it);
				it++;
			}
		}
	}
	return root;//���ظ����
}

SfCtainTreeNode * CreateSurfContainTree(const varray<varray<Spline>>& surf, const varray<Spline> conlines, varray<varray<int>> seg, varray<bool> genus, varray<int> sfnum)
{
	if (surf.size() == 0) return nullptr;
	//��genusΪ��,Ĭ�����ж��Ƿǿ����(�������涼Ҫִ���ʷ�)
	if (genus.size() == 0) genus.resize(surf.size(), false);
	//��segΪ��,Ĭ�����пɷָ�
	if (seg.size() == 0) {
		seg.resize(surf.size());
		for (int i = 0; i < seg.size(); ++i) {
			seg[i].resize(surf[i].size(), 1);
		}
	}
	//�������趨
	bool numflag = false;
	if (sfnum.size() == 0)
		numflag = true;		//����Ĭ��1/2/3��˳���Զ����
	else if (sfnum[0] == 0)
		numflag = true;
	else
		numflag = false;


	assert(surf.size() == genus.size());

	queue< SfCtainTreeNode *> nodes;
	//��һ���������������Χ����,������������������Ϊroot��������������ߣ������ٽ��д���
	SfCtainTreeNode* root = new SfCtainTreeNode(surf[0], conlines);
	root->isSeg = seg[0];
	root->isGenus = genus[0];

	//���
	if (numflag)
		root->num = 1;
	else
		root->num = sfnum[0];
	


	SfCtainTreeNode* cur = nullptr;
	for (int i = 1; i < surf.size(); ++i)
	{
		cur = new SfCtainTreeNode(surf[i]);
		cur->isSeg = seg[i];
		cur->isGenus = genus[i];
		root->childs.push_back(cur);
		if (numflag) {
			//���
			cur->num = i + 1;
		}
		else {
			cur->num = sfnum[i];
		}
		//����
		//cur->vertex = OrderLinesAntioclock(cur->outLines);
		//cur->isSort = true;
	}

	//��������������нڵ�
	nodes.push(root);
	while (!nodes.empty())
	{
		cur = nodes.front();
		nodes.pop();
		if (!cur->childs.empty()) {
			ReSetTreeNode(cur);
			//����cur�ڵ�������ӽڵ�
			list< SfCtainTreeNode *>::iterator it = cur->childs.begin();
			for (; it != cur->childs.end();) {
				nodes.push(*it);
				it++;
			}
		}
	}
	return root;//���ظ����
}

void QuadTreeNodeSurf(SfCtainTreeNode * root)
{
	if (root->quaded) return;							//�Ѿ��ʷֹ����򷵻�
	varray<Spline> outLine = root->outLines;			//���������߼���
	varray<Spline> inLine;								//���������߼���
	varray<Spline> conLine = root->conLines;			//���������߼���
	list<SfCtainTreeNode*>::iterator it = root->childs.begin();

	//���ӽڵ���������Ϊ����������
	cout << "�����ӽڵ���Ϊ��������" << endl;
	for (; it != root->childs.end();) {
		//��������������
		SfCtainTreeNode * cur = (*it);
		for (int i = 0; i < cur->outNumber.size(); i++) {
			inLine.push_back(cur->allLines[cur->outNumber[i]]);
		}
		it++;
	}

	//�������������㼯���߼�
	QuadPart* quad = new QuadPart(outLine, inLine, conLine, root->isSeg);

	//���Ӳ��������εĵ����꼰��Ӧ����εĿ�����
	for (it = root->childs.begin(); it != root->childs.end();) {
		SfCtainTreeNode * cur = (*it);
		varray<Vec3> points;
		for (auto p : cur->vertex) {
			points.push_back(p);
		}
		quad->removepoints_Pol.push_back(points);
		quad->remove_genus.push_back(cur->isGenus);
		quad->remove_seg.push_back(cur->isSeg);
		it++;
	}

	if (root->isGenus) {
		//�ǿ�������,����outNumber����
		quad->OrderOutPointsAntioclock();
		varray<Vec3> points;
		for (auto p : quad->all_Points) {
			points.push_back(p.vetx);
		}
		//���Լ���Ϊ�������ʷֵĶ����
		cout << "�������棬�������ʷ�" << endl;
		quad->removepoints_Pol.push_back(points);
		quad->remove_genus.push_back(true);
		quad->remove_seg.push_back(root->isSeg);
		quad->CreateALGraph();
		quad->CreatePolygon();
		quad->OrderPolygens();
		quad->ResetPolygens();
		quad->RemovePolygon();
		quad->ResetData();		
	}
	else {
		//ִ���ʷֲ���
		cout << "�ʷ���..." << endl;
		quad->ExecuteQuadPart();
	}


	//��ȡ�����ʷ���ɵĽ��
	root->allLines.resize(quad->all_NurbsLines.size());
	for (int i = 0; i < root->allLines.size(); i++) {
		//��ȡ��������
		root->allLines[i] = quad->all_NurbsLines[i];
	}
	root->vertex.resize(quad->out_Points.size());
	for (int i = 0; i < root->vertex.size(); i++) {
		//��ȡ��������������
		root->vertex[i] = quad->out_Points[i].vetx;
	}
	root->outNumber.resize(quad->out_Number.size());
	for (int i = 0; i < root->outNumber.size(); i++) {
		//��ȡ�������������
		root->outNumber[i] = quad->out_Number[i];
	}
	root->quadPolNumber.resize(quad->quad_Polygon.size(), varray<int>(4));
	for (int i = 0; i < root->quadPolNumber.size(); i++) {
		//��ȡ�����ı��ε�����
		for (int j = 0; j < 4; j++) {
			root->quadPolNumber[i][j] = quad->quad_Polygon[i].objectInAllLines[j];
		}
	}
	//ȥ��������������
	root->outLines.clear();

	//�ͷ��ʷֶ���
	delete quad;
	quad = nullptr;
	root->quaded = true;		//���������ʷ�
}

void QuadWithContainTree(SfCtainTreeNode * root)
{
	cout << "��ʼ�ʷ�..." << endl;
	stack<SfCtainTreeNode *> nodes;
	queue<SfCtainTreeNode *> mq;
	mq.push(root);
	//���нڵ���ջ ����ջ�ײ�Ϊ���ڵ㣬����ΪҶ�ӽڵ�
	cout << "���нڵ����ջ��..." << endl;
	while (!mq.empty()) {
		SfCtainTreeNode * cur = mq.front();
		mq.pop();
		nodes.push(cur);
		//��Ϊ�����Σ���ǰ�����ʷ�
		if (cur->childs.empty() && !cur->isGenus && !cur->quaded&&cur->outLines.size() == 3) {
			cout << "���һ���ڵ�����Ϊ�����Σ���ǰ�ʷ�..." << endl;
			QuadTreeNodeSurf(cur);
			continue;
		}
		//�ӽڵ㲻Ϊ�գ����ӽڵ�������
		if (!cur->childs.empty())
		{
			list< SfCtainTreeNode *>::iterator it = cur->childs.begin();
			for (; it != cur->childs.end();)
			{
				mq.push(*it);
				it++;
			}
		}
	}

	//�����е����߷ŵ�al��
	varray<Spline> al;//lambda���ʽ
	cout << "�����нڵ�����߷ŵ�al��..." << endl;
	auto bfs = [&al,root]() {
		al.clear();
		queue<SfCtainTreeNode *> qq;
		if (root) qq.push(root);
		while (!qq.empty()) {
			SfCtainTreeNode *tmp = qq.front();
			qq.pop();
			for (auto i : tmp->allLines) {
				al.push_back(i);
			}
			if (!tmp->childs.empty()) {
				list< SfCtainTreeNode *>::iterator it;
				for (it = tmp->childs.begin(); it != tmp->childs.end(); it++) {
					qq.push(*it);
				}
			}
		}
	};

	//�����һ���ڵ㿪ʼ��������ʷ�
	int i = 0;
	cout << "�����һ���ڵ���������ʷ�..." << endl;
	while (!nodes.empty()) {
		SfCtainTreeNode * cur = nodes.top();
		nodes.pop();
		//if (cur->num == 1) 
		//	system("pause");
		QuadTreeNodeSurf(cur);
	}
}

void GetQuadSurfaces(SfCtainTreeNode * root, varray<SplineSurface>& surf)
{
	assert(root != nullptr);

}

void StartTreeOperation(const varray<varray<Spline>>& surf, const varray<Spline>& conlines, const varray<varray<int>> seg, varray<bool> genus, varray<Spline>& allLines, varray<SplineSurface>& allSurf)
{
	allLines.clear();
	allSurf.clear();
	SfCtainTreeNode* root = CreateSurfContainTree(surf, conlines, seg, genus);

	QuadWithContainTree(root);

	queue<SfCtainTreeNode*> q;
	q.push(root);

	while (!q.empty())
	{
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
		cur->GetSurfs(allSurf);

		list<SfCtainTreeNode*>::iterator it = cur->childs.begin();
		for (; it != cur->childs.end(); it++) {
			q.push(*it);
		}
	}
	delete root;
	root = nullptr;
}

Spline CreatNurbsLine(const Vec3 & p1, const Vec3 & p2)
{
	Spline tempLine;
	tempLine.m_Degree = 1;
	tempLine.m_Knots.push_back(0);
	tempLine.m_Knots.push_back(0);
	tempLine.m_Knots.push_back(1);
	tempLine.m_Knots.push_back(1);
	tempLine.m_CtrlPts.push_back(Vec4(p1));
	tempLine.m_CtrlPts.push_back(Vec4(p2));
	tempLine.DegreeElevate(2);
	return tempLine;
}

void GetLinesByPoint(const MyPolygon & inputPolygen, int pointNumber, const varray<PartNurbsLine>& allLine, Spline & frLine, Spline & backLine)
{
	int n = inputPolygen.p_Points.size();
	int frNum = (pointNumber - 1) >= 0 ? (pointNumber - 1) : (n - 1);
	frLine = allLine[inputPolygen.objectInAllLines[frNum]];
	backLine = allLine[inputPolygen.objectInAllLines[pointNumber]];
}

void GetVecsByPoint(const MyPolygon & inputPolygen, int pointNumber, const varray<PartNurbsLine>& allLine, Vec3 & frVec, Vec3 & backVec)
{
	Spline fl, bl;
	GetLinesByPoint(inputPolygen, pointNumber, allLine, fl, bl);

	Vec3 p = inputPolygen.p_Points[pointNumber].vetx;
	frVec = (JudgeTwoPointsCoincide(p, fl.m_CtrlPts[0])) ? ((Vec3)(fl.m_CtrlPts[1] - fl.m_CtrlPts[0])).Normalize() :
		((Vec3)(*(fl.m_CtrlPts.end() - 2) - *(fl.m_CtrlPts.end() - 1))).Normalize();
	backVec = (JudgeTwoPointsCoincide(p, bl.m_CtrlPts[0])) ? ((Vec3)(bl.m_CtrlPts[1] - bl.m_CtrlPts[0])).Normalize() :
		((Vec3)(*(bl.m_CtrlPts.end() - 2) - *(bl.m_CtrlPts.end() - 1))).Normalize();
}

void GetFBcoord(const MyPolygon & inputPolygen, int pointNumber, Vec3 & frontCoord, Vec3 & behindCoord)
{
	int pLength = inputPolygen.p_Points.size();
	if (pointNumber == 0)
	{
		frontCoord = inputPolygen.p_Points[pLength - 1].vetx;
		behindCoord = inputPolygen.p_Points[pointNumber + 1].vetx;
	}
	else if (pointNumber == pLength - 1)
	{
		frontCoord = inputPolygen.p_Points[pointNumber - 1].vetx;
		behindCoord = inputPolygen.p_Points[0].vetx;
	}
	else
	{
		frontCoord = inputPolygen.p_Points[pointNumber - 1].vetx;
		behindCoord = inputPolygen.p_Points[pointNumber + 1].vetx;
	}
}

//����������
PartNurbsLine GetPartNurbsLine(const Vec3 &p1, const Vec3 &p2, const varray<PartNurbsLine> &allLine, MyPolygon pol)
{
	PartNurbsLine line;
	Spline tempLine;
	Vec3 vec1, vec2;							//p1��p2�������ߵ�����
	int num1 = -1, num2 = -1;					//p1��p2�ڶ�����е����
	Vec3 midPoint = Vec3((p1.x + p2.x)*0.5, (p1.y + p2.y)*0.5, 1);//�����е�

	vec1 = (p2 - p1).Normalize();
	vec2 = (p1 - p2).Normalize();
	for (int i = 0; i < pol.p_Points.size(); i++) {
		if (JudgeTwoPointsCoincide(p1, pol.p_Points[i].vetx)) {
			num1 = i;
			break;
		}
	}
	for (int i = 0; i < pol.p_Points.size(); i++) {
		if (JudgeTwoPointsCoincide(p2, pol.p_Points[i].vetx)) {
			num2 = i;
			break;
		}
	}
	assert(num1 != -1 && num2 != -1);

	Vec3 frV1, bkV1, frV2, bkV2;				//p1��p2ǰ�����ߵ�������
	GetVecsByPoint(pol, num1, allLine, frV1, bkV1);
	GetVecsByPoint(pol, num2, allLine, frV2, bkV2);

	Vec3 fstrV1, bstrV1, fstrV2, bstrV2;		//p1��p2ǰ��ֱ�ߵķ�������
	GetFBcoord(pol, num1, fstrV1, bstrV1);
	fstrV1 = (fstrV1 - p1).Normalize();
	bstrV1 = (bstrV1 - p1).Normalize();

	GetFBcoord(pol, num2, fstrV2, bstrV2);
	fstrV2 = (fstrV2 - p2).Normalize();
	bstrV2 = (bstrV2 - p2).Normalize();

	double agle1, stragle1, agle2, stragle2;	//�����������ֱ�߼�����������֮��ĽǶ�
	double pagle1, pstragle1, pagle2, pstragle2;//p1��p2���Ժ�����(ֱ��)�������������ֱ�����������Ƕ�

	agle1 = GetAngle(bkV1, frV1);
	agle2 = GetAngle(bkV2, frV2);
	stragle1 = GetAngle(bstrV1, fstrV1);
	stragle2 = GetAngle(bstrV2, fstrV2);

	pagle1 = GetAngle(bkV1, vec1);
	pstragle1 = GetAngle(bstrV1, vec1);
	pagle2 = GetAngle(bkV2, vec2);
	pstragle2 = GetAngle(bstrV2, vec2);

	int flag = 0, flag2 = 0;
	//�ж�p1�Ƿ����Ҫ��
	if (pagle1 <= agle1) {
		//p1,p2�����������߼нǷ�Χ��
		if (abs(pagle1 - 0.0) < 1e-2 || (pagle1<(15/180.0*PI))) {
			//����������
			flag = 1;
		}
		else if (abs(pagle1 - agle1) < 1e-2 ||(agle1-pagle1)< (15 / 180.0 * PI)) {
			//����ǰ����
			flag = 2;
		}
	}
	else{
		//�����������߼нǷ�Χ��
		double ang1 = (2 * PI - agle1)*0.5;
		double ang2 = pagle1 - agle1;
		if (ang2 >= ang1) {
			//����������
			flag = 1;
		}
		else {
			//����ǰ����
			flag = 2;
		}
	}

	if (flag == 0) {
		//p1�����Ҫ��,�ж�p2��
		if (pagle2 <= agle2) {
			//p2,p1�����������߼нǷ�Χ��
			if (abs(pagle2 - 0.0) < 1e-2 || (pagle2 < (15 / 180.0 * PI))) {
				//����������
				flag2 = 1;
			}
			else if (abs(pagle2 - agle2) < 1e-2 || (agle2 - pagle2) < (15 / 180.0 * PI)) {
				//����ǰ����
				flag2 = 2;
			}
		}
		else {
			//�����������߼нǷ�Χ��
			double ang1 = (2 * PI - agle2)*0.5;
			double ang2 = pagle2 - agle2;
			if (ang2 >= ang1) {
				//����������
				flag2 = 1;
			}
			else {
				//����ǰ����
				flag2 = 2;
			}
		}

		if (flag2 == 1) {
			//p2->p1���������ߣ�p2���м俿����p1��Ҫ������߿���
			double angle = GetAngle(bkV2, fstrV2);//����ǰֱ�����Ƕ�
			angle = (angle < agle2) ? angle : agle2;//ѡȡ��С�ĽǶ�
			angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
			Vec3 v2 = bkV2.RotateZ(angle).Normalize();//����p2�����߷�������
			
			//����p1�����߷�������
			angle = (pagle1 < pstragle1) ? pagle1 : pstragle1;
			angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
			Vec3 v1 = vec1.RotateZ(-angle).Normalize();

			double length = max(abs(p1.x - p2.x), abs(p1.y - p2.y));
			Vec3 tmp1 = p1 + v1 * 3 * length;
			Vec3 tmp2 = p2 + v2 * 3 * length;
			midPoint = GetTwoLinesCrossPoint(p1, tmp1, p2, tmp2);
		}

		else if (flag2 == 2) {
			//p2->p1����ǰ���ߣ�p2���м俿����p1��ǰ���߿���
			double angle = GetAngle(bstrV2, frV2);//��ֱǰ�������Ƕ�
			angle = (angle < agle2) ? angle : agle2;
			angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
			Vec3 v2 = frV2.RotateZ(-angle).Normalize();//����p2�����߷�������

			//����p1�����߷�������
			double ag1 = agle1 - pagle1;//�������ߵ�ǰ�����������ĽǶ�
			double strag1 = stragle1 - pstragle1;//�������ߵ�ǰֱ�������ĽǶ�
			assert(ag1 >= 0);
			assert(strag1 >= 0);
			angle = (ag1 < strag1) ? ag1 : strag1;
			angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
			Vec3 v1 = vec1.RotateZ(angle).Normalize();//��ʱ����ת

			double length = max(abs(p1.x - p2.x), abs(p1.y - p2.y));
			Vec3 tmp1 = p1 + v1 * 3 * length;
			Vec3 tmp2 = p2 + v2 * 3 * length;
			midPoint = GetTwoLinesCrossPoint(p1, tmp1, p2, tmp2);
		}
	}

	else if (flag == 1) {
		//p1->p2���������ߣ�p1���м俿����p2��Ҫ������߿���
		double angle = GetAngle(bkV1, fstrV1);//����ǰֱ�����Ƕ�
		angle = (angle < agle1) ? angle : agle1;//ѡȡ��С�ĽǶ�
		angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
		Vec3 v1 = bkV1.RotateZ(angle).Normalize();//����p1�����߷�������

		//����p2�����߷�������
		angle = (pagle2 < pstragle2) ? pagle2 : pstragle2;
		angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
		Vec3 v2 = vec2.RotateZ(-angle).Normalize();

		double length = max(abs(p1.x - p2.x), abs(p1.y - p2.y));
		Vec3 tmp1 = p1 + v1 * 3 * length;
		Vec3 tmp2 = p2 + v2 * 3 * length;
		midPoint = GetTwoLinesCrossPoint(p1, tmp1, p2, tmp2);
	}

	else if (flag == 2) {
		//p1->p2����ǰ���ߣ�p1���м俿����p2��Ҫ��ǰ���߿���
		double angle = GetAngle(bstrV1, frV1);//��ֱǰ�������Ƕ�
		angle = (angle < agle1) ? angle : agle1;
		angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
		Vec3 v1 = frV1.RotateZ(-angle).Normalize();//����p2�����߷�������

		//����p2�����߷�������
		double ag2 = agle2 - pagle2;//�������ߵ�ǰ�����������ĽǶ�
		double strag2 = stragle2 - pstragle2;//�������ߵ�ǰֱ�������ĽǶ�
		//assert(ag2 >= 0);
		if (ag2 < 0) {
			//p2��ǰ��������������ͬһ���ˣ�ֱ�Ӳ�Ҫ�������
			return {};
		}
		assert(strag2 >= 0);
		angle = (ag2 < strag2) ? ag2 : strag2;
		angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
		Vec3 v2 = vec2.RotateZ(angle).Normalize();//��ʱ����ת

		double length = max(abs(p1.x - p2.x), abs(p1.y - p2.y));
		Vec3 tmp1 = p1 + v1 * 3 * length;
		Vec3 tmp2 = p2 + v2 * 3 * length;
		midPoint = GetTwoLinesCrossPoint(p1, tmp1, p2, tmp2);
	}

	if (flag == 0 && flag2 == 0) {
		//ֱ�Ӳ���ֱ����ʽ��NURBS��������
		tempLine.m_Degree = 1;
		tempLine.m_Knots.push_back(0);
		tempLine.m_Knots.push_back(0);
		tempLine.m_Knots.push_back(1);
		tempLine.m_Knots.push_back(1);
		tempLine.m_CtrlPts.push_back(Vec4(p1));
		tempLine.m_CtrlPts.push_back(Vec4(p2));
		tempLine.DegreeElevate(2);
	}
	else {
		//����������ʽ��NURBS�������� ��������
		tempLine.m_Degree = 2;
		tempLine.m_Knots.push_back(0);
		tempLine.m_Knots.push_back(0);
		tempLine.m_Knots.push_back(0);
		tempLine.m_Knots.push_back(1);
		tempLine.m_Knots.push_back(1);
		tempLine.m_Knots.push_back(1);
		tempLine.m_CtrlPts.push_back(Vec4(p1));
		tempLine.m_CtrlPts.push_back(Vec4(midPoint.x, midPoint.y, midPoint.z, 1));
		tempLine.m_CtrlPts.push_back(Vec4(p2));
	}

	//�ж��½������Ƿ�Ͷ�������������н���
	int fnum1 = pol.GetFrontNumber(num1);
	int fnum2 = pol.GetFrontNumber(num2);
	for (int i = 0; i < pol.objectInAllLines.size(); ++i) {
		if (i != num1 && i != num2 && i != fnum1 && i != fnum2) {
			bool isIntersect = ISTwoNurbsLineIntersect(tempLine, allLine[pol.objectInAllLines[i]]);
			assert(!isIntersect);
		}
	}

	line.CreatPartNurbsLine(tempLine, 1);
	return line;
}

PartNurbsLine GetPartNurbsLine(const Vec3 & p1, const Vec3 & frontp, double u, const varray<PartNurbsLine>& allLine, MyPolygon pol)
{
	PartNurbsLine line;
	Spline tempLine;
	Vec3 vec1, vec2;							//p1��p2�������ߵ�����
	int num1 = -1, num2 = -1;					//p1��frontp�ڶ�����е����

	for (int i = 0; i < pol.p_Points.size(); i++) {
		if (JudgeTwoPointsCoincide(p1, pol.p_Points[i].vetx)) {
			num1 = i;
			break;
		}
	}
	for (int i = 0; i < pol.p_Points.size(); i++) {
		if (JudgeTwoPointsCoincide(frontp, pol.p_Points[i].vetx)) {
			num2 = i;
			break;
		}
	}
	assert(num1 != -1 && num2 != -1);

	varray<Vec4> der;						//u�����꼰������
	allLine[pol.objectInAllLines[num2]].PtsDerivs(u, 1, der);
	assert(der.size());
	Vec3 p2 = der[0];						//�ضϵ�����

	Vec3 midPoint = Vec3((p1.x + p2.x)*0.5, (p1.y + p2.y)*0.5, 1);

	vec1 = (p2 - p1).Normalize();
	vec2 = (p1 - p2).Normalize();



	Vec3 frV1, bkV1, frV2, bkV2;				//p1��p2ǰ�����ߵ�������
	Vec3 fstrV1, bstrV1, fstrV2, bstrV2;		//p1��p2ǰ��ֱ�ߵķ�������
	GetVecsByPoint(pol, num1, allLine, frV1, bkV1);
	//GetVecsByPoint(pol, num1, allLine, frV2, bkV2);

	GetFBcoord(pol, num1, fstrV1, bstrV1);
	fstrV1 = (fstrV1 - p1).Normalize();
	bstrV1 = (bstrV1 - p1).Normalize();
	//GetFBcoord(pol, num2, fstrV2, bstrV2);
	//fstrV2 = (fstrV2 - p1).Normalize();
	//bstrV2 = (bstrV2 - p1).Normalize();
	if (JudgeTwoPointsCoincide(frontp, allLine[pol.objectInAllLines[num2]].m_CtrlPts[0])) {
		//�������׸����Ƶ����ǰ��˵�
		bkV2 = der[1].Normalize();
		frV2 = der[1].Normalize();
		frV2.x *= -1.0;
		frV2.y *= -1.0;

		fstrV2 = (Vec3)((Vec3)allLine[pol.objectInAllLines[num2]].m_CtrlPts[0] - p2).Normalize();
		bstrV2 = (Vec3)(*(allLine[pol.objectInAllLines[num2]].m_CtrlPts.end() - 1) - (Vec4)p2).Normalize();
	}
	else {
		frV2 = der[1].Normalize();
		bkV2 = der[1].Normalize();
		bkV2.x *= -1.0;
		bkV2.y *= -1.0;
		bstrV2 = (Vec3)(allLine[pol.objectInAllLines[num2]].m_CtrlPts[0] - (Vec4)p2).Normalize();
		fstrV2 = (Vec3)(*(allLine[pol.objectInAllLines[num2]].m_CtrlPts.end() - 1) - (Vec4)p2).Normalize();
	}


	double agle1, stragle1, agle2, stragle2;	//�����������ֱ�߼�����������֮��ĽǶ�
	double pagle1, pstragle1, pagle2, pstragle2;//p1��p2���Ժ�����(ֱ��)�������������ֱ�����������Ƕ�

	agle1 = GetAngle(bkV1, frV1);
	agle2 = GetAngle(bkV2, frV2);
	stragle1 = GetAngle(bstrV1, fstrV1);
	stragle2 = GetAngle(bstrV2, fstrV2);

	pagle1 = GetAngle(bkV1, vec1);
	pstragle1 = GetAngle(bstrV1, vec1);
	pagle2 = GetAngle(bkV2, vec2);
	pstragle2 = GetAngle(bstrV2, vec2);

	int flag = 0, flag2 = 0;
	//�ж�p1�Ƿ����Ҫ��
	if (pagle1 <= agle1) {
		//p1,p2�����������߼нǷ�Χ��
		if (abs(pagle1 - 0.0) < 1e-2) {
			//����������
			flag = 1;
		}
		else if (abs(pagle1 - agle1) < 1e-2) {
			//����ǰ����
			flag = 2;
		}
	}
	else {
		//�����������߼нǷ�Χ��
		double ang1 = (2 * PI - agle1)*0.5;
		double ang2 = pagle1 - agle1;
		if (ang2 >= ang1) {
			//����������
			flag = 1;
		}
		else {
			//����ǰ����
			flag = 2;
		}
	}

	if (flag == 0) {
		//p1�����Ҫ��,�ж�p2��
		if (pagle2 <= agle2) {
			//p2,p1�����������߼нǷ�Χ��
			if (abs(pagle2 - 0.0) < 1e-2) {
				//����������
				flag2 = 1;
			}
			else if (abs(pagle2 - agle2) < 1e-2) {
				//����ǰ����
				flag2 = 2;
			}
		}
		else {
			//�����������߼нǷ�Χ��
			double ang1 = (2 * PI - agle2)*0.5;
			double ang2 = pagle2 - agle2;
			if (ang2 >= ang1) {
				//����������
				flag2 = 1;
			}
			else {
				//����ǰ����
				flag2 = 2;
			}
		}

		if (flag2 == 1) {
			//p2->p1���������ߣ�p2���м俿����p1��Ҫ������߿���
			double angle = GetAngle(bkV2, fstrV2);//����ǰֱ�����Ƕ�
			angle = (angle < agle2) ? angle : agle2;//ѡȡ��С�ĽǶ�
			angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
			Vec3 v2 = bkV2.RotateZ(angle).Normalize();//����p2�����߷�������

			//����p1�����߷�������
			angle = (pagle1 < pstragle1) ? pagle1 : pstragle1;
			angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
			Vec3 v1 = vec1.RotateZ(-angle).Normalize();

			double length = max(abs(p1.x - p2.x), abs(p1.y - p2.y));
			Vec3 tmp1 = p1 + v1 * 3 * length;
			Vec3 tmp2 = p2 + v2 * 3 * length;
			midPoint = GetTwoLinesCrossPoint(p1, tmp1, p2, tmp2);
		}

		else if (flag2 == 2) {
			//p2->p1����ǰ���ߣ�p2���м俿����p1��ǰ���߿���
			double angle = GetAngle(bstrV2, frV2);//��ֱǰ�������Ƕ�
			angle = (angle < agle2) ? angle : agle2;
			angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
			Vec3 v2 = frV2.RotateZ(-angle).Normalize();//����p2�����߷�������

			//����p1�����߷�������
			double ag1 = agle1 - pagle1;//�������ߵ�ǰ�����������ĽǶ�
			double strag1 = stragle1 - pstragle1;//�������ߵ�ǰֱ�������ĽǶ�
			assert(ag1 >= 0);
			assert(strag1 >= 0);
			angle = (ag1 < strag1) ? ag1 : strag1;
			angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
			Vec3 v1 = vec1.RotateZ(angle).Normalize();//��ʱ����ת

			double length = max(abs(p1.x - p2.x), abs(p1.y - p2.y));
			Vec3 tmp1 = p1 + v1 * 3 * length;
			Vec3 tmp2 = p2 + v2 * 3 * length;
			midPoint = GetTwoLinesCrossPoint(p1, tmp1, p2, tmp2);
		}
	}

	else if (flag == 1) {
		//p1->p2���������ߣ�p1���м俿����p2��Ҫ������߿���
		double angle = GetAngle(bkV1, fstrV1);//����ǰֱ�����Ƕ�
		angle = (angle < agle1) ? angle : agle1;//ѡȡ��С�ĽǶ�
		angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
		Vec3 v1 = bkV1.RotateZ(angle).Normalize();//����p1�����߷�������

		//����p2�����߷�������
		angle = (pagle2 < pstragle2) ? pagle2 : pstragle2;
		angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
		Vec3 v2 = vec2.RotateZ(-angle).Normalize();

		double length = max(abs(p1.x - p2.x), abs(p1.y - p2.y));
		Vec3 tmp1 = p1 + v1 * 3 * length;
		Vec3 tmp2 = p2 + v2 * 3 * length;
		midPoint = GetTwoLinesCrossPoint(p1, tmp1, p2, tmp2);
	}

	else if (flag == 2) {
		//p1->p2����ǰ���ߣ�p1���м俿����p2��Ҫ��ǰ���߿���
		double angle = GetAngle(bstrV1, frV1);//��ֱǰ�������Ƕ�
		angle = (angle < agle1) ? angle : agle1;
		angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
		Vec3 v1 = frV1.RotateZ(-angle).Normalize();//����p2�����߷�������

		//����p2�����߷�������
		double ag2 = agle2 - pagle2;//�������ߵ�ǰ�����������ĽǶ�
		double strag2 = stragle2 - pstragle2;//�������ߵ�ǰֱ�������ĽǶ�
		assert(ag2 >= 0);
		assert(strag2 >= 0);
		angle = (ag2 < strag2) ? ag2 : strag2;
		angle = (angle > (45 / 180)*PI) ? (0.3*angle) : (0.5*angle);
		Vec3 v2 = vec2.RotateZ(angle).Normalize();//��ʱ����ת

		double length = max(abs(p1.x - p2.x), abs(p1.y - p2.y));
		Vec3 tmp1 = p1 + v1 * 3 * length;
		Vec3 tmp2 = p2 + v2 * 3 * length;
		midPoint = GetTwoLinesCrossPoint(p1, tmp1, p2, tmp2);
	}

	if (flag == 0 && flag2 == 0) {
		//ֱ�Ӳ���ֱ����ʽ��NURBS��������
		tempLine.m_Degree = 1;
		tempLine.m_Knots.push_back(0);
		tempLine.m_Knots.push_back(0);
		tempLine.m_Knots.push_back(1);
		tempLine.m_Knots.push_back(1);
		tempLine.m_CtrlPts.push_back(Vec4(p1));
		tempLine.m_CtrlPts.push_back(Vec4(p2));
		tempLine.DegreeElevate(2);
	}
	else {
		//����������ʽ��NURBS��������
		tempLine.m_Degree = 2;
		tempLine.m_Knots.push_back(0);
		tempLine.m_Knots.push_back(0);
		tempLine.m_Knots.push_back(0);
		tempLine.m_Knots.push_back(1);
		tempLine.m_Knots.push_back(1);
		tempLine.m_Knots.push_back(1);
		tempLine.m_CtrlPts.push_back(Vec4(p1));
		tempLine.m_CtrlPts.push_back(Vec4(midPoint.x,midPoint.y,midPoint.z, 1));
		tempLine.m_CtrlPts.push_back(Vec4(p2));
	}

	//�ж��½������Ƿ�Ͷ�������������н���
	int fnum1 = pol.GetFrontNumber(num1);
	for (int i = 0; i < pol.objectInAllLines.size(); ++i) {
		if (i != num1 && i != num2 && i != fnum1) {
			bool isIntersect = ISTwoNurbsLineIntersect(tempLine, allLine[pol.objectInAllLines[i]]);
			assert(!isIntersect);
		}
	}

	line.CreatPartNurbsLine(tempLine, 1);
	return line;
}

double GetPartU(const Spline & line, const Vec3 & begin, const Vec3 & end, const Vec3 & crossp)
{
	Vec3 vec1, vec2;
	vec1 = crossp - begin;
	vec2 = end - begin;
	double d1, d2;
	d1 = vec1.Magnitude();
	d2 = vec2.Magnitude();
	double u = d1 / d2;

	Vec3 tempPoint = line.m_CtrlPts[0];
	if (!JudgeTwoPointsCoincide(begin, tempPoint))
	{
		u = 1 - u;
	}
	return u;
}

void CreatArc(const Vec3 & p1, const Vec3 & p2, const Vec3 & centre, Spline & line)
{
	line.m_CtrlPts.clear();
	line.m_Degree = 2;
	line.m_Knots.clear();
	Vec3 mid = (p1 + p2)*0.5;
	Vec3 vec1, vec2;
	vec1 = p1 - centre;
	vec2 = mid - centre;
	double w = vec1.Dot(vec2) / (vec1.Magnitude()*vec2.Magnitude());
	double lth1, lth2;
	lth1 = vec1.Magnitude();
	lth2 = lth1 / w;
	vec2 = vec2.Normalize();
	vec2 *= lth2;
	
	Vec3 midPoint = vec2 + centre;
	line.m_Knots.push_back(0);
	line.m_Knots.push_back(0);
	line.m_Knots.push_back(0);
	line.m_Knots.push_back(1);
	line.m_Knots.push_back(1);
	line.m_Knots.push_back(1);

	line.m_CtrlPts.push_back(Vec4(p1.x,p1.y,p1.z, 1));
	line.m_CtrlPts.push_back(Vec4(midPoint.x, midPoint.y, midPoint.z, w));
	line.m_CtrlPts.push_back(Vec4(p2.x, p2.y, p2.z, 1));
}

void QuadOrder(const MyPolygon &quad,  varray<PartNurbsLine> &allLine, varray<Spline> &coonsLine)
{
	coonsLine.clear();
	Spline tempLine;
	Vec3 tempPoint;
	if(quad.p_Points.size()!=4)
	{
		cout << "\n�ı�������ʱ��������Ϊ���ı���" << endl;
		return;
	}
	//��һ��
	allLine[quad.objectInAllLines[0]].CreatNurbsLine(tempLine);
	tempPoint = tempLine.m_CtrlPts[0];
	if (!JudgeTwoPointsCoincide(quad.p_Points[0].vetx, tempPoint))
	{
		tempLine.CruveReverse();
	}
	coonsLine.push_back(tempLine);
	//�ڶ���
	allLine[quad.objectInAllLines[3]].CreatNurbsLine(tempLine);
	tempPoint = tempLine.m_CtrlPts[0];
	if (!JudgeTwoPointsCoincide(quad.p_Points[0].vetx, tempPoint))
	{
		tempLine.CruveReverse();
	}
	coonsLine.push_back(tempLine);
	//������
	allLine[quad.objectInAllLines[2]].CreatNurbsLine(tempLine);
	tempPoint = tempLine.m_CtrlPts[0];
	if (!JudgeTwoPointsCoincide(quad.p_Points[3].vetx, tempPoint))
	{
		tempLine.CruveReverse();
	}
	coonsLine.push_back(tempLine);
	//������
	allLine[quad.objectInAllLines[1]].CreatNurbsLine(tempLine);
	tempPoint = tempLine.m_CtrlPts[0];
	if (!JudgeTwoPointsCoincide(quad.p_Points[1].vetx, tempPoint))
	{
		tempLine.CruveReverse();
	}
	coonsLine.push_back(tempLine);
}

void QuadOrder(const varray<int> &quadNums, const varray<Spline>& allLine, varray<Spline>& coonsLine)
{
	coonsLine.clear();
	Spline tempLine;
	Vec3 tempPoint;
	assert(quadNums.size() == 4);
	//if (quadNums.size() != 4)
	//{
	//	cout << "\n�ı�������ʱ��������Ϊ���ı���" << endl;
	//	return;
	//}

	bool nextLine = false;
	map<int, bool> founded;
	for (int i = 0; i < 4; ++i) {
		founded[i] = false;//��ʼ��
	}
	//��һ��
	tempLine = allLine[quadNums[0]];
	tempPoint = tempLine.m_CtrlPts[0];
	coonsLine.push_back(tempLine);
	founded[0] = true;
	//�ڶ���
	for (int i = 1; i < 4; ++i) {
		nextLine = false;
		if (founded[i]) continue;
		tempLine = allLine[quadNums[i]];
		//������� �����׺�β ����ͬ��
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
	coonsLine.push_back(tempLine);
	tempPoint = *(tempLine.m_CtrlPts.end() - 1);
	//������
	for (int i = 1; i < 4; ++i) {
		nextLine = false;
		if (founded[i]) continue;
		tempLine = allLine[quadNums[i]];
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
	coonsLine.push_back(tempLine);
	tempPoint = *(tempLine.m_CtrlPts.end() - 1);
	//������
	for (int i = 1; i < 4; ++i) {
		nextLine = false;
		if (founded[i]) continue;
		tempLine = allLine[quadNums[i]];
		if (JudgeTwoPointsCoincide(tempLine.m_CtrlPts[0], tempPoint))
		{
			tempLine.CruveReverse();
			nextLine = true;
			founded[i] = true;
		}
		if (nextLine) break;
		if (JudgeTwoPointsCoincide(*(tempLine.m_CtrlPts.end() - 1), tempPoint)) {//*�Ե�ַʹ�ã�ȡֵ����� 
			nextLine = true;
			founded[i] = true;
		}
		if (nextLine) break;
	}
	coonsLine.push_back(tempLine);
}

bool JudgeTwoPointsCoincide(const Vec3 & p1, const Vec3 & p2)
{
	if ((fabs(p1.x - p2.x) < 1e-4) && (fabs(p1.y - p2.y) < 1e-4) && (fabs(p1.z - p2.z) < 1e-4))
		return true;
	else
		return false;
}

bool JudgeTwoLinesCoincide(const Spline & L1, const Spline & L2)
{
	Vec3 p1, p2, q1, q2;
	p1 = L1.m_CtrlPts[0];
	p2 = *(L1.m_CtrlPts.end() - 1);
	q1 = L2.m_CtrlPts[0];
	q2 = *(L2.m_CtrlPts.end() - 1);

	if (JudgeTwoPointsCoincide(p1, q1))
	{
		if(JudgeTwoPointsCoincide(p2,q2))
		{
			return true;
		}
	}
	else
	{
		if (JudgeTwoPointsCoincide(p1, q2))
		{
			if(JudgeTwoPointsCoincide(p2, q1))
			{
				return true;
			}
		}
	}
	return false;	
}


//��������
double GetAngle(Vec3 vec1, Vec3 vec2)
{
	//������λ�����н�
	assert(abs(vec1.z - 0.0) < 1e-10);
	assert(abs(vec2.z - 0.0) < 1e-10);
	double angle;
	vec1 = vec1.Normalize();
	vec2 = vec2.Normalize();
	double dot = vec1.Dot(vec2);//������
	if (fabs(fabs(dot) - 1) < 1e-15)
	{
		return dot > 0 ? 0 : PI;
	}
	angle = (double)acos(dot);

	Vec3 cross = vec1.Cross(vec2);
	return cross.z > 0 ? angle : 2 * PI - angle;
}

//�õ����߽����
Vec3 GetTwoLinesCrossPoint(const Vec3 & p1, const Vec3 & p2, const Vec3 & p3, const Vec3 & p4)
{
	Vec3 crossp;
	double a1, b1, c1, a2, b2, c2, D;
	a1 = p1.y - p2.y;
	b1 = p2.x - p1.x;
	c1 = p1.x*p2.y - p1.y*p2.x;

	a2 = p3.y - p4.y;
	b2 = p4.x - p3.x;
	c2 = p3.x*p4.y - p3.y*p4.x;

	D = a1 * b2 - a2 * b1;
	double crossx = (b1*c2 - b2 * c1) / D;
	double crossy = (c1*a2 - c2 * a1) / D;
	crossp = Vec3(crossx, crossy, p1.z);
	return crossp;
}

//��ö����ָ����ŵ��ǰ���
void QuadPart::GetVecs(const MyPolygon &inputPolygen, int pointNumber, Vec3 & frontVec, Vec3 & behindVec)
{
	int pLength = inputPolygen.p_Points.size();
	if (pointNumber == 0)
	{
		frontVec = inputPolygen.p_Points[pLength - 1].vetx - inputPolygen.p_Points[pointNumber].vetx;
		behindVec = inputPolygen.p_Points[pointNumber + 1].vetx - inputPolygen.p_Points[pointNumber].vetx;
	}
	else if (pointNumber == pLength - 1)
	{
		frontVec = inputPolygen.p_Points[pointNumber - 1].vetx - inputPolygen.p_Points[pointNumber].vetx;
		behindVec = inputPolygen.p_Points[0].vetx - inputPolygen.p_Points[pointNumber].vetx;
	}
	else
	{
		frontVec = inputPolygen.p_Points[pointNumber - 1].vetx - inputPolygen.p_Points[pointNumber].vetx;
		behindVec = inputPolygen.p_Points[pointNumber + 1].vetx - inputPolygen.p_Points[pointNumber].vetx;
	}
}

//�����ָ����ű�
double QuadPart::GetFuncSita(const MyPolygon &inputPolygen, int pointNumber, Vec3 vecAB)
{
	//����Ȩ�ؼ���
	if (inputPolygen.convexHull[pointNumber] == 1) //����
	{
		Vec3 fVec, bVec; //ǰ������
		GetVecs(inputPolygen, pointNumber, fVec, bVec);
		double sita = GetAngle(bVec, vecAB);
		double sita1 = GetAngle(bVec, -fVec);
		double sita3 = PI;
		double sita4 = GetAngle(bVec, fVec);
		double sita2 = sita4 / 2.0;
		if (sita > 0 && sita <= sita1)
		{
			return sita4 * sita / sita1;
		}
		else if (sita >= sita1 && sita <= sita2)
		{
			double up = sita - sita1;
			double down = sita2 - sita1;
			return sita4 * (1.0 + (up / down));
		}
		else if (sita >= sita2 && sita <= sita3)
		{
			double up = sita - sita2;
			double down = sita3 - sita2;
			return sita4 * (2.0 - (up / down));
		}
		else if (sita >= sita3 && sita <= sita4)
		{
			double up = sita - sita3;
			double down = sita4 - sita3;
			return sita4 * (1.0 - (up / down));
		}
		else
		{
			return 0;
		}
	}
	//͹��Ȩ�ؼ���
	else if (inputPolygen.convexHull[pointNumber] == 0)
	{
		Vec3 fVec, bVec; //ǰ������
		GetVecs(inputPolygen, pointNumber, fVec, bVec);
		double sita = GetAngle(bVec, vecAB);
		double sita2 = GetAngle(bVec, fVec);
		double sita1 = sita2 / 2.0;
		double index = 3.0*sita2 / PI - 1.0;
		double coefficient = pow(2.0, index);//����2��index����	

		if (sita > 0 && sita < sita1)
		{
			return coefficient * sita2*sita / sita1;
		}
		else if (sita >= sita1 && sita <= sita2)
		{
			double up = sita - sita1;
			double down = sita2 - sita1;
			return coefficient * sita2*(1 - up / down);
		}
		else
		{
			return 0;
		}
	}
	return 0;
}


//��ʼ�㡢�����㡢�е㡢���㣨�����Ϸ���
double QuadPart::GetFuncSita(const Vec3 & begin, const Vec3 & end, const Vec3 & mid, const Vec3 & pt)
{
	Vec3 vec1, vec2, vec3;
	double ang;
	vec1 = begin - mid;
	vec2 = end - mid;
	vec3 = pt - mid;
	vec1 = vec1.Normalize();
	vec2 = vec2.Normalize();
	vec3 = vec3.Normalize();
	ang = GetAngle(vec2, vec1);

	if (ang>=PI) //����
	{
		Vec3 fVec, bVec; //ǰ������
		fVec = vec1;
		bVec = vec2;
		double sita = GetAngle(bVec, vec3);
		double sita1 = GetAngle(bVec, -fVec);
		double sita3 = PI;
		double sita4 = GetAngle(bVec, fVec);
		double sita2 = sita4 / 2.0;
		if (sita > 0 && sita <= sita1)
		{
			return sita4 * sita / sita1;
		}
		else if (sita >= sita1 && sita <= sita2)
		{
			double up = sita - sita1;
			double down = sita2 - sita1;
			return sita4 * (1.0 + (up / down));
		}
		else if (sita >= sita2 && sita <= sita3)
		{
			double up = sita - sita2;
			double down = sita3 - sita2;
			return sita4 * (2.0 - (up / down));
		}
		else if (sita >= sita3 && sita <= sita4)
		{
			double up = sita - sita3;
			double down = sita4 - sita3;
			return sita4 * (1.0 - (up / down));
		}
		else
		{
			return 0;
		}
	}

	else if (ang < PI)
	{
		Vec3 fVec, bVec; //ǰ������
		fVec = vec1;
		bVec = vec2;
		double sita = GetAngle(bVec, vec3);
		double sita2 = GetAngle(bVec, fVec);
		double sita1 = sita2 / 2.0;
		double index = 3.0*sita2 / PI - 1.0;
		double coefficient = pow(2.0, index);

		if (sita > 0 && sita < sita1)
		{
			return coefficient * sita2*sita / sita1;
		}
		else if (sita >= sita1 && sita <= sita2)
		{
			double up = sita - sita1;
			double down = sita2 - sita1;
			return coefficient * sita2*(1 - up / down);
		}
		else
		{
			return 0;
		}
	}
	return 0;
}

void QuadPart::OrderVisiblePoints(const MyPolygon &inputPolygen, varray<VisiblePoints>& visiblePointPair)
{
	varray<VisiblePoints> concavePair, convexPair, tempPair; //��������͹����
	DandI temp; 
	varray<DandI> funcSita, completeFuncSita; //��������
	int vppLength = visiblePointPair.size();
	if (vppLength != 0)
	{
		//��ð���͹����
		for (int i = 0; i < vppLength; i++)
		{
			if (visiblePointPair[i].ifConvex == 1)
			{
				concavePair.push_back(visiblePointPair[i]);
			}
			else if (visiblePointPair[i].ifConvex == 0)
			{
				convexPair.push_back(visiblePointPair[i]);
			}
		}
		int caveLength = concavePair.size();
		int vexLength = convexPair.size();
		
		if (caveLength != 0) //�����ϲ�Ϊ��
		{
			tempPair.clear();
			funcSita.clear();
			for (int i = 0; i < caveLength; i++)
			{
				int firstP = concavePair[i].firstPoint;
				int secondP = concavePair[i].secondPoint;
				Vec3 AB = inputPolygen.p_Points[secondP].vetx - inputPolygen.p_Points[firstP].vetx;
				Vec3 BA = inputPolygen.p_Points[firstP].vetx - inputPolygen.p_Points[secondP].vetx;
				double funcA = GetFuncSita(inputPolygen, firstP, AB);
				double funcB = GetFuncSita(inputPolygen, secondP, BA);
				double funcAB = funcA + funcB;
				temp.funcVar = funcAB;
				temp.order = i;
				funcSita.push_back(temp);
				concavePair[i].w = funcAB;
			}
			//�԰����Ͻ��н�������
			int fsLength = funcSita.size();
			if (fsLength > 1)
			{
				for (int i = 0; i < fsLength -1; i++)
				{
					for (int j = 0; j < fsLength - 1 - i;j++)
					{
						if (funcSita[j].funcVar < funcSita[j + 1].funcVar)
						{
							temp.funcVar = funcSita[j].funcVar;
							temp.order = funcSita[j].order;
							funcSita[j].funcVar = funcSita[j + 1].funcVar;
							funcSita[j].order = funcSita[j + 1].order;
							funcSita[j + 1].funcVar = temp.funcVar;
							funcSita[j + 1].order = temp.order;
						}
					}
				}
			}	
			for (int i = 0; i < fsLength; i++)
			{
				tempPair.push_back(concavePair[funcSita[i].order]);
			}
			concavePair.clear();

			concavePair = tempPair;
			//for (int i = 0; i < tempPair.size(); i++)
			//{
			//	concavePair.push_back(tempPair[i]);
			//}
			tempPair.clear();
		}
	
		//����͹����
		if (vexLength > 1)
		{
			tempPair.clear();
			funcSita.clear();
			varray<double> allS; //����S��ֵ����
			varray<double> allD; //���п��ӵ�֮��ľ��뼯��
			double DMax; //allD�����ֵ
			double cs = 1.0, ct = 0.55; //���Ʋ���
			for (int i = 0; i < vexLength; i++)
			{

				varray<double> allAngle;
				double sita;
				Vec3 fVec, bVec;
				int firstP = convexPair[i].firstPoint;
				int secondP = convexPair[i].secondPoint;
				Vec3 AB = inputPolygen.p_Points[secondP].vetx - inputPolygen.p_Points[firstP].vetx;
				Vec3 BA = inputPolygen.p_Points[firstP].vetx - inputPolygen.p_Points[secondP].vetx;
				double DistanceAB = AB.Magnitude(); //AB�ĳ���
				//A�����нǶȼ���
				GetVecs(inputPolygen, firstP, fVec, bVec);
				sita = GetAngle(AB, fVec);
				allAngle.push_back(sita);
				sita = GetAngle(bVec, AB);
				allAngle.push_back(sita);
				//B�����нǶȼ���
				GetVecs(inputPolygen, secondP, fVec, bVec);
				sita = GetAngle(bVec, BA);
				allAngle.push_back(sita);
				sita = GetAngle(BA, fVec);
				allAngle.push_back(sita);
				//ȡ4���Ƕ���Сֵ
				for (int i = 0; i < 4; i++)
				{
					sita = min(sita, allAngle[i]);
				}
				double S = sita * 2 / PI;
				allS.push_back(S);
				allD.push_back(DistanceAB);
			}
			DMax = allD[0];
			for (int i = 0; i < allD.size(); i++)
			{
				DMax = max(DMax, allD[i]);
			}
			//���Ȩ����ֵ
			for (int i = 0; i < vexLength; i++)
			{
				double S = allS[i];
				double T = (DMax - allD[i]) / DMax;
				temp.funcVar = cs * S + ct * T;
				temp.order = i;
				funcSita.push_back(temp);
				convexPair[i].w = temp.funcVar;
			}
			//��͹���Ͻ�������
			int fsLength = funcSita.size();
			for (int i = 0; i < fsLength - 1; i++)
			{
				for (int j = 0; j < fsLength - 1 - i; j++)
				{
					if (funcSita[j].funcVar < funcSita[j + 1].funcVar)
					{
						temp.funcVar = funcSita[j].funcVar;
						temp.order = funcSita[j].order;
						funcSita[j].funcVar = funcSita[j + 1].funcVar;
						funcSita[j].order = funcSita[j + 1].order;
						funcSita[j + 1].funcVar = temp.funcVar;
						funcSita[j + 1].order = temp.order;
					}
				}
			}
			for (int i = 0; i < fsLength; i++)
			{
				tempPair.push_back(convexPair[funcSita[i].order]);
			}
			convexPair.clear();
			convexPair = tempPair;
			//for (int i = 0; i < tempPair.size(); i++)
			//{
			//	convexPair.push_back(tempPair[i]);
			//}
			tempPair.clear();
		}
		//������������Ŀ��ӵ�Լ��ϣ������ɴ�С(Ȩ����)���Ȱ���͹˳��
		visiblePointPair.clear();
		if (caveLength != 0)
		{
			for (int i = 0; i < caveLength; i++)
			{
				visiblePointPair.push_back(concavePair[i]);
			}
		}
		if (vexLength != 0)
		{
			for (int i = 0; i < vexLength; i++)
			{
				visiblePointPair.push_back(convexPair[i]);
			}
		}
	}
	else
	{
		return;
	}
}

void QuadPart::SortVisiblePoints(MyPolygon &inputPolygen, varray<VisiblePoints>& visiblePointPair, varray<VisiblePoints>& sharpVisiblePointPair)
{
	int vppLength = visiblePointPair.size();
	//VisiblePoints tempPair;
	varray<VisiblePoints> tempPairs1, tempPairs2;
	if (vppLength != 0)
	{
		//��������
		for (int i = 0; i < vppLength; i++)
		{
			if (visiblePointPair[i].ifSharp == 0)
			{
				tempPairs1.push_back(visiblePointPair[i]);
			}
			else
			{
				tempPairs2.push_back(visiblePointPair[i]);
			}
		}

		//�Ǽ��������
		OrderVisiblePoints(inputPolygen, tempPairs1);
		//���������
		OrderVisiblePoints(inputPolygen, tempPairs2);

		//���÷Ǽ����㼯�ϼ���
		visiblePointPair.clear();
		sharpVisiblePointPair.clear();
		if (tempPairs1.size() != 0)
		{
			for (int i = 0; i < tempPairs1.size(); i++)
			{
				visiblePointPair.push_back(tempPairs1[i]);
			}
		}
		if (tempPairs2.size() != 0)
		{
			for (int i = 0; i < tempPairs2.size(); i++)
			{
				sharpVisiblePointPair.push_back(tempPairs2[i]);
			}
		}
	}
	else
	{
		return;
	}
}

int QuadPart::GetNextPoint(const MyPolygon & inputPolygen, const int & pNumber)
{
	int length = inputPolygen.p_Points.size();
	if(pNumber==length-1)
	{
		return 0;
	}
	else
	{
		return pNumber + 1;
	}
}

double QuadPart::GetPartPointU(const MyPolygon & inputPolygen, const varray<PartNurbsLine>& allLines, const Vec3 & begin, const Vec3 & end, const Vec3 & crossp)
{
	Vec3 vec1, vec2;
	vec1 = crossp - begin;
	vec2 = end - begin;
	double d1, d2;
	d1 = vec1.Magnitude();//��������
	d2 = vec2.Magnitude();
	double u = d1 / d2;
	PartNurbsLine tempLine;

	int flag=-1;
	for (int i = 0; i < inputPolygen.p_Points.size(); i++)
	{
		if (JudgeTwoPointsCoincide(begin, inputPolygen.p_Points[i].vetx))
			flag = i;//�ҳ���ʼ�����������±�
	}

	if (flag >= 0)
	{
		tempLine = allLines[inputPolygen.objectInAllLines[flag]];//�������
	}

	Vec3 tempPoint = tempLine.m_CtrlPts[0];
	if (!JudgeTwoPointsCoincide(begin, tempPoint))//�����غϣ���һ��ȡu
	{
		u = 1 - u;
	}
	return u;
}

void QuadPart::CalPartLines(const MyPolygon & inputPolygen,  const varray<PartNurbsLine> &allLines, varray<VisibleLine>& partLine, const int &mode)
{
	partLine.clear();
	varray<int> part;//�ɷָ��ߵ����
	varray<VisiblePoints> vp;//���п��ӵ�
	VisibleLine vbL;

	FindAnotherVp(inputPolygen, vp);
	int vpsLth = vp.size();
	for (int i = 0; i < inputPolygen.p_Points.size(); i++)
	{
		if (allLines[inputPolygen.objectInAllLines[i]].ifSegLine == 1)
			part.push_back(i);
	}
	//Ѱ�ҿɷָ��ߵĿ��ӵ�
	if (part.size() == 0 )
	{
		return;
	}
	else
	{
		for (int i = 0; i < part.size(); i++)
		{
			int p1 = part[i];
			int p2 = GetNextPoint(inputPolygen, p1);//�õ���һ����
			varray<int> temp1, temp2;
			for (int j = 0; j < vp.size(); j++)
			{
				if (vp[j].firstPoint == p1)
				{
					temp1.push_back(vp[j].secondPoint);
				}
				if (vp[j].secondPoint == p1)
				{
					temp1.push_back(vp[j].firstPoint);
				}
				if (vp[j].firstPoint == p2)
				{
					temp2.push_back(vp[j].secondPoint);
				}
				if (vp[j].secondPoint == p2)
				{
					temp2.push_back(vp[j].firstPoint);
				}
			}
			if (temp1.size() == 0)
			{
				continue;
			}
			if (temp2.size() == 0)
			{
				continue;
			}

			for (int m = 0; m < temp1.size(); m++)
			{
				for (int n = 0; n < temp2.size(); n++)
				{
					if (temp1[m] == temp2[n])
					{
						vbL.pt = temp1[m];
						vbL.beginP = p1;
						vbL.endP = p2;
						partLine.push_back(vbL);
					}
				}
			}
		}
	}

	//����ָ��ߵ�u ������Ȩֵ
	if (partLine.size() == 0)
	{
		return;
	}
	else
	{
		for (int i = 0; i < partLine.size(); i++)
		{
			Vec3 p1, p2, pT, p4, tempP;
			p1 = inputPolygen.p_Points[partLine[i].beginP].vetx;
			p2 = inputPolygen.p_Points[partLine[i].endP].vetx;
			pT = inputPolygen.p_Points[partLine[i].pt].vetx;

			//�����ж��Ƿ���ڴ���,(����pT�ı������߽���֣���������ݽ������߸ĳ�����)
			p4 = GetFootPerpendicular(p1, p2, pT);//�����
			if (JudegeFootPerpendicular(p1, p2, p4))
			{
				if (mode == 1)
				{
					Vec3 frontVec, lastVec;
					GetVecs(inputPolygen, partLine[i].pt, frontVec, lastVec);

					bool sflag=LinkLineSharpJudge(frontVec, lastVec, pT, p4);
					if (sflag == false) //���ڼ��
					{
						partLine[i].u = -1;
						continue;
					}
				}
				partLine[i].u = GetPartPointU(inputPolygen, allLines, p1, p2, p4);
				tempP = p4 - pT;
				tempP = tempP.Normalize();
				double w1 = GetFuncSita(inputPolygen, partLine[i].pt, tempP);
				double w2 = GetFuncSita(p1, p2, p4, pT);
				partLine[i].w = w1 + w2;
				continue;
			}
			else
			{
				//���н�ƽ�����ж�
				Vec3 froVec, befVec, midVec;
				GetVecs(inputPolygen, partLine[i].pt, froVec, befVec);
				midVec=GetAngularBisectorVec(froVec, befVec);
				double angle = GetAngle(befVec, froVec);
				midVec = midVec * 10;
				tempP = pT + midVec;
				p4 = GetTwoLinesCrossPoint(p1, p2, pT, tempP);
				if (JudegeFootPerpendicular(p1, p2, p4))
				{
					if (mode == 1)
					{
						Vec3 frontVec, lastVec;
						GetVecs(inputPolygen, partLine[i].pt, frontVec, lastVec);

						bool sflag = LinkLineSharpJudge(frontVec, lastVec, pT, p4);
						if (sflag == false) //���ڼ��
						{
							partLine[i].u = -1;
							continue;
						}
					}
					partLine[i].u = GetPartPointU(inputPolygen, allLines, p1, p2, p4);
					tempP = p4 - pT;
					tempP = tempP.Normalize();
					double w1 = GetFuncSita(inputPolygen, partLine[i].pt, tempP);
					double w2 = GetFuncSita(p1, p2, p4, pT);
					partLine[i].w = w1 + w2;
					
					cout << "��ע���Ȩ��:" << partLine[i].w << endl;//������ע���Ȩ��

					continue;
				}
				else
				{
					//ѡ���е�
					partLine[i].u = 0.5;
					p4 = p1 + p2;
					p4 = 0.5*p4;//ȡ�е�
					if (mode == 1)
					{
						Vec3 frontVec, lastVec;
						GetVecs(inputPolygen, partLine[i].pt, frontVec, lastVec);

						bool sflag = LinkLineSharpJudge(frontVec, lastVec, pT, p4);
						if (sflag == false) //���ڼ��
						{
							partLine[i].u = -1;
							continue;
						}
					}

					tempP = p4 - pT;
					tempP = tempP.Normalize();
					double w1 = GetFuncSita(inputPolygen, partLine[i].pt, tempP);
					double w2 = GetFuncSita(p1, p2, p4, pT);
					partLine[i].w = w1 + w2;

					cout << "��ע���Ȩ��:"<<partLine[i].w << endl;//��ע���Ȩ��

					continue;
				}
			}
		}

		
		int partLth = partLine.size();
		if (partLth > 0)
		{
			varray<VisibleLine> tempVbls;
			for (int i = 0; i < partLth; i++)
			{
				if (partLine[i].u != -1)
				{
					tempVbls.push_back(partLine[i]);
				}
			}
			partLine.clear();
			partLine = tempVbls;
			tempVbls.clear();
		}

		//��partLine�������򣬰���w�Ӵ�С
		partLth = partLine.size();
		if (partLth > 1)
		{
			for (int i = 0; i < partLth - 1; i++)
			{
				for (int j = 0; j < partLth - i - 1; j++)
				{
					cout << "��partline��Ȩֵ����..." << endl;
					if (partLine[j].w < partLine[j + 1].w)
					{
						vbL = partLine[j];
						partLine[j] = partLine[j + 1];
						partLine[j + 1] = vbL;
					}
				}
			}
		}
	}

}

void QuadPart::CreatPartNurbsLine(const Vec3 & p1, const Vec3 & p2, PartNurbsLine & line)
{
	Spline tempLine;
	tempLine.m_Degree = 1;
	tempLine.m_Knots.push_back(0);
	tempLine.m_Knots.push_back(0);
	tempLine.m_Knots.push_back(1);
	tempLine.m_Knots.push_back(1);
	tempLine.m_CtrlPts.push_back(Vec4(p1));
	tempLine.m_CtrlPts.push_back(Vec4(p2));
	tempLine.DegreeElevate(2);//����
	line.CreatPartNurbsLine(tempLine, 1);//�������ߣ�Ĭ������Ϊ�ɷָ�	

	cout << "�������ߣ�Ĭ������Ϊ�ɷָ�" << endl;
}

bool QuadPart::ConncetTwoPoints(MyPolygon &pendingPol, varray<MyPolygon>& quadPol, varray<MyPolygon>& morePol, varray<PartNurbsLine>& allLines, VisiblePoints & visibleP)
{
	int qLth, mLth, alLth, pengingLth;
	qLth = quadPol.size(); //�ı��θ���
	mLth = morePol.size(); //����θ���
	alLth = allLines.size(); //�����ߵ���Ŀ
	pengingLth = pendingPol.p_Points.size();//����������

	PartNurbsLine tempLine;
	MyPolygon tempPol;
	varray<MyPolygon> outputPol;


	//����������NURBS����
	tempLine = GetPartNurbsLine(pendingPol.p_Points[visibleP.firstPoint].vetx, pendingPol.p_Points[visibleP.secondPoint].vetx, allLines, pendingPol);

	if (tempLine.m_CtrlPts.size() == 0) {
		//��ʾ���ص��ǿ�����
		return false;
	}
	allLines.push_back(tempLine);
	alLth = allLines.size();

	//��һ������0�ŵ�
	if (visibleP.firstPoint == 0)
	{
		//�ʷֳ��ĵ�һ�������
		tempPol.p_Points.push_back(pendingPol.p_Points[0]);
		tempPol.objectInAllLines.push_back(alLth-1);
		for (int i = visibleP.secondPoint; i < pengingLth; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		outputPol.push_back(tempPol);
		//�ڶ��������
		tempPol.convexHull.clear();
		tempPol.objectInAllLines.clear();
		tempPol.p_NurbsLines.clear();
		tempPol.p_Points.clear();
		for (int i = 0; i < visibleP.secondPoint; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		tempPol.p_Points.push_back(pendingPol.p_Points[visibleP.secondPoint]);
		tempPol.objectInAllLines.push_back(alLth - 1);
		outputPol.push_back(tempPol);
	}

	else
	{
		for (int i = 0; i < visibleP.firstPoint; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		tempPol.p_Points.push_back(pendingPol.p_Points[visibleP.firstPoint]);
		tempPol.objectInAllLines.push_back(alLth - 1);
		for (int i = visibleP.secondPoint; i < pengingLth; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		outputPol.push_back(tempPol);

		//�ڶ��������
		tempPol.convexHull.clear();
		tempPol.objectInAllLines.clear();
		tempPol.p_NurbsLines.clear();
		tempPol.p_Points.clear();
		for (int i = visibleP.firstPoint; i < visibleP.secondPoint; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		tempPol.p_Points.push_back(pendingPol.p_Points[visibleP.secondPoint]);
		tempPol.objectInAllLines.push_back(alLth - 1);
		outputPol.push_back(tempPol);
	}

	//ɾ�����һ�������
	morePol.pop_back();
	//�ж��ʷֳ��Ķ�����Ƿ�����ı���
	for (int i = 0; i < 2; i++)
	{
		int pLth = outputPol[i].p_Points.size();
		if (pLth == 4)
		{
			//���ж��Ƿ���ڰ���
			JudgeConcavePoint(outputPol[i]);
			for (int j = 0; j < 4; j++)
			{
				if (outputPol[i].convexHull[j] == 1)
				{
					return false;
				}
			}
			if (largeAngle) {
				//�ж��Ƿ��нǶȹ���
				Vec3 vec1, vec2;
				for (int j = 0; j < 4; j++) {
					GetVecs(outputPol[i], j, vec1, vec2);
					double angle = GetAngle(vec2, vec1);
					if (angle > lag)
						return false;
				}
			}
			for (int j = 0; j < outputPol[i].objectInAllLines.size(); j++) {
				//20210205�޸�
				allLines[outputPol[i].objectInAllLines[j]].ifSegLine = 0;
			}
			//allLines[alLth - 1].ifSegLine = 0;
			quadPol.push_back(outputPol[i]);
		}
		else
		{
			morePol.push_back(outputPol[i]);
		}
	}
	return true;
}

bool QuadPart::ConnectLineWithPoint(MyPolygon & pendingPol, varray<MyPolygon>& quadPol, varray<MyPolygon>& morePol, varray<PartNurbsLine>& allLines, VisibleLine & visibleL)
{
	int qLth, mLth, alLth, pengingLth;
	qLth = quadPol.size(); //�ı��θ���
	mLth = morePol.size(); //����θ���
	alLth = allLines.size(); //�����ߵ���Ŀ
	pengingLth = pendingPol.p_Points.size();

	PartNurbsLine tempLine;
	MyPolygon tempPol;
	varray<MyPolygon> outputPol;
	varray<PartNurbsLine> partLines;
	SubPoint tempSubP;

	if (visibleL.beginP == 0 && visibleL.endP == (pengingLth-1))
	{
		int pn = visibleL.beginP;
		visibleL.beginP = visibleL.endP;
		visibleL.endP = pn;
	}
	else if(!(visibleL.endP == 0 && visibleL.beginP == (pengingLth - 1)))
	{
		if (visibleL.beginP > visibleL.endP)
		{
			int pn = visibleL.beginP;
			visibleL.beginP = visibleL.endP;
			visibleL.endP = pn;
		}
	}


	tempLine = allLines[pendingPol.objectInAllLines[visibleL.beginP]];
	tempLine.PartSegmentation(visibleL.u, partLines);
	if (partLines.size() == 0)
	{
		std::cout << "δ�ɹ��ָ�����" << std::endl;
		//system("pause");
		return false;
	}
	int cptLth1 = partLines[0].m_CtrlPts.size();
	Vec3 p1 = partLines[0].m_CtrlPts[0];
	Vec3 p2 = partLines[0].m_CtrlPts[cptLth1 - 1];

	if (!(JudgeTwoPointsCoincide(pendingPol.p_Points[visibleL.beginP].vetx, p1) || JudgeTwoPointsCoincide(pendingPol.p_Points[visibleL.beginP].vetx, p2)))
	{
		tempLine = partLines[0];
		partLines[0] = partLines[1];
		partLines[1] = tempLine;
	}

	cptLth1 = partLines[0].m_CtrlPts.size();
	p1 = partLines[0].m_CtrlPts[0];
	p2 = partLines[0].m_CtrlPts[cptLth1 - 1];
	Vec3 segPoint;
	if (JudgeTwoPointsCoincide(pendingPol.p_Points[visibleL.beginP].vetx, p1))
	{
		segPoint = p2;
	}
	else
	{
		segPoint = p1;
	}

	if (pendingPol.p_Points[visibleL.beginP].ifOutPoint == 1 && pendingPol.p_Points[visibleL.endP].ifOutPoint == 1)
	{
		tempSubP.ifOutPoint = 1;
	}
	else
	{
		tempSubP.ifOutPoint = -1; //�趨Ϊ�Ȳ�������������Ҳ����ԭ����������
	}	
	tempSubP.vetx = segPoint;

	allLines.push_back(partLines[0]);
	allLines.push_back(partLines[1]);

	//��ȡ�µ�����NURBS����
	tempLine = GetPartNurbsLine(pendingPol.p_Points[visibleL.pt].vetx, pendingPol.p_Points[visibleL.beginP].vetx, visibleL.u, allLines, pendingPol);

	allLines.push_back(tempLine);
	alLth = allLines.size();

	//��pt��beginP��endPλ�÷��������
	int tempFlag;
	if (visibleL.pt == 0)
	{
		tempFlag = 0;
	}
	else if (visibleL.beginP == 0)
	{
		tempFlag = 1;
	}
	else if (visibleL.endP == 0)
	{
		tempFlag = 2;
	}
	else
	{
		if (visibleL.pt > visibleL.endP)
		{
			tempFlag = 3;
		}
		else if (visibleL.pt < visibleL.beginP)
		{
			tempFlag = 4;
		}
	}

	tempPol.convexHull.clear();
	tempPol.objectInAllLines.clear();
	tempPol.p_NurbsLines.clear();
	tempPol.p_Points.clear();
	outputPol.clear();
	switch (tempFlag)
	{
	case 0:
		for (int i = 0; i < visibleL.beginP; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		tempPol.p_Points.push_back(pendingPol.p_Points[visibleL.beginP]);
		tempPol.objectInAllLines.push_back(alLth-3);
		tempPol.p_Points.push_back(tempSubP);
		tempPol.objectInAllLines.push_back(alLth - 1);
		outputPol.push_back(tempPol);

		tempPol.convexHull.clear();
		tempPol.objectInAllLines.clear();
		tempPol.p_NurbsLines.clear();
		tempPol.p_Points.clear();

		//�ڶ������������
		tempPol.p_Points.push_back(pendingPol.p_Points[0]);
		tempPol.objectInAllLines.push_back(alLth - 1);
		tempPol.p_Points.push_back(tempSubP);
		tempPol.objectInAllLines.push_back(alLth - 2);
		for (int i = visibleL.endP; i < pengingLth; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		outputPol.push_back(tempPol);
		break;

	case 1:
		tempPol.p_Points.push_back(pendingPol.p_Points[0]);
		tempPol.objectInAllLines.push_back(alLth - 3);
		tempPol.p_Points.push_back(tempSubP);
		tempPol.objectInAllLines.push_back(alLth - 1);
		for (int i = visibleL.pt; i < pengingLth; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		outputPol.push_back(tempPol);

		tempPol.convexHull.clear();
		tempPol.objectInAllLines.clear();
		tempPol.p_NurbsLines.clear();
		tempPol.p_Points.clear();

		tempPol.p_Points.push_back(tempSubP);
		tempPol.objectInAllLines.push_back(alLth - 2);
		for (int i = visibleL.endP; i < visibleL.pt; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		tempPol.p_Points.push_back(pendingPol.p_Points[visibleL.pt]);
		tempPol.objectInAllLines.push_back(alLth - 1);
		outputPol.push_back(tempPol);
		break;

	case 2:
		tempPol.p_Points.push_back(tempSubP);
		tempPol.objectInAllLines.push_back(alLth - 1);
		for (int i = visibleL.pt; i < pengingLth - 1; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		tempPol.p_Points.push_back(pendingPol.p_Points[visibleL.beginP]);
		tempPol.objectInAllLines.push_back(alLth - 3);
		outputPol.push_back(tempPol);

		tempPol.convexHull.clear();
		tempPol.objectInAllLines.clear();
		tempPol.p_NurbsLines.clear();
		tempPol.p_Points.clear();

		tempPol.p_Points.push_back(tempSubP);
		tempPol.objectInAllLines.push_back(alLth - 2);
		for (int i = visibleL.endP; i < visibleL.pt; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		tempPol.p_Points.push_back(pendingPol.p_Points[visibleL.pt]);
		tempPol.objectInAllLines.push_back(alLth - 1);
		outputPol.push_back(tempPol);
		break;

	case 3:
		for (int i = 0; i < visibleL.beginP; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		tempPol.p_Points.push_back(pendingPol.p_Points[visibleL.beginP]);
		tempPol.objectInAllLines.push_back(alLth - 3);
		tempPol.p_Points.push_back(tempSubP);
		tempPol.objectInAllLines.push_back(alLth - 1);
		for (int i = visibleL.pt; i < pengingLth; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		outputPol.push_back(tempPol);

		tempPol.convexHull.clear();
		tempPol.objectInAllLines.clear();
		tempPol.p_NurbsLines.clear();
		tempPol.p_Points.clear();

		tempPol.p_Points.push_back(tempSubP);
		tempPol.objectInAllLines.push_back(alLth - 2);
		for (int i = visibleL.endP; i < visibleL.pt; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		tempPol.p_Points.push_back(pendingPol.p_Points[visibleL.pt]);
		tempPol.objectInAllLines.push_back(alLth - 1);
		outputPol.push_back(tempPol);
		break;

	case 4:
		for (int i = 0; i < visibleL.pt; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		tempPol.p_Points.push_back(pendingPol.p_Points[visibleL.pt]);
		tempPol.objectInAllLines.push_back(alLth - 1);
		tempPol.p_Points.push_back(tempSubP);
		tempPol.objectInAllLines.push_back(alLth - 2);
		for (int i = visibleL.endP; i < pengingLth; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		outputPol.push_back(tempPol);
		
		tempPol.convexHull.clear();
		tempPol.objectInAllLines.clear();
		tempPol.p_NurbsLines.clear();
		tempPol.p_Points.clear();

		for (int i = visibleL.pt; i < visibleL.beginP; i++)
		{
			tempPol.p_Points.push_back(pendingPol.p_Points[i]);
			tempPol.objectInAllLines.push_back(pendingPol.objectInAllLines[i]);
		}
		tempPol.p_Points.push_back(pendingPol.p_Points[visibleL.beginP]);
		tempPol.objectInAllLines.push_back(alLth - 3);
		tempPol.p_Points.push_back(tempSubP);
		tempPol.objectInAllLines.push_back(alLth - 1);
		outputPol.push_back(tempPol);
		break;

	default:
		break;
	}

	//ɾ�����һ�������
	morePol.pop_back();

	//�ж������Ķ���ε����ͱ����Ƿ�һ��
	for (const auto &pol : outputPol) {
		if (pol.p_Points.size() != pol.objectInAllLines.size())
			system("pause");
	}

	//�޸Ķ�Ӧ����һ�����������(�������)
	mLth = morePol.size();
	if (tempSubP.ifOutPoint != 1 && mLth!=0)
	{
		int lNum = pendingPol.objectInAllLines[visibleL.beginP];
		int objPolCondition = -1;
		int pNum = -1;
		MyPolygon tempPol2;
		for (int i = 0; i < mLth; i++)
		{
			int lth = morePol[i].objectInAllLines.size();
			for (int j = 0; j < lth; j++)
			{
				if (morePol[i].objectInAllLines[j] == lNum)
				{
					objPolCondition = i;
					pNum = j;
					tempPol2 = morePol[i];
					break;
				}
			}
			if (objPolCondition != -1)
			{
				break;
			}
		}
		if(pNum>-1)
		{
			if (pNum == tempPol2.objectInAllLines.size() - 1)
			{
				//ɾ�����һ����
				tempPol2.objectInAllLines.pop_back();
				if (JudgeTwoPointsCoincide(tempPol2.p_Points[tempPol2.p_Points.size() - 1].vetx, pendingPol.p_Points[visibleL.beginP].vetx))
				{
					tempPol2.p_Points.push_back(tempSubP);
					tempPol2.objectInAllLines.push_back(alLth - 3);
					tempPol2.objectInAllLines.push_back(alLth - 2);
				}
				else
				{
					tempPol2.p_Points.push_back(tempSubP);
					tempPol2.objectInAllLines.push_back(alLth - 2);
					tempPol2.objectInAllLines.push_back(alLth - 3);
				}
			}
			else
			{
				if (JudgeTwoPointsCoincide(tempPol2.p_Points[pNum].vetx, pendingPol.p_Points[visibleL.beginP].vetx))
				{
					tempPol2.p_Points.insert(&tempPol2.p_Points[pNum + 1], tempSubP);
					tempPol2.objectInAllLines[pNum] = alLth - 3;
					tempPol2.objectInAllLines.insert(&tempPol2.objectInAllLines[pNum + 1], alLth - 2);
				}
				else
				{
					tempPol2.p_Points.insert(&tempPol2.p_Points[pNum + 1], tempSubP);
					tempPol2.objectInAllLines[pNum] = alLth - 2;
					tempPol2.objectInAllLines.insert(&tempPol2.objectInAllLines[pNum + 1], alLth - 3);
				}
			}
			morePol[objPolCondition] = tempPol2;
		}
		
	}
	
	//�ж��ʷֳ��Ķ�����Ƿ�����ı���
	for (int i = 0; i < 2; i++)
	{
		int pLth = outputPol[i].p_Points.size();
		if (pLth == 4)
		{
			JudgeConcavePoint(outputPol[i]);
			for (int j = 0; j < 4; j++)
			{
				allLines[outputPol[i].objectInAllLines[j]].ifSegLine = 0;
				if (outputPol[i].convexHull[j] == 1)
				{
					return false;
				}
			}
			if (largeAngle) {
				//�ж��Ƿ��нǶȹ���
				Vec3 vec1, vec2;
				for (int j = 0; j < 4; j++) {
					GetVecs(outputPol[i], j, vec1, vec2);
					double angle = GetAngle(vec2, vec1);
					if (angle > lag)
						return false;
				}
			}
			quadPol.push_back(outputPol[i]);
		}
		else
		{
			morePol.push_back(outputPol[i]);
		}
	}
	return true;

}


bool QuadPart::QuadOperation(const varray<MyPolygon>& quadPolygen, const varray<MyPolygon>& morePolygen, const varray<PartNurbsLine>& allLines, varray<MyPolygon>& allPolygen, varray<PartNurbsLine>& endLines, int mode)
{
	cout << "QuadOperation()�ʷֲ�����......" << endl;
	bool flag = false;
	/*varray<MyPolygon> Quad, More;
	varray<PartNurbsLine> AllLines;
	Quad = quadPolygen;
	More = morePolygen;
	AllLines = allLines;*/
	int moreLth = morePolygen.size();
	if (moreLth == 0)
	{
		cout << "ȫ���ʷ����" << endl;
		endLines = allLines;
		allPolygen = quadPolygen;
		return true;
	}
	if (moreLth == 1)
	{
		cout << "ֻʣһ�������δ�ʷ�" << endl;
	}

		
	//ȡ���һ������ν��д���
	MyPolygon tempPolygen = morePolygen[moreLth - 1];
	int pnumber = tempPolygen.p_Points.size();
	
	//������(���ҽ���ֻ����һ������������ʱ)
	if (pnumber == 3) {
		cout << "���һ�������Ϊ�����Σ������ʷ�..." << endl;
		assert(quadPolygen.size() == 0);
		assert(moreLth == 1);
		Vec3 gravity = CalGravity(tempPolygen.p_Points);	//�����������������
		varray<Spline> templines;
		varray<Vec3> newPoints;//�ָ�������������
		varray<double> U;//�ضϴ��Ĳ���uֵ
		for (auto i : tempPolygen.objectInAllLines) {
			//������߽ضϴ����꼰uֵ
			double u;
			Vec3 p;
			//�õ������
			Vec3 p2 = GetFootPerpendicular(allLines[i].m_CtrlPts[0], *(allLines[i].m_CtrlPts.end() - 1), gravity);
			if (JudegeFootPerpendicular(allLines[i].m_CtrlPts[0], *(allLines[i].m_CtrlPts.end() - 1), p2)) {
				//�õ���������u
				u = GetPartU(allLines[i], allLines[i].m_CtrlPts[0], *(allLines[i].m_CtrlPts.end() - 1), p2);
			}
			else {
				//����ȡ0.5
				u = 0.5;
			}
			U.push_back(u);
			p= allLines[i].GetLinePoint(u);//�õ��ڵ�u���ĵ�����
			newPoints.push_back(p);

			PartNurbsLine tmpline;//���������ĵ����е�����
			CreatPartNurbsLine(gravity, p, tmpline);//���ӳ���
			tmpline.ifSegLine = 0;
			endLines.push_back(tmpline);
		}
		for (int i = 0; i < tempPolygen.objectInAllLines.size();i++) {
			//���������ضϺ������
			varray<PartNurbsLine> tmpline;
			PartNurbsLine curLine = allLines[tempPolygen.objectInAllLines[i]];
			curLine.PartSegmentation(U[i], tmpline);//�ض��ߴ���
			assert(tmpline.size() == 2);
			tmpline[0].ifSegLine = 0;//��Ϊ���ɷָ�
			tmpline[1].ifSegLine = 0;
			//����ضϺ�����
			if (JudgeTwoPointsCoincide(tmpline[0].m_CtrlPts[0], tempPolygen.p_Points[i].vetx) ||
				JudgeTwoPointsCoincide(*(tmpline[0].m_CtrlPts.end() - 1), tempPolygen.p_Points[i].vetx)) {
				endLines.push_back(tmpline[0]);
				endLines.push_back(tmpline[1]);
			}
			else
			{
				endLines.push_back(tmpline[1]);
				endLines.push_back(tmpline[0]);
			}
		}
		MyPolygon tmpPol;
		//�����һ���ı���
		tmpPol.p_Points.clear();
		tmpPol.objectInAllLines.clear();
		tmpPol.p_Points.push_back(SubPoint(1,gravity));
		tmpPol.p_Points.push_back(SubPoint(1, newPoints[0]));
		tmpPol.p_Points.push_back(SubPoint(1, tempPolygen.p_Points[1].vetx));
		tmpPol.p_Points.push_back(SubPoint(1, newPoints[1]));

		tmpPol.objectInAllLines.push_back(0);
		tmpPol.objectInAllLines.push_back(4);
		tmpPol.objectInAllLines.push_back(5);
		tmpPol.objectInAllLines.push_back(1);

		allPolygen.push_back(tmpPol);

		//����ڶ����ı���
		tmpPol.p_Points.clear();
		tmpPol.objectInAllLines.clear();
		tmpPol.p_Points.push_back(SubPoint(1, gravity));
		tmpPol.p_Points.push_back(SubPoint(1, newPoints[1]));
		tmpPol.p_Points.push_back(SubPoint(1, tempPolygen.p_Points[2].vetx));
		tmpPol.p_Points.push_back(SubPoint(1, newPoints[2]));

		tmpPol.objectInAllLines.push_back(1);
		tmpPol.objectInAllLines.push_back(6);
		tmpPol.objectInAllLines.push_back(7);
		tmpPol.objectInAllLines.push_back(2);

		allPolygen.push_back(tmpPol);

		//����������ı���
		tmpPol.p_Points.clear();
		tmpPol.objectInAllLines.clear();
		tmpPol.p_Points.push_back(SubPoint(1, gravity));
		tmpPol.p_Points.push_back(SubPoint(1, newPoints[2]));
		tmpPol.p_Points.push_back(SubPoint(1, tempPolygen.p_Points[0].vetx));
		tmpPol.p_Points.push_back(SubPoint(1, newPoints[0]));

		tmpPol.objectInAllLines.push_back(2);
		tmpPol.objectInAllLines.push_back(8);
		tmpPol.objectInAllLines.push_back(3);
		tmpPol.objectInAllLines.push_back(0);

		allPolygen.push_back(tmpPol);

		return true;
	}

	//�ı���
	if (!preProcessQuadPol&&pnumber == 4) {
		cout << "���һ�������Ϊ�ı���,��û����ǰ������ȡ���ı߼����У�������д���..." << endl;
		varray<MyPolygon> Quad, More;
		varray<PartNurbsLine> AllLines;
		Quad = quadPolygen;
		More = morePolygen;
		AllLines = allLines;
		//���ж��Ƿ���ڰ���
		JudgeConcavePoint(tempPolygen);
		for (int j = 0; j < 4; j++)
		{
			if (tempPolygen.convexHull[j] == 1)
			{
				return false;
			}
		}
		if (largeAngle) {
			//�ж��Ƿ��нǶȹ���
			Vec3 vec1, vec2;
			for (int j = 0; j < 4; j++) {
				GetVecs(tempPolygen, j, vec1, vec2);
				double angle = GetAngle(vec2, vec1);
				if (angle > lag)
					return false;
			}
		}
		for (int j = 0; j < tempPolygen.objectInAllLines.size(); j++) {
			AllLines[tempPolygen.objectInAllLines[j]].ifSegLine = 0;
		}
		Quad.push_back(tempPolygen);
		More.pop_back();
		flag = QuadOperation(Quad, More, AllLines, allPolygen, endLines, mode);//�ݹ�
		return flag;
	}

	//�����μ�����
	if (pnumber > 5)
	{
		cout << "�����μ����϶���ν����ʷ�..." << endl;
		varray<VisiblePoints> vsp, svsp;	//�޼��ʹ��ڼ�����
		varray<VisibleLine> partLine;		//�ض�ĳһ���ߵ��������
		JudgeConcavePoint(tempPolygen);		//�ж��Ƿ���ڰ���
		CalPartLines(tempPolygen, allLines, partLine, mode);
		FindVisiblePoint(tempPolygen, vsp);
		VisiblePointConcaveJudge(tempPolygen, vsp);
		VisiblePointSharpJudeg(tempPolygen, vsp, allLines);
		SortVisiblePoints(tempPolygen, vsp, svsp);
		if (mode == 1)
		{
			svsp.clear();
		}
		int L1, L2, L3;
		L1 = vsp.size();
		L2 = partLine.size();
		L3 = svsp.size();

		if (byW) {
			//��vsp��partLine����ͳһ����ѡ��
			int i = 0, j = 0, k = 0;
			vector<int> seq(L1 + L2);		//0��ʾ����L1����,1��ʾ����L2����
			while (i < L1 &&j < L2) {
				if (vsp[i].w >= partLine[j].w) {
					seq[k] = 0;
					++i;
				}
				else {
					seq[k] = 1;
					++j;
				}
				++k;
			}
			while (i < L1) {
				seq[k++] = 0;
				++i;
			}
			while (j < L2) {
				seq[k++] = 1;
				++j;
			}
			assert(k == seq.size());

			//����seq���������ʷ�
			i = 0, j = 0;
			for (k = 0; k < seq.size(); ++k) {
				if (seq[k] == 0) {
					//���������
					varray<MyPolygon> Quad, More;
					varray<PartNurbsLine> AllLines;
					Quad = quadPolygen;
					More = morePolygen;
					AllLines = allLines;
					bool tempFlag = ConncetTwoPoints(tempPolygen, Quad, More, AllLines, vsp[i]);
					if (!tempFlag)
					{
						++i;
						continue;
					}
					else
					{
						flag = QuadOperation(Quad, More, AllLines, allPolygen, endLines, mode);
						if (flag)
						{
							return true;
						}
						else
						{
							++i;
							continue;
						}
					}
				}

				else if (seq[k] == 1) {
					//�ض�ʽ����
					varray<MyPolygon> Quad, More;
					varray<PartNurbsLine> AllLines;
					Quad = quadPolygen;
					More = morePolygen;
					AllLines = allLines;
					bool tempFlag = ConnectLineWithPoint(tempPolygen, Quad, More, AllLines, partLine[j]);
					if (!tempFlag)
					{
						++j;
						continue;
					}
					else
					{
						flag = QuadOperation(Quad, More, AllLines, allPolygen, endLines, mode);
						if (flag)
						{
							return true;
						}
						else
						{
							++j;
							continue;
						}
					}
				}
			}
		}

		else {
			//vsp��partLine���Խ��а������
			if (L1 != 0)
			{
				for (int i = 0; i < L1; i++)
				{
					varray<MyPolygon> Quad, More;
					varray<PartNurbsLine> AllLines;
					Quad = quadPolygen;
					More = morePolygen;
					AllLines = allLines;


					bool tempFlag = ConncetTwoPoints(tempPolygen, Quad, More, AllLines, vsp[i]);
					if (!tempFlag)
					{
						continue;
					}
					else
					{
						flag = QuadOperation(Quad, More, AllLines, allPolygen, endLines, mode);
						if (flag)
						{
							return true;
						}
						else
						{
							continue;
						}
					}
				}
			}

			if (L2 != 0)
			{
				for (int i = 0; i < L2; i++)
				{
					varray<MyPolygon> Quad, More;
					varray<PartNurbsLine> AllLines;
					Quad = quadPolygen;
					More = morePolygen;
					AllLines = allLines;
					bool tempFlag = ConnectLineWithPoint(tempPolygen, Quad, More, AllLines, partLine[i]);
					if (!tempFlag)
					{
						continue;
					}
					else
					{
						flag = QuadOperation(Quad, More, AllLines, allPolygen, endLines, mode);
						if (flag)
						{
							return true;
						}
						else
						{
							continue;
						}
					}
				}
			}
		}
		

		//������svsp��㼯�ϵı���
		if (L3 != 0)
		{
			for (int i = 0; i < L3; i++)
			{
				varray<MyPolygon> Quad, More;
				varray<PartNurbsLine> AllLines;
				Quad = quadPolygen;
				More = morePolygen;
				AllLines = allLines;
				bool tempFlag = ConncetTwoPoints(tempPolygen, Quad, More, AllLines, svsp[i]);
				if (!tempFlag)
				{
					continue;
				}
				else
				{
					flag = QuadOperation(Quad, More, AllLines, allPolygen, endLines, mode);
					if (flag)
					{
						return true;
					}
					else
					{
						continue;
					}
				}
			}
		}	
	}
	
	//�����
	if (pnumber == 5)
	{
		cout << "����ν����ʷ�..." << endl;
		for (int i = 0; i < 5; i++)
		{
			if (allLines[tempPolygen.objectInAllLines[i]].ifSegLine == 1)
			{
				break;
			}
			if (i == 4)
			{
				return false;
			}	
		}
		varray<VisibleLine> partLine;
		JudgeConcavePoint(tempPolygen);
		CalPartLines(tempPolygen, allLines, partLine, mode);
		if (partLine.size() == 0)
		{
			return false;
		}
		else
		{
			for (int i = 0; i < partLine.size(); i++)
			{
				varray<MyPolygon> Quad, More;
				varray<PartNurbsLine> AllLines;
				Quad = quadPolygen;
				More = morePolygen;
				AllLines = allLines;
				bool tempFlag = ConnectLineWithPoint(tempPolygen, Quad, More, AllLines, partLine[i]);
				if (!tempFlag)
				{
					continue;
				}
				else
				{
					flag = QuadOperation(Quad, More, AllLines, allPolygen, endLines, mode);
					if (flag)
					{
						return true;
					}
					else
					{
						continue;
					}
				}
			}
		}
	}

	return false;
}




//�������߲����ÿɷָ���
void PartNurbsLine::CreatPartNurbsLine(Spline inputLine, int ifSeg)
{
	m_Degree = inputLine.m_Degree;
	m_CtrlPts = inputLine.m_CtrlPts;
	m_Knots = inputLine.m_Knots;
	ifSegLine = ifSeg;
}

//����Nurbsline����
void PartNurbsLine::CreatNurbsLine(Spline & newLine)
{
	newLine.m_CtrlPts = m_CtrlPts;
	newLine.m_Degree = m_Degree;
	newLine.m_Knots = m_Knots;
	return;
}

void PartNurbsLine::PartSegmentation(const double u, varray<PartNurbsLine>& lines)
{
	lines.clear();
	if (u <= 0 || u >= 1)
		return;

	lines.resize(2);

	Spline l0(*this);
	int insN = m_Degree;
	int span = 0;
	if (IsInSet(u, m_Knots))
	{
		span = FindSpan(u, m_Degree, m_CtrlPts.size(), m_Knots);
		for (int i = span; i > 0; --i)
		{
			if (u == m_Knots[i])
				--insN;
			else
				break;
		}
	}
	if (insN > 0)
		l0.KnotInsert(u, insN);
	lines[0].m_Degree = l0.m_Degree;
	lines[1].m_Degree = l0.m_Degree;
	//�����½ڵ�ʸ��
	for (int i = 0; i < l0.m_Degree + 1; ++i)
	{
		lines[0].m_Knots.push_back(0);
		lines[1].m_Knots.push_back(0);
	}
	span = l0.FindSpan(u, l0.m_Degree, l0.m_CtrlPts.size(), l0.m_Knots);
	for (int i = l0.m_Degree + 1; i <= span - l0.m_Degree; ++i)
	{
		double newu = l0.m_Knots[i] / u;
		lines[0].m_Knots.push_back(newu);
	}
	for (int i = span + 1; i < l0.m_CtrlPts.size(); ++i)
	{
		double newu = 1 - (1 - l0.m_Knots[i]) / (1 - u);
		lines[1].m_Knots.push_back(newu);
	}
	for (int i = 0; i < l0.m_Degree + 1; ++i)
	{
		lines[0].m_Knots.push_back(1);
		lines[1].m_Knots.push_back(1);
	}
	//���Ƶ�
	int i = 0;
	for (i; i < lines[0].m_Knots.size() - lines[0].m_Degree - 1; ++i)
		lines[0].m_CtrlPts.push_back(l0.m_CtrlPts[i]);
	for (--i; i < l0.m_CtrlPts.size(); ++i)
		lines[1].m_CtrlPts.push_back(l0.m_CtrlPts[i]);
	lines[0].ifSegLine = 1;
	lines[1].ifSegLine = 1;
}

Vec3 PartNurbsLine::CalBeginDirecVec()
{
	Vec3 vec;
	vec = Vec3(m_CtrlPts[1].x - m_CtrlPts[0].x, m_CtrlPts[1].y - m_CtrlPts[0].y, m_CtrlPts[1].z - m_CtrlPts[0].z);
	return vec;
}

Vec3 PartNurbsLine::CalBeginDirecVec() const
{
	Vec3 vec;
	vec = Vec3(m_CtrlPts[1].x - m_CtrlPts[0].x, m_CtrlPts[1].y - m_CtrlPts[0].y, m_CtrlPts[1].z - m_CtrlPts[0].z);
	return vec;
}

Vec3 PartNurbsLine::CalEndDirecVec()
{
	Vec3 vec;
	int ctpLength = this->m_CtrlPts.size();
	vec = Vec3(m_CtrlPts[ctpLength - 2].x - m_CtrlPts[ctpLength - 1].x, m_CtrlPts[ctpLength - 2].y - m_CtrlPts[ctpLength - 1].y, m_CtrlPts[ctpLength - 2].z - m_CtrlPts[ctpLength - 1].z);
	return vec;
}

Vec3 PartNurbsLine::CalEndDirecVec() const
{
	Vec3 vec;
	int ctpLength = this->m_CtrlPts.size();
	vec = Vec3(m_CtrlPts[ctpLength - 2].x - m_CtrlPts[ctpLength - 1].x, m_CtrlPts[ctpLength - 2].y - m_CtrlPts[ctpLength - 1].y, m_CtrlPts[ctpLength - 2].z - m_CtrlPts[ctpLength - 1].z);
	return vec;
}

//�õ������
Vec3 GetFootPerpendicular(const Vec3 & begin, const Vec3 & end, const Vec3 &pt)
{
	double length, cosSita;
	Vec3  retVal,vec1, vec2;
	//vec1 = end - begin;
	//length = vec1.Magnitude();
	/*pt = begin + end;
	pt.x = pt.x*0.5;
	pt.y = pt.y*0.5;
	pt.z = pt.z*0.5;

	pt.x += length;
	vec2 = pt - begin;
	vec1 = vec1.Normalize();
	vec2 = vec2.Normalize();
	cosSita = vec1.Dot(vec2);
	if (abs(cosSita - 1) < 1e-8 || abs(cosSita + 1) < 1e-8)
	{
		pt.x -= length;
		pt.y += length;
	}*/

	double dx = begin.x - end.x;
	double dy = begin.y - end.y;
	double dz = begin.z - end.z;

	//�غ�
	if (abs(dx) < 1e-8 && abs(dy) < 1e-8 && abs(dz) < 1e-8)
	{
		retVal = begin;
		return retVal;
	}

	double u = (pt.x - begin.x)*(begin.x - end.x) +
		(pt.y - begin.y)*(begin.y - end.y) + (pt.z - begin.z)*(begin.z - end.z);
	u = u / ((dx*dx) + (dy*dy) + (dz*dz));
	retVal.x = begin.x + u * dx;
	retVal.y = begin.y + u * dy;
	retVal.z = begin.z + u * dz;

	return retVal;
}

//�ж��Ƿ�ֱ
bool JudegeFootPerpendicular(const Vec3 & begin, const Vec3 & end, const Vec3 &crossPoint)
{
	if (isinf(crossPoint.x) || isinf(crossPoint.y)|| isnan(crossPoint.x) || isnan(crossPoint.y)) {
		//�ⲿ���ж�������
		return false;
	}
	Vec3 vec1, vec2;
	vec1 = crossPoint - begin;
	vec2 = crossPoint - end;
	if (abs(vec1.x - 0) < 1e-8 && abs(vec1.y - 0) < 1e-8 && abs(vec1.z - 0) < 1e-8)
		return false;
	if (abs(vec2.x - 0) < 1e-8 && abs(vec2.y - 0) < 1e-8 && abs(vec2.z - 0) < 1e-8)
		return false;
	vec1 = vec1.Normalize();
	vec2 = vec2.Normalize();
	double cosSita = vec1.Dot(vec2);

	if (abs(cosSita - 0) < 1e-8)
		return false;
	if (cosSita + 1 < 1e-8)
		return true;
	return false;//��ʲô���������,����false
}
//�õ���ƽ�����е�����
Vec3 GetAngularBisectorVec(const Vec3 & leftVec, const Vec3 & rigthtVec)
{
	Vec3 vec1, vec2, middleVec;
	double leftLength, rightLength;
	leftLength = leftVec.Magnitude();
	rightLength = rigthtVec.Magnitude();
	vec1 = leftVec * rightLength;
	vec2 = rigthtVec * leftLength;
	middleVec = vec1 + vec2;
	middleVec = middleVec.Normalize();
	if (middleVec.Magnitude() - 0 < 1e-8)
	{
		middleVec = Vec3(rigthtVec.y*(-1.0), rigthtVec.x, rigthtVec.z);
		middleVec = middleVec.Normalize();
	}
	else
	{
		double angle = GetAngle(vec2, vec1);
		if (angle > PI)
		{
			middleVec = -middleVec;
		}
	}
	return middleVec;
}

Vec3 GetAngularBisectorPoint(const Vec3 & begin, const Vec3 & end, const Vec3 & pt, const Vec3 & vec)
{
	Vec3 pt2 = pt + vec;
	Vec3 crossp= GetTwoLinesCrossPoint(begin, end, pt, pt2);
	return crossp;
}

void SurfaceConverse(varray<Spline>& lines, Vec3 targetPts, Vec3 targetVec, Matrix4d & conMx, Matrix4d & reconMx)
{
	if (lines.size() < 2)
	{
		cout << "�ռ�任����������������" << endl;
		return;
	}
	else
	{
		varray<Spline> tmpLines;
		Vec3 vec1, vec2, norVec; 
		vec1 = lines[0].m_CtrlPts[0] - *(lines[0].m_CtrlPts.end()-1);
		vec1 = vec1.Normalize();
		int beg = 1;
		while (beg < lines.size())
		{
			vec2= lines[beg].m_CtrlPts[0] - *(lines[beg].m_CtrlPts.end() - 1);
			vec2 = vec2.Normalize();
			double cosst = vec1.Dot(vec2);
			if (abs(cosst - 1) > 1e-4)
			{
				//������ķ�����
				norVec = vec1.Cross(vec2);
				norVec = norVec.Normalize();
				break;
			}
			if (beg == lines.size() - 1)
			{
				cout << "δ��ú��ʷ�����" << endl;
			}
		}
		CalTransMat(targetPts, targetVec, lines[0].m_CtrlPts[0], norVec, conMx);
		reconMx = conMx.inverse();
		for (int i = 0; i < lines.size(); i++)
		{
			TransByMat(lines[i].m_CtrlPts, conMx);
		}
	}
}

void LineConverse(varray<Spline>& lines, Matrix4d & conMx)
{
	for (int i = 0; i < lines.size(); i++)
	{
		TransByMat(lines[i].m_CtrlPts, conMx);
	}
}

int CalU(double u)
{
	if (fabs(u - 0) < 1e-5)
	{
		return 0;
	}
	else if (fabs(u - 1) < 1e-5)
	{
		return 1;
	}
	return -1;
}

varray<Vec3> OrderLinesAntioclock(varray<Spline> &polLines)
{
	int line1 = -1, line2 = -1;										//��ѡȡ�ĵ�һ���������
	Vec3 tmpPoint;
	Vec3 secondPoint1, secondPoint2;								//��ѡȡ�ĵڶ����㣬����ѡȡ��Ҫ��
	varray<Vec3> orderedPoints;
	if (polLines.size() < 4)										//���������߼�����������4������ӡerror������
	{
		cout << "OrderLinesAntioclock��LinesSize<4" << endl;
		//system("pause");
		if (polLines.size() == 1) {
			cout << "OrderLinesAntioclock-Only One Line" << endl;
			Vec3 p = (Vec3)polLines[0].m_CtrlPts[0];
			orderedPoints.push_back(p);
			p = (Vec3)(*(polLines[0].m_CtrlPts.end() - 1));
			orderedPoints.push_back(p);
		}
		return orderedPoints;
	}
	
	tmpPoint = polLines[0].m_CtrlPts[0];
	for (auto i : polLines)											//�ҵ�y����ֵ��С�ĵ���Ϊ��ʼ��
	{
		if (i.m_CtrlPts[0].y < tmpPoint.y)
		{
			tmpPoint = i.m_CtrlPts[0];
		}
		if ((i.m_CtrlPts.end() - 1)->y < tmpPoint.y)
		{
			tmpPoint = *(i.m_CtrlPts.end() - 1);
		}
	}
	orderedPoints.push_back(tmpPoint);								//����ʼ���������õļ���

	for (int i = 0; i < polLines.size();i++)						//�ҵ���ʼ������������������
	{
		if (JudgeTwoPointsCoincide(tmpPoint, polLines[i].m_CtrlPts[0]))
		{
			if (line1 == -1 && line2 == -1)
			{
				secondPoint1 = *(polLines[i].m_CtrlPts.end() - 1);
				line1 = i;
			}
			else
			{
				secondPoint2 = *(polLines[i].m_CtrlPts.end() - 1);
				line2 = i;
			}
		}
		else if (JudgeTwoPointsCoincide(tmpPoint, *(polLines[i].m_CtrlPts.end() - 1)))
		{
			if (line1 == -1 && line2 == -1)
			{
				secondPoint1 = polLines[i].m_CtrlPts[0];
				line1 = i;
			}
			else
			{
				secondPoint2 = polLines[i].m_CtrlPts[0];
				line2 = i;
			}
		}
		if (line1 != -1 && line2 != -1)
		{
			break;
		}
	}

	Vec3 vec1 = secondPoint1 - tmpPoint;							//��ʼ����������1
	Vec3 vec2 = secondPoint2 - tmpPoint;							//��ʼ����������2
	vec1 = vec1.Normalize();
	vec2 = vec2.Normalize();
	Vec3 cross = vec1.Cross(vec2);//���������
	if (cross.z > 0)												//secondPoint1Ϊ�����
	{
		tmpPoint = secondPoint1;
		orderedPoints.push_back(tmpPoint);
		//��������˳��
		Spline tmpLine = polLines[line1];
		polLines[line1] = polLines[0];
		polLines[0] = tmpLine;
	}
	else if (cross.z < 0)											//secondPoint2Ϊ�����
	{
		tmpPoint = secondPoint2;
		orderedPoints.push_back(tmpPoint);
		Spline tmpLine = polLines[line2];
		polLines[line2] = polLines[0];
		polLines[0] = tmpLine;
	}
	else if (fabs(cross.z - 0)<1e-5)								//��������180��
	{
		if (vec1.x > vec2.x)
		{
			tmpPoint = secondPoint1;
			orderedPoints.push_back(tmpPoint);
			Spline tmpLine = polLines[line1];
			polLines[line1] = polLines[0];
			polLines[0] = tmpLine;
		}
		else
		{
			tmpPoint = secondPoint2;
			orderedPoints.push_back(tmpPoint);
			Spline tmpLine = polLines[line2];
			polLines[line2] = polLines[0];
			polLines[0] = tmpLine;
		}
	}

	//ѭ�������ʱ��˳���
	while (orderedPoints.size() < polLines.size())
	{
		for (int i = orderedPoints.size() - 1; i < polLines.size(); i++)
		{
			if (JudgeTwoPointsCoincide(tmpPoint, polLines[i].m_CtrlPts[0]))
			{
				tmpPoint = *(polLines[i].m_CtrlPts.end() - 1);
				//��������
				Spline tmpLine = polLines[i];
				polLines[i] = polLines[orderedPoints.size() - 1];
				polLines[orderedPoints.size() - 1] = tmpLine;
				//�����µĵ�
				orderedPoints.push_back(tmpPoint);
				break;
			}
			else if (JudgeTwoPointsCoincide(tmpPoint, *(polLines[i].m_CtrlPts.end() - 1)))
			{
				tmpPoint = polLines[i].m_CtrlPts[0];
				//��������
				Spline tmpLine = polLines[i];
				polLines[i] = polLines[orderedPoints.size() - 1];
				polLines[orderedPoints.size() - 1] = tmpLine;
				//�����µĵ�
				orderedPoints.push_back(tmpPoint);
				break;
			}
			if (i == polLines.size() - 1)
			{
				cout << "OrderLinesAntioclock��δ�ҵ���������" << endl;	//������������δ�ҵ���������һ������
				//system("pause");
			}
		}
	}

	//�������һ�����Ƿ���ȷ
	if (JudgeTwoPointsCoincide(orderedPoints[0], (polLines.end() - 1)->m_CtrlPts[0])
		|| JudgeTwoPointsCoincide(orderedPoints[0], *((polLines.end() - 1)->m_CtrlPts.end() - 1)))
	{
		cout << "���һ������Ϊ��Ӧ����" << endl;
	}

	return orderedPoints;
}

bool CalRightu(double u, double y, Spline Nline, vector<double>& us)
{
	const int N = 1000;//���߷ֳ�1000�߶Σ���ɢ��

	Vec3 temp = Nline.GetLinePoint(u);
	if (abs(temp.y - y) < 0.00001)
	{
		us.push_back(u);
	}


	auto function = [&](float begin1, float end2)
	{
		int Nn = 10000;
	lable:
		float a = abs(begin1 - end2);
		float step = a / Nn;
		double distance = 0.00001;
		vector<float> u;
		for (int i = 0; i != Nn; ++i)
		{
			Vec3 begin_point = Nline.GetLinePoint(begin1 + i * step);
			Vec3 end_point = Nline.GetLinePoint(begin1 + (i + 1)*step);
			if (abs(y - begin_point.y) <= distance)
			{
				u.push_back(begin1 + i * step);
			}
			else if (abs(y - end_point.y) <= distance)
			{
				u.push_back(begin1 + (i + 1) * step);
			}
		}
		if (u.size() == 0)
		{
			Nn = Nn * 10;
			goto lable;
		}
		else
		{
			us.push_back(u[0]);
		}
	};

	float step = 1.0 / N;
	float begin;
	begin = 0;
	for (int i = 0; i != N; ++i)
	{
		Vec3 begin_point = Nline.GetLinePoint(begin + i * step);
		Vec3 end_point = Nline.GetLinePoint(begin + (i + 1)*step);
		if (abs(y - begin_point.y) < 0.00001)
		{
			us.push_back(begin + i * step);
			continue;
		}
		if (abs(y - end_point.y) < 0.00001)
		{
			us.push_back(begin + (i + 1)*step);
			continue;
		}
		if ((y > begin_point.y&&end_point.y > y) || (y < begin_point.y&&end_point.y < y))
		{
			function(begin + i * step, begin + (i + 1)*step);
		}
	}

	//�������u����ȥ��
	for (vector<double>::iterator it1 = us.begin(); it1 != us.end();it1++) {
		for (vector<double>::iterator it2 = it1 + 1; it2 != us.end();) {
			if (abs(*it1 - *it2) < 1e-5) {
				it2 = us.erase(it2);
			}
			else
			{
				it2++;
			}
		}
	}

	if (us.size() == 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}

bool CalRightu(vector<double> u, double y, Spline Nline, vector<double>& us)
{
	const int N = 1000;//���߷ֳ�һ�ٸ��߶Σ���ɢ��
	for (int i = 0; i != u.size(); ++i)
	{
		Vec3 temp = Nline.GetLinePoint(u[i]);
		if (abs(temp.y - y) < 0.00001)
		{
			us.push_back(u[i]);
		}
	}

	auto function = [&](double begin1, double end2)
	{
		int Nn = 10000;
	lable:
		double a = abs(begin1 - end2);
		double step = a / Nn;
		double distance = 0.00001;
		vector<float> u;
		for (int i = 0; i != Nn; ++i)
		{
			Vec3 begin_point = Nline.GetLinePoint(begin1 + i * step);
			Vec3 end_point = Nline.GetLinePoint(begin1 + (i + 1)*step);
			if (abs(y - begin_point.y) <= distance)
			{
				u.push_back(begin1 + i * step);
			}
			else if (abs(y - end_point.y) <= distance)
			{
				u.push_back(begin1 + (i + 1) * step);
			}
		}
		if (u.size() == 0)
		{
			Nn = Nn * 10;
			goto lable;
		}
		else
		{
			us.push_back(u[0]);
		}
	};

	double step = 1.0 / N;
	double begin;
	begin = 0;
	for (int i = 0; i != N; ++i)
	{
		Vec3 begin_point = Nline.GetLinePoint(begin + i * step);
		Vec3 end_point = Nline.GetLinePoint(begin + (i + 1)*step);
		if (abs(y - begin_point.y) < 1e-10)
		{
			us.push_back(begin + i * step);
			continue;
		}
		if (abs(y - end_point.y) < 1e-10)
		{
			us.push_back(begin + (i + 1)*step);
			continue;
		}
		if ((y > begin_point.y&&end_point.y > y) || (y < begin_point.y&&end_point.y < y))
		{
			function(begin + i * step, begin + (i + 1)*step);
		}
	}

	//�������u����ȥ��
	for (vector<double>::iterator it1 = us.begin(); it1 != us.end(); it1++) {
		for (vector<double>::iterator it2 = it1 + 1; it2 != us.end();) {
			if (abs(*it1 - *it2) < 1e-5) {
				it2 = us.erase(it2);
			}
			else
			{
				it2++;
			}
		}
	}

	if (us.size() == 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}

bool CalRightu(vector<float> u, float y, Spline Nline, vector<float>& us)
{
	const int N = 1000;//���߷ֳ�һ�ٸ��߶Σ���ɢ��
	for (int i = 0; i != u.size(); ++i)
	{
		Vec3 temp = Nline.GetLinePoint(u[i]);
		if (abs(temp.y - y) < 0.00001)
		{
			us.push_back(u[i]);
		}
	}

	auto function = [&](float begin1, float end2)
	{
		int Nn = 10000;
	lable:
		float a = abs(begin1 - end2);
		float step = a / Nn;
		double distance = 0.00001;
		vector<float> u;
		for (int i = 0; i != Nn; ++i)
		{
			Vec3 begin_point = Nline.GetLinePoint(begin1 + i * step);
			Vec3 end_point = Nline.GetLinePoint(begin1 + (i + 1)*step);
			if (abs(y - begin_point.y) <= distance)
			{
				u.push_back(begin1 + i * step);
			}
			else if (abs(y - end_point.y) <= distance)
			{
				u.push_back(begin1 + (i + 1) * step);
			}
		}
		if (u.size() == 0)
		{
			Nn = Nn * 10;
			goto lable;
		}
		else
		{
			us.push_back(u[0]);
		}
	};

	float step = 1.0 / N;
	float begin;
	begin = 0;
	for (int i = 0; i != N; ++i)
	{
		Vec3 begin_point = Nline.GetLinePoint(begin + i * step);
		Vec3 end_point = Nline.GetLinePoint(begin + (i + 1)*step);
		if (abs(y - begin_point.y) < 0.00001)
		{
			us.push_back(begin + i * step);
			continue;
		}
		if (abs(y - end_point.y) < 0.00001)
		{
			us.push_back(begin + (i + 1)*step);
			continue;
		}
		if ((y > begin_point.y&&end_point.y > y) || (y < begin_point.y&&end_point.y < y))
		{
			function(begin + i * step, begin + (i + 1)*step);
		}
	}

	//�������u����ȥ��
	for (vector<float>::iterator it1 = us.begin(); it1 != us.end(); it1++) {
		for (vector<float>::iterator it2 = it1 + 1; it2 != us.end();) {
			if (abs(*it1 - *it2) < 1e-5) {
				it2 = us.erase(it2);
			}
			else
			{
				it2++;
			}
		}
	}

	if (us.size() == 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}

int LineIntersectNurbs(SISLCurve *& curve, const varray<Vec3>& strLine, const Spline * nurbsCurve)
{
	int count = 0;
	double *intpar = nullptr;			//�������ֵ
	int numintpt = 0;					//������
	int numintcu = 0;					//�غ�������
	vector<double> itp;					//��ȡ�Ľ������
	StraLineIntersectNurbsLine(curve, strLine, intpar, numintpt, numintcu);
	//�Կ��ܵĽ�������ж�
	if (numintpt != 0)
	{
		vector<double> tmpItp;
		itp.resize(numintpt);
		for (int i = 0; i < numintpt; i++) {
			itp[i] = intpar[i];
		}
		CalRightu(itp, strLine[0].y, *nurbsCurve, tmpItp);
		for (vector<double>::iterator it1 = itp.begin(); it1 != itp.end();) {
			bool flag = false;
			for (vector<double>::iterator it2 = tmpItp.begin(); it2 != tmpItp.end();) {
				if (abs(*it1 - *it2) < 1e-5) {
					it2 = tmpItp.erase(it2);
					flag = true;
					break;
				}
				it2++;
			}
			if (flag == false) {
				//��ʾԭ�����������ȷ
				it1 = itp.erase(it1);
			}
			else
				it1++;
		}
		for (int i = 0; i < tmpItp.size(); i++) {
			itp.push_back(tmpItp[i]);
		}
		numintpt = itp.size();

		//for (int i = 0; i < numintpt;) {
		//	if (CalU(intpar[i]) != -1)
		//	{
		//		//uΪ0��1
		//		continue;
		//		i++;
		//	}
		//	else {
		//		//u������֮��,��Ҫ�ж�
		//		double interVal = intpar[i];
		//		bool interPoint = false;
		//		interPoint = isInterPoint(interVal, curve, nurbsCurve);
		//		if (interPoint == true)
		//		{
		//			intpar[i] = interVal;
		//		}
		//		else
		//		{
		//			for (int j = i + 1; j < numintpt; j++)
		//			{
		//				intpar[j - 1] = intpar[j];
		//			}
		//			numintpt--;
		//		}
		//	}
		//}
		////ѭ���ж������ǲ������ظ���
		//for (int i = 0; i < numintpt-1; i++)
		//{
		//	for (int j = i + 1; j < numintpt; j++)
		//	{
		//		if (fabs(intpar[i] - intpar[j]) < 1e-20)
		//		{
		//			//��ʾ��������ͬһ������
		//			for (int k = j + 1; k < numintpt; k++)
		//			{
		//				intpar[k - 1] = intpar[k];
		//			}
		//			numintpt--;
		//		}
		//	}
		//}
	}

	if (numintpt != 0) {
		//���ڽ������
		//ֱ����itp���������ݽ��м���
		for (int i = 0; i < numintpt; i++) {
			double tmpU = itp[i];
			//double tmpU = intpar[i];	//���㴦�Ĳ���ֵ
			int flagU = CalU(tmpU);
			//������м�Ĳ���
			varray<Vec4> Der;//CNurbs�����ʽ���Ľ��������������
			if (0 == flagU) {
				//����Ϊ���߲���Ϊ0���ĵ�
				if (nurbsCurve) {
					//CNurbs�����ʽ�����߲�Ϊ��
					Der.clear();
					nurbsCurve->PtsDerivs(0, 1, Der);
				}
				if (Der.size() < 2) assert(0);
				if (Der[0].x < strLine[0][0]) continue;	//��������ߵ��ӳ�����

				if (fabs(Der[1].y - 0.0) < 1e-5) {
					//�ô�������ˮƽ
					Vec3 p = nurbsCurve->GetLinePoint(0.01);
					point2d tmpVec = point2d(p.x - Der[0].x, p.y - Der[0].y);
					tmpVec = tmpVec.Normalize();
					if (tmpVec.y < 0.0) {
						//��������3��4���ޱ�ʾ���õ���Ϊ����
						count++;
					}
				}
				else {
					//�ô���������ˮƽ
					//Vec3 p = nurbsCurve->GetLinePoint(0.05);
					//point2d tmpVec = point2d(p.x - Der[0].x, p.y - Der[0].y);
					//tmpVec = tmpVec.Normalize();
					//if (tmpVec.y < 0.0) {
					//	//��������3��4���ޱ�ʾ���õ���Ϊ����
					//	count++;
					//}
					if (Der[1].y < 0.0) {
						//�õ�����������3/4���ޣ�������Ϊ����
						count++;
					}
				}
			}

			else if (1 == flagU) {
				//����Ϊ���߲���Ϊ1���ĵ�
				if (nurbsCurve) {
					//CNurbs�����ʽ�����߲�Ϊ��
					Der.clear();
					nurbsCurve->PtsDerivs(1, 1, Der);
				}
				if (Der.size() < 2) assert(0);
				if (Der[0].x < strLine[0][0]) continue;

				if (abs(Der[1].y - 0.0) < 1e-5) {
					//�ô�������ˮƽ
					Vec3 p = nurbsCurve->GetLinePoint(0.99);
					point2d tmpVec = point2d(p.x - Der[0].x, p.y - Der[0].y);
					tmpVec = tmpVec.Normalize();
					if (tmpVec.y < 0.0) {
						//��������3��4���ޱ�ʾ���õ���Ϊ����
						count++;
					}
				}
				else {
					//�ô���������ˮƽ
					//Vec3 p = nurbsCurve->GetLinePoint(0.95);
					//point2d tmpVec = point2d(p.x - Der[0].x, p.y - Der[0].y);
					//tmpVec = tmpVec.Normalize();
					//if (tmpVec.y < 0.0) {
					//	//��������3��4���ޱ�ʾ���õ���Ϊ����
					//	count++;
					//}
					if (Der[1].y > 0.0) {
						//�õ�����������3/4���ޣ�������Ϊ����(��Ϊ�ǲ���Ϊ1����������Ӧ�÷���)
						count++;
					}
				}
			}

			else {
				//����Ϊ������һ��
				Der.clear();
				if (nurbsCurve) {
					//CNurbs�����ʽ�����߲�Ϊ��
					nurbsCurve->PtsDerivs(tmpU, 1, Der);
				}
				if (Der.size() < 2) assert(0);
				if (Der[0].x < strLine[0][0]) continue;

				if (abs(Der[1].y - 0.0) < 1e-5) {
					//�ô�������ˮƽ
					assert((tmpU - 0.01) > 0.0);
					assert((tmpU + 0.01) < 1.0);

					point2d tmpVec1, tmpVec2;
					//CNurbs��õĽ����������
					varray<Vec4> Der_1, Der_2;
					nurbsCurve->PtsDerivs(tmpU - 0.01, 1, Der_1);
					tmpVec1 = point2d(Der_1[0].x - Der[0].x, Der_1[0].y - Der[0].y);
					nurbsCurve->PtsDerivs(tmpU + 0.01, 1, Der_2);
					tmpVec2 = point2d(Der_2[0].x - Der[0].x, Der_2[0].y - Der[0].y);

					tmpVec1 = tmpVec1.Normalize();
					tmpVec2 = tmpVec2.Normalize();

					if (tmpVec1.y < 0.0 && tmpVec2.y < 0.0) {
						//�õ��������0.05������������������Ϊ1��2����ʱ�����õ���Ϊ����
						count++;
					}
				}

				else {
					//�ô���������ˮƽ
					count++;
				}
			}
		}
	}

	//if (numintpt != 0) {
	//	//���ڽ������
	//	for (int i = 0; i < numintpt; i++) {
	//		double tmpU = intpar[i];	//���㴦�Ĳ���ֵ
	//		int flagU = CalU(tmpU);
	//		//������м�Ĳ���
	//		varray<Vec4> Der;//CNurbs�����ʽ���Ľ��������������
	//		if (0 == flagU) {
	//			//����Ϊ���߲���Ϊ0���ĵ�
	//			std::vector<std::vector<double>> mderive;
	//			int dim = CalLinePoint(curve, 0, mderive);
	//			if (mderive[0][0] < strLine[0][0]) {
	//				//��������ߵ��ӳ�����
	//				continue;
	//			}
	//			if (nurbsCurve) {
	//				//CNurbs�����ʽ�����߲�Ϊ��
	//				nurbsCurve->PtsDerivs(0, 1, Der);
	//			}
	//			if (Der.size() < 2) assert(0);
	//			if (dim == 2) {
	//				//ά����2
	//				if (fabs(mderive[1][1] - 0.0) < 1e-5|| fabs(Der[1].y-0.0) < 1e-5) {
	//					//�ô�������ˮƽ
	//					std::vector<std::vector<double>> mderive_2;
	//					CalLinePoint(curve, 0.05, mderive_2);
	//					point2d tmpVec = point2d(mderive_2[0][0] - mderive[0][0], mderive_2[0][1] - mderive[0][1]);
	//					tmpVec = tmpVec.Normalize();
	//					if (tmpVec.y < 0.0) {
	//						//��������3��4���ޱ�ʾ���õ���Ϊ����
	//						count++;
	//					}
	//				}
	//				else {
	//					//�ô���������ˮƽ
	//					if (mderive[1][1] < 0.0) {
	//						//�õ�����������3/4���ޣ�������Ϊ����
	//						count++;
	//					}
	//				}
	//			}
	//		}
	//		else if (1 == flagU) {
	//			//����Ϊ���߲���Ϊ1���ĵ�
	//			std::vector<std::vector<double>> mderive;
	//			int dim = CalLinePoint(curve, 1, mderive);
	//			if (mderive[0][0] < strLine[0][0]) {
	//				//��������ߵ��ӳ�����
	//				continue;
	//			}
	//			if (nurbsCurve) {
	//				//CNurbs�����ʽ�����߲�Ϊ��
	//				nurbsCurve->PtsDerivs(1, 1, Der);
	//			}
	//			if (Der.size() < 2) assert(0);
	//			if (dim == 2) {
	//				//ά����2
	//				if (fabs(mderive[1][1] - 0.0) < 1e-5 || fabs(Der[1].y - 0.0) < 1e-5) {
	//					//�ô�������ˮƽ
	//					std::vector<std::vector<double>> mderive_2;
	//					CalLinePoint(curve, 0.95, mderive_2);
	//					point2d tmpVec = point2d(mderive_2[0][0] - mderive[0][0], mderive_2[0][1] - mderive[0][1]);
	//					tmpVec = tmpVec.Normalize();
	//					if (tmpVec.y < 0.0) {
	//						//��������3��4���ޱ�ʾ���õ���Ϊ����
	//						count++;
	//					}
	//				}
	//				else {
	//					//�ô���������ˮƽ
	//					if (mderive[1][1] < 0.0) {
	//						//�õ�����������3/4���ޣ�������Ϊ����
	//						count++;
	//					}
	//				}
	//			}
	//		}

	//		else {
	//			//����Ϊ������һ��
	//			std::vector<std::vector<double>> mderive;
	//			Der.clear();
	//			int dim = CalLinePoint(curve, tmpU, mderive);
	//			if (mderive[0][0] < strLine[0][0]) {
	//				//��������ߵ��ӳ�����
	//				continue;
	//			}
	//			if (nurbsCurve) {
	//				//CNurbs�����ʽ�����߲�Ϊ��
	//				nurbsCurve->PtsDerivs(tmpU, 1, Der);
	//			}
	//			if (Der.size() < 2) assert(0);
	//			if (dim == 2) {
	//				//ά����2
	//				if (fabs(mderive[1][1] - 0.0) < 1e-5 || fabs(Der[1].y - 0.0) < 1e-5) {
	//					//�ô�������ˮƽ
	//					assert((tmpU - 0.05) > 0.0);
	//					assert((tmpU + 0.05) < 1.0);

	//					point2d tmpVec1, tmpVec2;

	//					if (fabs(mderive[1][1] - 0.0) < 1e-5) {
	//						//��sisl��õĽ����������
	//						std::vector<std::vector<double>> mderive_2, mderive_3;//�õ����0.05�����������꼰������

	//						CalLinePoint(curve, tmpU - 0.05, mderive_2);
	//						tmpVec1 = point2d(mderive_2[0][0] - mderive[0][0], mderive_2[0][1] - mderive[0][1]);
	//						
	//						CalLinePoint(curve, tmpU + 0.05, mderive_3);
	//						tmpVec2 = point2d(mderive_3[0][0] - mderive[0][0], mderive_3[0][1] - mderive[0][1]);	
	//					}
	//					else if(fabs(Der[1].y - 0.0) < 1e-5){
	//						//CNurbs��õĽ����������
	//						varray<Vec4> Der_1, Der_2;
	//						nurbsCurve->PtsDerivs(tmpU-0.05, 1, Der_1);
	//						tmpVec1 = point2d(Der_1[0].x - Der[0].x, Der_1[0].y - Der[0].y);
	//						nurbsCurve->PtsDerivs(tmpU+0.05, 1, Der_2);
	//						tmpVec2 = point2d(Der_2[0].x - Der[0].x, Der_2[0].y - Der[0].y);
	//					}

	//					tmpVec1 = tmpVec1.Normalize();
	//					tmpVec2 = tmpVec2.Normalize();
	//					
	//					if (tmpVec1.y < 0.0 && tmpVec2.y < 0.0) {
	//						//�õ��������0.05������������������Ϊ1��2����ʱ�����õ���Ϊ����
	//						count++;
	//					}						
	//				}
	//				else {
	//					//�ô���������ˮƽ
	//					count++;
	//				}
	//			}
	//		}
	//	}
	//}

	if (intpar != nullptr)
		delete[] intpar;

	return count;
}

int PointRelatePolygon(Vec3 p,  varray<Spline>& pol)
{
	int count = 0;															//���߷������εĽ�������
	int flag = -1;															//��������ż�ж���־��0��ʾż�����ⲿ��1��ʾ�������ڲ�
	double xMin, xMax;														//�����x������Сֵ�����ֵ
	double vecLength;														//���߷��߶γ���
	Spline ray;															//�����߶�
	SISLCurve * sislRay;													//�����߶�SISL��ʾ
	varray<SISLCurve *> sislLines;											//SISL��ʽ���߱�ʾ
	if (pol.size() < 4)
	{
		cout << "PointRelatePolygon��Size<4" << endl;
		//return -1;
	}
	if (pol.size() < 3)
	{
		cout << "PointRelatePolygon��Size<3" << endl;
		return -1;
	}

	//������ν�����ʱ������,��������ж˵�����
	varray<Vec3> pointsCord = OrderLinesAntioclock(pol);

	//���������ߵ�x�������
	xMin = pol[0].m_CtrlPts[0].x;
	xMax = pol[0].m_CtrlPts[0].x;
	for (auto i : pol)
	{
		xMin = min(xMin, i.m_CtrlPts[0].x);
		xMin = min(xMin, (i.m_CtrlPts.end() - 1)->x);
		xMax = max(xMax, i.m_CtrlPts[0].x);
		xMax = max(xMax, (i.m_CtrlPts.end() - 1)->x);
	}
	vecLength = 3 * (xMax - xMin);									//���߶γ�����Ϊ�����x���������������Ա�֤�н���

	//�ж��õ��Ƿ�Ϊ����εĶ˵�,���ǣ�����0
	for (auto i : pointsCord)
	{
		if (JudgeTwoPointsCoincide(p, i))
		{
			cout << "PointRelatePolygon���õ�Ϊ����ζ˵�" << endl;
			return 0;
		}
	}

	//����������ת��ΪSISL��ʽ
	for (int i = 0; i < pol.size(); i++)
	{
		SISLCurve * tmpCurve = nullptr;
		tmpCurve = NurbsLineToSislLine(pol[i],2);
		sislLines.push_back(tmpCurve);
		//freeCurve(tmpCurve);										//����ɾ������ʱSISL����
		tmpCurve = nullptr;
	}

	//�ж��õ��Ƿ���������(�Ƕ˵㴦)
	for (int i = 0; i < sislLines.size(); i++)
	{
		double* pt = Point3dToSislPoint(p);
		double u = PointIntersectNurbsLine(pt, sislLines[i], 2);
		int endP = CalU(u);
		if (abs(u + 1)>DBL_EPSILON && endP == -1)
		{
			cout << "PointRelatePolygon���õ����������Ҳ�Ϊ�˵�" << endl;
			return 3;
		}
	}

	////�������㴴��NURBS����(�˴�Ϊ���߷��߶�)
	//ray.m_Degree = 1;
	//ray.m_Knots.push_back(0);
	//ray.m_Knots.push_back(0);
	//ray.m_Knots.push_back(1);
	//ray.m_Knots.push_back(1);
	//ray.m_CtrlPts.push_back(Vec4(p, 1));
	//ray.m_CtrlPts.push_back(Vec4(p.x + vecLength, p.y, p.z, 1));
	//sislRay = NurbsLineToSislLine(ray);

	//����ֱ���߶Σ�������
	varray<Vec3> strLine;
	strLine.push_back(p);
	strLine.push_back(Vec3(p.x + vecLength, p.y, p.z));

	//����sislRay����������������
	for (int i = 0; i < sislLines.size(); i++)
	{
		count += LineIntersectNurbs(sislLines[i], strLine, &pol[i]);
	}

	//ɾ��sisl����
	for (int i = 0; i < sislLines.size(); i++)
	{
		freeCurve(sislLines[i]);
	}

	flag = count % 2;
	if (0 == flag) {
		return 2;
	}
	else if (1 == flag) {
		return 1;
	}
	else {
		assert(0);
		return -1;//��������
	}
}

bool PolygonsRelationship(varray<Spline> pola, varray<Spline> polb)
{
	varray<Vec3>p = OrderLinesAntioclock(pola);
	for (auto &i : pola) {
		//����u=0.5�����겢����
		Vec3 q = i.GetLinePoint(0.5);
		p.push_back(q);
	}
	for (auto &i : p) {
		int res = PointRelatePolygon(i, polb);
		if (res == -1) assert(0);
		if (res == 2) return false;
	}
	return true;
}

SfCtainTreeNode::SfCtainTreeNode() : isSort(false), quaded(false)
{
}

SfCtainTreeNode::SfCtainTreeNode(const varray<Spline>& outlines) : isSort(false), quaded(false)
{
	this->outLines = outlines;
}

SfCtainTreeNode::SfCtainTreeNode(const varray<Spline>& outlines, const varray<Spline>& conlines) :  isSort(false), quaded(false)
{
	this->outLines = outlines;
	this->conLines = conlines;
}

SfCtainTreeNode::SfCtainTreeNode(const varray<Spline>& outlines, const varray<Spline>& conlines, varray<int> seg): quaded(false)
{
	this->outLines = outlines;
	this->conLines = conlines;
	this->isSeg = seg;
}

SfCtainTreeNode::SfCtainTreeNode(const varray<Spline>& outlines, const varray<Spline>& conlines, varray<int> seg, bool genus): quaded(false)
{
	this->outLines = outlines;
	this->conLines = conlines;
	this->isSeg = seg;
	this->isGenus = genus;
}

SfCtainTreeNode::~SfCtainTreeNode()
{
	//�ͷ������ڵĶ�����
	list<SfCtainTreeNode*>::iterator it = this->childs.begin();
	for (; it != this->childs.end(); it++) {
		delete (*it);
		(*it) = nullptr;
	}
}

void SfCtainTreeNode::GetOutPoints(varray<Vec3>& ps)
{
	ps.clear();
	ps = this->vertex;
}

void SfCtainTreeNode::GetOutLines(varray<Spline>& outline)
{
	outline.clear();
	outline.resize(this->outNumber.size());
	for (int i = 0; i < this->outNumber.size(); ++i) {
		outline[i] = this->allLines[outNumber[i]];
	}
}

void SfCtainTreeNode::GetSurfs(varray<SplineSurface>& surfs)
{
	surfs.clear();
	if (this->isGenus) {
		//���ýڵ�Ϊ��������
		return;
	}
	for (auto i : this->quadPolNumber) {
		varray<Spline> sfLines;
		SplineSurface tmpSurf;
		QuadOrder(i, this->allLines, sfLines);
		assert(sfLines.size() == 4);//���� assert������:�����ֵΪ�٣���Ϊ0������ abort ����ֹ��������

		tmpSurf.CoonsInterpolate(sfLines);
		surfs.push_back(tmpSurf);
	}
}


int MyPolygon::GetFrontNumber(int pointNumber)
{
	int pLength = this->p_Points.size();
	if (pointNumber == 0)
		return pLength - 1;
	else
		return pointNumber - 1;
}
