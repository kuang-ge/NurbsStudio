#include "stdafx.h"
#include "PSO.h"
#include <cstdlib>

//���캯��
PSO::PSO()
{
	m_Dim = 1;
	m_c1 = 0.2;
	m_c2 = 0.8;
	m_r = 0.8;
	m_PNum = 16;
	m_iteraNum = 100;
	m_error = 0.00001;

	m_GloItNum = 0;
}

/*��������
double c1: ������֪ϵ��
double c2: ���֪ʶϵ��
double r: �ٶȸ���Լ������
int PNum: ��������
int iteraNum: ����������
varray<double> error: �޶����*/
void PSO::SetPSOPara(const int PNum, const int iteraNum, const double c1, const double c2, const double r,
	const double error)
{
	m_c1 = c1;
	m_c2 = c2;
	m_r = r;
	m_PNum = PNum;
	m_iteraNum = iteraNum;
	m_error = error;
}

/*���÷�Χ
Pmin��������Сλ��
Pmax���������λ��
Vmin�������ٶ���Сλ��
Vmax�������ٶ����λ��
*/
void PSO::SetRange(const varray<double>& Pmin, const varray<double>& Pmax, const varray<double>& Vmin, const varray<double>& Vmax)
{
	m_Pmin = Pmin;
	m_Pmax = Pmax;
	m_Vmin = Vmin;
	m_Vmax = Vmax;
}

//error: �޶����
void PSO::SetError(const double error)
{
	m_error = error;
}

//��ȡȫ���������Ӽ���Ӧ����Ӧ��
void PSO::GetResult(varray<double> & PGlobal, double & DisGol)
{
	PGlobal = m_PGlobal;
	DisGol = m_DisGol;
}

//������������Ⱥ
void PSO::PSOUpdate()
{
	m_Dim = m_Pmin.size();

	InitialPV();
	InitialDis();

	//ѭ��1����������ѭ������
	for (int t = 0; t < m_iteraNum; ++t)
	{
		//�ٶȸ���Ȩ��
		double w = 0.55*(1 - 1.0*t / m_iteraNum) + 0.4;

		//ѭ��2����ÿ�������ٶȺ�λ�ý��и���
		for (int i = 0; i < m_PNum; ++i)
		{
			//�����ٶ�
			varray<double> tV;
			for (int j = 0; j < m_Dim; ++j)
			{
				double epsilon = (double)std::rand() / RAND_MAX;
				double yita = (double)std::rand() / RAND_MAX;

				double V = w*m_velocity[i][j] + m_c1*epsilon*(m_PLocal[i][j] - m_P[i][j]) + m_c2*yita*(m_PGlobal[j] - m_P[i][j]);

				if (V < m_Vmin[j])
					V = m_Vmin[j];
				else if (V > m_Vmax[j])
					V = m_Vmax[j];

				tV.push_back(V);
			}
			m_velocity[i] = tV;

			//����λ��
			varray<double> tP;
			for (int j = 0; j < m_Dim; ++j)
			{
				double P = m_r*m_velocity[i][j] + m_P[i][j];

				if (P < m_Pmin[j])
					P = m_Pmin[j];
				else if (P > m_Pmax[j])
					P = m_Pmax[j];

				tP.push_back(P);
			}
			m_P[i] = tP;
		}//!ѭ��2

		 //��������ֵ
		if (UpdateLG())
			break;
	}//!ѭ��1
}

//��ʼ��λ�ü��ٶ�,����ɷ�Χ�޶�
void PSO::InitialPV()
{
	m_P.clear();
	m_velocity.clear();
	m_PLocal.clear();

	for (int i = 0; i < m_PNum; ++i)
	{
		varray<double> P, V;
		for (int j = 0; j < m_Dim; ++j)
		{
			//λ��
			double tP = ((double)std::rand() / RAND_MAX) * (m_Pmax[j] - m_Pmin[j]) + m_Pmin[j];
			P.push_back(tP);

			//�ٶ�
			double tV = ((double)std::rand() / RAND_MAX) * (m_Vmax[j] - m_Vmin[j]) + m_Vmin[j];
			V.push_back(tV);
		}
		m_P.push_back(P);
		m_velocity.push_back(V);
	}
}

//��ʼ������ֵ��ȫ������ֵ
inline void PSO::InitialDis()
{
	m_PLocal = m_P;
	CalLocDis(m_PLocal, m_DisLoc);
	m_PGlobal = m_P[0];
	m_DisGol = m_DisLoc[0];

	//����ȫ������ֵ
	for (int i = 0; i < m_DisLoc.size(); ++i)
	{
		if (m_DisLoc[i] < m_DisGol)
		{
			m_DisGol = m_DisLoc[i];
			m_PGlobal = m_PLocal[i];
		}
	}
}

//��������ֵ,�ж��Ƿ�ﵽ����
bool PSO::UpdateLG()
{
	bool flag = false;
	//��һ��ȫ������
	double sg = m_DisGol;
	//��ǰ������Ӧ��
	varray<double> dis;
	CalLocDis(m_P, dis);
	//��������ֵ
	bool isNoChange = true;
	for (int i = 0; i < dis.size() && isNoChange; ++i)
	{
		if (i != 0)
		{
			double dsi = abs(long(dis[i] - dis[i - 1]));
			if (dsi == 0)
				continue;
			isNoChange = false;
		}
		if (dis[i] < m_DisLoc[i])
		{
			m_DisLoc[i] = dis[i];
			m_PLocal[i] = m_P[i];
			if (dis[i] < m_DisGol)
			{
				m_DisGol = dis[i];
				m_PGlobal = m_P[i];
			}
		}
	}
	double ds = abs(long(sg - m_DisGol));
	if (ds != 0)
		++m_GloItNum;
	if ((ds < m_error && m_GloItNum > m_iteraNum / 4) || isNoChange)
		flag = true;
	return flag;
}

//�������и�����Ӧ��
void PSO::CalLocDis(const varray<varray<double>>& P, varray<double>& dis)
{
	dis.clear();
	for (int i = 0; i < P.size(); ++i)
	{
		double s = CalDis(P[i]);
		dis.push_back(s);
	}
}


