//��ά��PSO�㷨
//@ USST ���ı� 2018
#ifndef CPSO
#define CPSO
#include "varray.h"

class PSO
{
protected:
	int m_Dim;//ά��
	int m_PNum;//��������
	int m_iteraNum;//����������
	double m_c1;//������֪ϵ��
	double m_c2;//���֪ʶϵ��
	double m_r;//�ٶȸ���Լ������
	double m_error;//�޶����

	varray<varray<double>> m_P;//����Ⱥλ��
	varray<varray<double>> m_velocity;//����Ⱥ�ٶ�

	varray<double> m_PGlobal;//ȫ������ֵ����λ��
	double m_DisGol;//ȫ��������Ӧ��
	
	varray<double> m_Pmin;//������Сλ��
	varray<double> m_Pmax;//�������λ��
	varray<double> m_Vmin;//�����ٶ���Сλ��
	varray<double> m_Vmax;//�����ٶ����λ��
	
private:
	int m_GloItNum;
	varray<double> m_DisLoc;//����������Ӧ��
	varray<varray<double>> m_PLocal;//ÿ�����Ӹ�������ֵ

public:
	//���캯��
	PSO();

	//��������
	virtual ~PSO() {}

	/*��������
	c1: ������֪ϵ��
	c2: ���֪ʶϵ��
	r: �ٶȸ���Լ������
	PNum: ��������
	iteraNum: ����������
	error: �޶����*/
	void SetPSOPara(const int PNum, const int iteraNum, const double c1, const double c2, const double r,
		const double error);

	/*���÷�Χ
	Pmin��������Сλ��
	Pmax���������λ��
	Vmin�������ٶ���Сλ��
	Vmax�������ٶ����λ��
	*/
	void SetRange(const varray<double>& Pmin, const varray<double>& Pmax, const varray<double>& Vmin, const varray<double>& Vmax);

    //error: �޶����
	void SetError(const double error);

	//��ȡȫ���������Ӽ���Ӧ����Ӧ��
	void GetResult(varray<double>& PGlobal, double& DisGol);

	//������������Ⱥ
	void PSOUpdate();

protected:
	//��Ӧ�Ⱥ���
	//����������д
	virtual double CalDis(const varray<double>& P) = 0;

	//��ʼ��λ�ü��ٶ�,����ɷ�Χ�޶�
	virtual void InitialPV();

private:
	//��ʼ������ֵ��ȫ������ֵ
	void InitialDis();

	//��������ֵ,�ж��Ƿ�ﵽ����
	bool UpdateLG();

	//�������и�����Ӧ��
	void CalLocDis(const varray<varray<double>>& P, varray<double>& dis);
};

#endif // !CPSO
