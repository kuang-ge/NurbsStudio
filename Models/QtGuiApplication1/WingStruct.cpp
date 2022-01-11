#include "WingStruct.h"

using YN::NurbsVol;
using YN::NurbsSurface;
using YN::NurbsLine;
using YN::Cube;
using YN::Surface;

YN::WingStruct::WingStruct()
{
	_root = NULL;
}

YN::WingStruct::WingStruct(string path, bool tag1 = true)
{
	vector<Cube *> m_HexVolumes=ReadVolumeTxt(path,tag1);
	abstractInfo(m_HexVolumes);
}

YN::WingStruct::WingStruct(vector<NurbsVol> nurbsVols)
{
	vector<Cube *> m_HexVolumes;//������ļ���
	m_HexVolumes = transNurbsVoltoCube(nurbsVols);//���ת��
	abstractInfo(m_HexVolumes);
}

YN::WingStruct::~WingStruct()
{
	vector<Cube *> a= bordFirstSearch();
	for (int i = 0; i != a.size(); ++i)
	{
		//���������Ƭ����ȫ��ɾ������ֹ�ڴ�й¶
		delete a[i]->_u0v0;
		a[i]->_u0v0 = NULL;
		delete a[i]->_u1v1;
		a[i]->_u1v1 = NULL;
		delete a[i]->_v0w0;
		a[i]->_v0w0 = NULL;
		delete a[i]->_v1w1;
		a[i]->_v1w1 = NULL;
		delete a[i]->_u0w0;
		a[i]->_u0w0 = NULL;
		delete a[i]->_u1w1;
		a[i]->_u1w1 = NULL;
		delete a[i];
		a[i] = NULL;
	}
}

vector<Cube *> YN::WingStruct::ReadVolumeTxt(string path,bool tag=true)
{
	vector<Cube *> m_HexVolumes;
	ifstream ifs;
	ifs.open(path);
	assert(ifs.is_open());

	int i;

	while (!ifs.eof())
	{
		string str;
		ifs >> str;
	contin:
		if (str == "PN")
		{
			int count;
			ifs >> count;
			for (int i = 0; i != count; ++i)
			{
				Cube *t = new Cube;
				m_HexVolumes.push_back(t);
			}
		}
		if (str == "PI")
		{
			ifs >> i;
		}
		if (str == "OD")
		{
			ifs >> m_HexVolumes.at(i)->_Vol_Data._u_Degree;
			ifs >> m_HexVolumes.at(i)->_Vol_Data._v_Degree;
			ifs >> m_HexVolumes.at(i)->_Vol_Data._w_Degree;
			//�����������Ĵ�������Ϊ��ȡ�ļ���ʱ����㲻��ȷ������ֶ�������ȷ�ˣ����԰������ע��
			m_HexVolumes.at(i)->_Vol_Data._u_Degree--;
			m_HexVolumes.at(i)->_Vol_Data._v_Degree--;
			m_HexVolumes.at(i)->_Vol_Data._w_Degree--;
		}
		if (str == "UK")
		{
			int count;
			ifs >> count;
			m_HexVolumes.at(i)->_Vol_Data._u_Knots->resize(count);
			for (int j = 0; j < count; j++)
			{
				double temp;
				ifs >> temp;
				(*m_HexVolumes.at(i)->_Vol_Data._u_Knots)[j] = temp;
			}
		}
		if (str == "VK")
		{
			int count;
			ifs >> count;
			m_HexVolumes.at(i)->_Vol_Data._v_Knots->resize(count);
			for (int j = 0; j < count; j++)
			{
				double temp;
				ifs >> temp;
				(*m_HexVolumes.at(i)->_Vol_Data._v_Knots)[j] = temp;
			}
		}
		if (str == "WK")
		{
			int count;
			ifs >> count;
			m_HexVolumes.at(i)->_Vol_Data._w_Knots->resize(count);
			for (int j = 0; j < count; j++)
			{
				double temp;
				ifs >> temp;
				(*m_HexVolumes.at(i)->_Vol_Data._w_Knots)[j] = temp;
			}
		}
		if (str == "CP")
		{
			ifs >> m_HexVolumes.at(i)->_Vol_Data._u_Num;
			ifs >> m_HexVolumes.at(i)->_Vol_Data._v_Num;
			ifs >> m_HexVolumes.at(i)->_Vol_Data._w_Num;
			ifs >> str;
			while (str != "PI")
			{
				point4d v;
				v.x = stof(str);
				ifs >> v.y;
				ifs >> v.z;
				if (tag == true)
				{
					ifs >> v.w;
				}
				m_HexVolumes.at(i)->_Vol_Data._ControlPts->push_back(v);
				ifs >> str;
				if (ifs.eof())
				{
					return m_HexVolumes;
				}
			}
			goto contin;
		}
	}
	
	//��ȡ����֮���������Ƿ���ȷ
	for (auto iter = m_HexVolumes.begin(); iter != m_HexVolumes.end(); ++iter)
	{
		//�����Ƶ����
		int num = (*iter)->_Vol_Data._u_Num*(*iter)->_Vol_Data._v_Num*(*iter)->_Vol_Data._w_Num;
		if (num != (*iter)->_Vol_Data._ControlPts->size())
		{
			m_HexVolumes.clear();
		}
		//��� �ڵ�ʸ��=���Ƶ����+����+1
		if((*iter)->_Vol_Data._u_Knots->size()!= (*iter)->_Vol_Data._u_Num+ 
			(*iter)->_Vol_Data._u_Degree+1)
			m_HexVolumes.clear();
		if ((*iter)->_Vol_Data._v_Knots->size() != (*iter)->_Vol_Data._v_Num +
			(*iter)->_Vol_Data._v_Degree + 1)
			m_HexVolumes.clear();
		if ((*iter)->_Vol_Data._w_Knots->size() != (*iter)->_Vol_Data._w_Num +
			(*iter)->_Vol_Data._w_Degree + 1)
			m_HexVolumes.clear();
	}

	return m_HexVolumes;
}

vector<point4d>  YN::WingStruct::OutputParaVolumeDataTxt(string filename, string Cp)
{
	vector<Cube *> m_HexVolumes;
	m_HexVolumes = bordFirstSearch();
	//���ģ����Ϣ
	ofstream ofs(filename);
	ofs << "PN " << " " << m_HexVolumes.size() << "\n";  //patch number
	for (int i = 0; i < m_HexVolumes.size(); i++)
	{
		ofs << "PI" << " " << i << "\n";   //��ǰƬ��id�š�

		ofs << "OD" << " " << m_HexVolumes.at(i)->_Vol_Data._u_Degree + 1 << " " << m_HexVolumes.at(i)->_Vol_Data._v_Degree + 1 << " " << m_HexVolumes.at(i)->_Vol_Data._w_Degree + 1 << "\n"; //Order

		ofs << "UK" << " " << m_HexVolumes.at(i)->_Vol_Data._u_Knots->size() << "\n";
		for (int i = 0; i < m_HexVolumes.at(i)->_Vol_Data._u_Knots->size(); i++)
			ofs << m_HexVolumes.at(i)->_Vol_Data._u_Knots->at(i) << " ";   //�ڵ�����
		ofs << "\n";
		ofs << "VK" << " " << m_HexVolumes.at(i)->_Vol_Data._v_Knots->size() << "\n";
		for (int i = 0; i < m_HexVolumes.at(i)->_Vol_Data._v_Knots->size(); i++)
			ofs << m_HexVolumes.at(i)->_Vol_Data._v_Knots->at(i) << " ";   //�ڵ�����
		ofs << "\n";
		ofs << "WK" << " " << m_HexVolumes.at(i)->_Vol_Data._w_Knots->size() << "\n";
		for (int i = 0; i < m_HexVolumes.at(i)->_Vol_Data._w_Knots->size(); i++)
			ofs << m_HexVolumes.at(i)->_Vol_Data._w_Knots->at(i) << " ";   //�ڵ�����
		ofs << "\n";

		ofs << "CP" << " " << m_HexVolumes.at(i)->_Vol_Data._u_Num << " " << m_HexVolumes.at(i)->_Vol_Data._v_Num << " " << m_HexVolumes.at(i)->_Vol_Data._w_Num << "\n"; //control point number
		for (int i = 0; i < m_HexVolumes.at(i)->_Vol_Data._ControlPts->size(); i++)
			ofs << m_HexVolumes.at(i)->_Vol_Data._ControlPts->at(i).x << " " 
			<< m_HexVolumes.at(i)->_Vol_Data._ControlPts->at(i).y << " " 
			<< m_HexVolumes.at(i)->_Vol_Data._ControlPts->at(i).z << " " 
			<< m_HexVolumes.at(i)->_Vol_Data._ControlPts->at(i).w << "\n";
	}
	ofs.close();
	//������Ƶ�ȫ�ֱ��
	ofstream idofs(Cp);
	idofs << "PN " << " " << m_HexVolumes.size() << "\n";  //patch number
	for (int i = 0; i < m_HexVolumes.size(); i++)
	{
		int idxnum = m_HexVolumes.at(i)->_volGloblControlPointID.size();
		idofs << "PN" << " " << i << "\n";
		idofs << "PI" << " " << idxnum << "\n";
		for (int j = 0; j < idxnum; j++)
		{
			idofs << m_HexVolumes.at(i)->_volGloblControlPointID[j] << " ";
		}
		idofs << "\n";
	}

	//���ʩ��Լ�������Ŀ��Ƶ�
	vector<point4d> Poin;
	point4d po;
	//��������֮��ľ���
	auto function = [](point4d p, float A, float B, float C, float D) {
		float temp = sqrt(A*A + B * B + C * C);
		return 1.0*abs(p.x*A + p.y*B + p.z*C + D) / temp;
	};
	//����ƽ��������������ƽ������Ŀ��Ƶ�
	auto Function = [&](float A, float B, float C, float D)
	{
		vector <int> result;
		float distance;
		for (int i = 0; i != m_HexVolumes.size(); ++i)
		{
			for (int j = 0; j != m_HexVolumes[i]->_Vol_Data._ControlPts->size(); ++j)
			{
				if (result.size() == 0)
				{
					distance = function(m_HexVolumes[i]->_Vol_Data._ControlPts->at(j) , A, B, C, D);
					result.push_back(m_HexVolumes[i]->_volGloblControlPointID[j]);
					po = m_HexVolumes[i]->_Vol_Data._ControlPts->at(j);
					Poin.push_back(po);
					continue;
				}
				else
				{
					float temp = function(m_HexVolumes[i]->_Vol_Data._ControlPts->at(j), A, B, C, D);
					if (abs(temp - distance) <= 0.1)
					{
						result.push_back(m_HexVolumes[i]->_volGloblControlPointID[j]);
						po = m_HexVolumes[i]->_Vol_Data._ControlPts->at(j);
						Poin.push_back(po);
					}
					else
					{
						if (distance - temp > 0)
						{
							result.clear();
							Poin.clear();
							distance = function(m_HexVolumes[i]->_Vol_Data._ControlPts->at(j), A, B, C, D);
							result.push_back(m_HexVolumes[i]->_volGloblControlPointID[j]);
							po = m_HexVolumes[i]->_Vol_Data._ControlPts->at(j);
							Poin.push_back(po);
						}
						else
						{
							continue;
						}
					}
				}
			}
		}
		return result;
	};
	//ȥ���ظ����Ƶ�
	auto single = [](vector<int> &a)
	{
		for (auto iter = a.begin(); iter != a.end(); ++iter)
		{
			for (auto iter2 = iter + 1; iter2 != a.end(); )
			{
				if (*iter2 == *iter)
				{
					iter2 = a.erase(iter2);
				}
				else
				{
					iter2++;
				}
			}
		}
	};

	if (m_HexVolumes.size() != 0)
	{
		vector<int> a = Function(0, 0, 1, -50);
		single(a);
		vector<point4d> aaaa = Poin;
		idofs << "WC" << " ";
		for (int i = 0; i != a.size(); ++i)
		{
			idofs << a[i] << " ";
		}
		idofs << "\n";
		a.clear();
		a = Function(0, 0, 1, 50);
		Poin.insert(Poin.end(), aaaa.begin(), aaaa.end());
		single(a);
		idofs << "WF" << " ";
		for (int i = 0; i != a.size(); ++i)
		{
			idofs << a[i] << " ";
		}
		idofs << "\n";
	}

	return Poin;
}

void  YN::WingStruct::OutputParaVolumeDataVTK(string filename)
{
	vector<Cube *> m_HexVolumes;
	m_HexVolumes = bordFirstSearch();
	ofstream file; //�������ʽ���ļ�
	file.open(filename);
	
	int uSegnum, vSegnum, wSegnum, patchPtNumCount, patchNumCount;
	uSegnum = vSegnum = wSegnum = 10;
	patchPtNumCount = (uSegnum + 1)*(vSegnum + 1)*(wSegnum + 1);  //�����Ŀ
	int patchNum = m_HexVolumes.size();
	point4d pt;
	int started = 0;

	file << "# vtk DataFile Version 2.0" << "\n" << "TET" << "\n" << "ASCII" << "\n";
	file << "\n" << "DATASET UNSTRUCTURED_GRID" << "\n";
	file << "POINTS " << patchNum * patchPtNumCount << " float" << "\n";

	for (int patID = 0; patID < patchNum; patID++)
	{
		for (int k = 0; k <= wSegnum; k++)
		{
			for (int j = 0; j <= vSegnum; j++)
			{
				for (int i = 0; i <= uSegnum; i++)
				{
					pt = m_HexVolumes[patID]->_Vol_Data.GetVolPoint(i*1.f / uSegnum, j*1.f / vSegnum, k*1.f / wSegnum);
					file << pt.x << " " << pt.y << " " << pt.z << " " << "\n";
				}
			}
		}
	}

	patchNumCount = uSegnum * vSegnum*wSegnum;   //��Ԫ��Ŀ��
	file << "CELLS " << patchNum * patchNumCount << " " << 9 * patchNum * patchNumCount << "\n";
	for (int patID = 0; patID < patchNum; patID++)
	{
		started = patID * patchPtNumCount;
		for (int k = 0; k < wSegnum; k++)
		{
			for (int j = 0; j < vSegnum; j++)
			{
				for (int i = 0; i < uSegnum; i++)
				{
					file << "8 " << started + k * (uSegnum + 1)*(vSegnum + 1) + (uSegnum + 1)*j + i << " "
						<< started + k * (uSegnum + 1)*(vSegnum + 1) + (uSegnum + 1)*j + i + 1 << " "
						<< started + k * (uSegnum + 1)*(vSegnum + 1) + (uSegnum + 1)*(j + 1) + i + 1 << " "
						<< started + k * (uSegnum + 1)*(vSegnum + 1) + (uSegnum + 1)*(j + 1) + i << " "
						<< started + (k + 1)*(uSegnum + 1)*(vSegnum + 1) + (uSegnum + 1)*j + i << " "
						<< started + (k + 1)*(uSegnum + 1)*(vSegnum + 1) + (uSegnum + 1)*j + i + 1 << " "
						<< started + (k + 1)*(uSegnum + 1)*(vSegnum + 1) + (uSegnum + 1)*(j + 1) + i + 1 << " "
						<< started + (k + 1)*(uSegnum + 1)*(vSegnum + 1) + (uSegnum + 1)*(j + 1) + i << "\n";
				}
			}
		}
	}

	file << "CELL_TYPES " << patchNum * patchNumCount << "\n";
	for (int l = 0; l < patchNum * patchNumCount; l++)
	{
		file << "12" << "\n";
	}
	file.close();
}

void YN::WingStruct::abstractInfo(vector<Cube *>& m_HexVolumes)
{
	vector<Surface *> m_all4sidePatch;//������Ƭ�ļ���
	//��ȡ������Ϣ
	m_all4sidePatch = upDateSurface(m_HexVolumes);//��ĸ���

	//��ȡ������Ϣ
	//�˷�ʱ��!�ٶ�̫��
	setTwinPatch(m_HexVolumes, m_all4sidePatch);

	//��ȡ��������Ϣ,Ҫ�����е�ƬUVW����Ŀ��Ƶ����ȫ����ͬ�ſ��Խ���
	Order(m_HexVolumes);

	//���ø��ڵ�
	_root = m_HexVolumes[0];

	//����ȫ�ֿ��Ƶ����
	setGloblControlPointId();
}

bool YN::WingStruct::setTwinPatch(vector<Cube *>& m_HexVolumes, vector<Surface *>& m_all4sidePatch)
{
	int patchnum = m_all4sidePatch.size();
	if (patchnum < 2)
		return false;
	for (int i = 0; i < patchnum; i++)
	{
		m_all4sidePatch.at(i)->_haveNextSuf = false;
		m_all4sidePatch.at(i)->_nextSurface = nullptr;
		m_all4sidePatch.at(i)->_pointOther = nullptr;
	}
	vector<int> idxarr;
	for (int i = 0; i < patchnum - 1; i++)
	{
		idxarr.clear();
		for (int j = i + 1; j < patchnum; j++)
		{
			if (m_all4sidePatch.at(j)->_haveNextSuf != false)
				continue;
			if (m_all4sidePatch.at(i)->_surfaceData.isTwoSurfaceSame(m_all4sidePatch.at(j)->_surfaceData))
			{
				idxarr.push_back(j);
				break;
			}
		}
		if (idxarr.size() < 1)
			continue;
		int minidx = idxarr.at(0);
		if (idxarr.size() > 1)  //�ҵ����治ֹһ����Ŀǰ��������ġ�����ҵ�����Ƭ��ֹһ����������Ծ�����С���档
		{
			float minlen = m_all4sidePatch.at(i)->_surfaceData.GetDistanceBetweenTwoSurface(m_all4sidePatch.at(idxarr.at(0))->_surfaceData);
			for (int k = 1; k < idxarr.size(); k++)
			{
				float templen = m_all4sidePatch.at(i)->_surfaceData.GetDistanceBetweenTwoSurface(m_all4sidePatch.at(idxarr.at(k))->_surfaceData);
				if (templen < minlen)
				{
					minidx = idxarr.at(k);
					minlen = templen;
				}
			}
		}
		m_all4sidePatch.at(i)->_haveNextSuf = true;
		m_all4sidePatch.at(i)->_nextSurface = m_all4sidePatch.at(minidx);
		m_all4sidePatch.at(i)->_pointOther = m_all4sidePatch.at(minidx)->_pointThis;
		m_all4sidePatch.at(minidx)->_haveNextSuf = true;
		m_all4sidePatch.at(minidx)->_nextSurface = m_all4sidePatch.at(i);
		m_all4sidePatch.at(minidx)->_pointOther = m_all4sidePatch.at(i)->_pointThis;
	}
	return true;
}

vector<Cube *> YN::WingStruct::transNurbsVoltoCube(vector<NurbsVol> nurbsVols)
{
	vector<Cube *> m_HexVolumes;
	for (int i = 0; i != nurbsVols.size(); ++i)
	{
		Cube * c = new Cube;
		m_HexVolumes.push_back(c);
	}
	for (int i = 0; i != nurbsVols.size(); i++)
	{
		//u�ڵ�����
		for (int j = 0; j != nurbsVols.at(i)._u_Knots->size(); ++j)
		{
			m_HexVolumes.at(i)->_Vol_Data._v_Knots->push_back(nurbsVols.at(i)._u_Knots->at(j));
		}
		//v�ڵ�����
		for (int j = 0; j != nurbsVols.at(i)._v_Knots->size(); ++j)
		{
			m_HexVolumes.at(i)->_Vol_Data._u_Knots->push_back(nurbsVols.at(i)._v_Knots->at(j));
		}
		//w�ڵ�����
		for (int j = 0; j != nurbsVols.at(i)._w_Knots->size(); ++j)
		{
			m_HexVolumes.at(i)->_Vol_Data._w_Knots->push_back(nurbsVols.at(i)._w_Knots->at(j));
		}
		//uvw����
		m_HexVolumes.at(i)->_Vol_Data._v_Degree = nurbsVols.at(i)._u_Degree;
		m_HexVolumes.at(i)->_Vol_Data._u_Degree = nurbsVols.at(i)._v_Degree;
		m_HexVolumes.at(i)->_Vol_Data._w_Degree = nurbsVols.at(i)._w_Degree;

		m_HexVolumes.at(i)->_Vol_Data._v_Num = nurbsVols.at(i)._u_Num;
		m_HexVolumes.at(i)->_Vol_Data._u_Num = nurbsVols.at(i)._v_Num;
		m_HexVolumes.at(i)->_Vol_Data._w_Num = nurbsVols.at(i)._w_Num;

		//���Ƶ�
		for (int j = 0; j != nurbsVols.at(i)._ControlPts->size(); ++j)
		{
			point3d v;
			v.x = nurbsVols.at(i)._ControlPts->at(j).x;
			v.y = nurbsVols.at(i)._ControlPts->at(j).y;
			v.z = nurbsVols.at(i)._ControlPts->at(j).z;
			m_HexVolumes.at(i)->_Vol_Data._ControlPts->push_back(v);
		}
	}
	return m_HexVolumes;
}

vector<Surface *> YN::WingStruct::upDateCubeSurface(Cube * cube)
{
	vector<Surface *> m_all4sidePatch;
	//u0v0��
	//�����һ��ʹ�����ȴ��� new
	if (cube->_u0v0== NULL)
	{
		cube->_u0v0 = new Surface;
	}
	else
	{
		cube->_u0v0->_surfaceData._ControlPts->clear();
	}
	cube->_u0v0->_surfaceData._v_Knots = cube->_Vol_Data._v_Knots;
	cube->_u0v0->_surfaceData._u_Knots = cube->_Vol_Data._u_Knots;
	cube->_u0v0->_surfaceData._u_Num = cube->_Vol_Data._u_Num;
	cube->_u0v0->_surfaceData._v_Num = cube->_Vol_Data._v_Num;
	cube->_u0v0->_surfaceData._u_Degree = cube->_Vol_Data._u_Degree;
	cube->_u0v0->_surfaceData._v_Degree = cube->_Vol_Data._v_Degree;
	//u0v0�Ŀ��Ƶ�
	for (int j = 0; j != cube->_Vol_Data._v_Num; ++j)
	{
		for (int k = 0; k != cube->_Vol_Data._u_Num; ++k)
		{
			cube->_u0v0->_surfaceData._ControlPts->push_back(cube->_Vol_Data.getControlPoint(k, j, 0));
		}
	}
	m_all4sidePatch.push_back(cube->_u0v0);
	cube->_u0v0->_pointThis= cube;
	cube->_u0v0->_pointThisNumber = u0v0;

	//u1v1��
	if (cube->_u1v1 == NULL)
	{
		cube->_u1v1 = new Surface;
	}
	else
	{
		cube->_u1v1->_surfaceData._ControlPts->clear();
	}
	cube->_u1v1->_surfaceData._v_Knots = cube->_Vol_Data._v_Knots;
	cube->_u1v1->_surfaceData._u_Knots = cube->_Vol_Data._u_Knots;
	cube->_u1v1->_surfaceData._u_Num = cube->_Vol_Data._u_Num;
	cube->_u1v1->_surfaceData._v_Num = cube->_Vol_Data._v_Num;
	cube->_u1v1->_surfaceData._u_Degree = cube->_Vol_Data._u_Degree;
	cube->_u1v1->_surfaceData._v_Degree = cube->_Vol_Data._v_Degree;
	for (int j = 0; j != cube->_Vol_Data._v_Num; ++j)
	{
		for (int k = 0; k != cube->_Vol_Data._u_Num; ++k)
		{
			cube->_u1v1->_surfaceData._ControlPts->push_back(cube->_Vol_Data.getControlPoint(k, j, cube->_Vol_Data._w_Num - 1));
		}
	}
	m_all4sidePatch.push_back(cube->_u1v1);
	cube->_u1v1->_pointThis = cube;
	cube->_u1v1->_pointThisNumber = u1v1;

	//v0w0��
	if (cube->_v0w0 == NULL)
	{
		cube->_v0w0 = new Surface;
	}
	else
	{
		cube->_v0w0->_surfaceData._ControlPts->clear();
	}
	cube->_v0w0->_surfaceData._u_Knots = cube->_Vol_Data._v_Knots;
	cube->_v0w0->_surfaceData._v_Knots = cube->_Vol_Data._w_Knots;
	cube->_v0w0->_surfaceData._u_Num = cube->_Vol_Data._v_Num;
	cube->_v0w0->_surfaceData._v_Num = cube->_Vol_Data._w_Num;
	cube->_v0w0->_surfaceData._u_Degree = cube->_Vol_Data._v_Degree;
	cube->_v0w0->_surfaceData._v_Degree = cube->_Vol_Data._w_Degree;
	//v0w0�Ŀ��Ƶ�
	for (int j = 0; j != cube->_Vol_Data._w_Num; ++j)
	{
		for (int k = 0; k != cube->_Vol_Data._v_Num; ++k)
		{
			cube->_v0w0->_surfaceData._ControlPts->push_back(cube->_Vol_Data.getControlPoint(0, k, j));
		}
	}
	m_all4sidePatch.push_back(cube->_v0w0);
	cube->_v0w0->_pointThis = cube;
	cube->_v0w0->_pointThisNumber = v0w0;

	//v1w1��
	if (cube->_v1w1 == NULL)
	{
		cube->_v1w1 = new Surface;
	}
	else
	{
		cube->_v1w1->_surfaceData._ControlPts->clear();
	}
	cube->_v1w1->_surfaceData._u_Knots = cube->_Vol_Data._v_Knots;
	cube->_v1w1->_surfaceData._v_Knots = cube->_Vol_Data._w_Knots;
	cube->_v1w1->_surfaceData._u_Num = cube->_Vol_Data._v_Num;
	cube->_v1w1->_surfaceData._v_Num = cube->_Vol_Data._w_Num;
	cube->_v1w1->_surfaceData._u_Degree = cube->_Vol_Data._v_Degree;
	cube->_v1w1->_surfaceData._v_Degree = cube->_Vol_Data._w_Degree;
	for (int j = 0; j != cube->_Vol_Data._w_Num; ++j)
	{
		for (int k = 0; k != cube->_Vol_Data._v_Num; ++k)
		{
			cube->_v1w1->_surfaceData._ControlPts->push_back(cube->_Vol_Data.getControlPoint(cube->_Vol_Data._u_Num - 1, k, j));
		}
	}
	m_all4sidePatch.push_back(cube->_v1w1);
	cube->_v1w1->_pointThis = cube;
	cube->_v1w1->_pointThisNumber = v1w1;

	//u0w0��
	if (cube->_u0w0 == NULL)
	{
		cube->_u0w0 = new Surface;
	}
	else
	{
		cube->_u0w0->_surfaceData._ControlPts->clear();
	}
	cube->_u0w0->_surfaceData._u_Knots = cube->_Vol_Data._u_Knots;
	cube->_u0w0->_surfaceData._v_Knots = cube->_Vol_Data._w_Knots;
	cube->_u0w0->_surfaceData._u_Num = cube->_Vol_Data._u_Num;
	cube->_u0w0->_surfaceData._v_Num = cube->_Vol_Data._w_Num;
	cube->_u0w0->_surfaceData._u_Degree = cube->_Vol_Data._u_Degree;
	cube->_u0w0->_surfaceData._v_Degree = cube->_Vol_Data._w_Degree;
	for (int j = 0; j != cube->_Vol_Data._w_Num; ++j)
	{
		for (int k = 0; k != cube->_Vol_Data._u_Num; ++k)
		{
			cube->_u0w0->_surfaceData._ControlPts->push_back(cube->_Vol_Data.getControlPoint(k, 0, j));
		}
	}
	m_all4sidePatch.push_back(cube->_u0w0);
	cube->_u0w0->_pointThis = cube;
	cube->_u0w0->_pointThisNumber = u0w0;

	//u1w1��
	if (cube->_u1w1 == NULL)
	{
		cube->_u1w1 = new Surface;
	}
	else
	{
		cube->_u1w1->_surfaceData._ControlPts->clear();
	}
	cube->_u1w1->_surfaceData._u_Knots = cube->_Vol_Data._u_Knots;
	cube->_u1w1->_surfaceData._v_Knots = cube->_Vol_Data._w_Knots;
	cube->_u1w1->_surfaceData._u_Num = cube->_Vol_Data._u_Num;
	cube->_u1w1->_surfaceData._v_Num = cube->_Vol_Data._w_Num;
	cube->_u1w1->_surfaceData._u_Degree = cube->_Vol_Data._u_Degree;
	cube->_u1w1->_surfaceData._v_Degree = cube->_Vol_Data._w_Degree;
	for (int j = 0; j != cube->_Vol_Data._w_Num; ++j)
	{
		for (int k = 0; k != cube->_Vol_Data._u_Num; ++k)
		{
			cube->_u1w1->_surfaceData._ControlPts->push_back(cube->_Vol_Data.getControlPoint(k, cube->_Vol_Data._v_Num - 1, j));
		}
	}
	m_all4sidePatch.push_back(cube->_u1w1);
	cube->_u1w1->_pointThis = cube;
	cube->_u1w1->_pointThisNumber = u1w1;

	return m_all4sidePatch;
}

vector<Surface *> YN::WingStruct::upDateSurface(vector<Cube *>& m_HexVolumes)
{
	//��������Ϣ
	vector<Surface *> m_all4sidePatch;
	for (int i = 0; i != m_HexVolumes.size(); ++i)
	{
		vector<Surface *> temp = upDateCubeSurface(m_HexVolumes[i]);
		m_all4sidePatch.insert(m_all4sidePatch.end(),temp.begin(),temp.end());
	}
	return m_all4sidePatch;
}

bool YN::WingStruct::Order(vector<Cube *>& m_HexVolumes)
{
	if (m_HexVolumes.empty())
	{
		return false;
	}
	std::list<Cube *> List;
	m_HexVolumes[0]->_isOrdered = true;
	List.push_back(m_HexVolumes[0]);
	while (!List.empty())
	{
		Cube* a = List.front();
		List.pop_front();
		//u0v0
		if (a->_u0v0->_pointOther != NULL && a->_u0v0->_pointOther->_isOrdered != true)
		{
			setNearCubeWithSameUVW(a, a->_u0v0->_pointOther, u0v0);
			a->_u0v0->_pointOther->_isOrdered = true;
			List.push_back(a->_u0v0->_pointOther);
		}
		//u1v1
		if (a->_u1v1->_pointOther != NULL && a->_u1v1->_pointOther->_isOrdered != true)
		{
			setNearCubeWithSameUVW(a, a->_u1v1->_pointOther, u1v1);
			a->_u1v1->_pointOther->_isOrdered = true;
			List.push_back(a->_u1v1->_pointOther);
		}
		//v0w0
		if (a->_v0w0->_pointOther != NULL && a->_v0w0->_pointOther->_isOrdered != true)
		{
			setNearCubeWithSameUVW(a, a->_v0w0->_pointOther, v0w0);
			a->_v0w0->_pointOther->_isOrdered = true;
			List.push_back(a->_v0w0->_pointOther);
		}
		//v1w1
		if (a->_v1w1->_pointOther != NULL && a->_v1w1->_pointOther->_isOrdered != true)
		{
			setNearCubeWithSameUVW(a, a->_v1w1->_pointOther, v1w1);
			a->_v1w1->_pointOther->_isOrdered = true;
			List.push_back(a->_v1w1->_pointOther);
		}
		//u0w0
		if (a->_u0w0->_pointOther != NULL && a->_u0w0->_pointOther->_isOrdered != true)
		{
			setNearCubeWithSameUVW(a, a->_u0w0->_pointOther, u0w0);
			a->_u0w0->_pointOther->_isOrdered = true;
			List.push_back(a->_u0w0->_pointOther);
		}
		//u1w1
		if (a->_u1w1->_pointOther != NULL && a->_u1w1->_pointOther->_isOrdered != true)
		{
			setNearCubeWithSameUVW(a, a->_u1w1->_pointOther, u1w1);
			a->_u1w1->_pointOther->_isOrdered = true;
			List.push_back(a->_u1w1->_pointOther);
		}
	}
	//����������Ƿ񶼾�������,����û�������Ƭ������false
	for (auto iter = m_HexVolumes.begin(); iter != m_HexVolumes.end(); ++iter)
	{
		if ((*iter)->_isOrdered == false)
		{
			return false;
		}
	}
	return true;
}

bool YN::WingStruct::setNearCubeWithSameUVW(Cube * StandarHexCell, Cube *WaitHexCell,int index)
{
	if (StandarHexCell == NULL || WaitHexCell == NULL || index < 0 || index>5)
	{
		return false;
	}
	Surface temp,temp2;
	vector<vector<int>> continor, continor2;
	vector<int> IndxOrder, IndxOrder_old, IndxOrder_orient, FinalOrder;
	vector<point4d> tt, tt1, tt2, tt3;
	int PatchIndx_r;

	//����������ʱ����temp��temp2�������
	switch (index)
	{
	case u0v0:
		temp = *StandarHexCell->_u0v0;
		temp2 = *StandarHexCell->_u0v0->_nextSurface;
		PatchIndx_r = 1;
		break;
	case u1v1:
		temp = *StandarHexCell->_u1v1;
		temp2 = *StandarHexCell->_u1v1->_nextSurface;
		PatchIndx_r = 0;
		break;
	case v0w0:
		temp = *StandarHexCell->_v0w0;
		temp2 = *StandarHexCell->_v0w0->_nextSurface;
		PatchIndx_r = 3;
		break;
	case v1w1:
		temp = *StandarHexCell->_v1w1;
		temp2 = *StandarHexCell->_v1w1->_nextSurface;
		PatchIndx_r = 2;
		break;
	case u0w0:
		temp = *StandarHexCell->_u0w0;
		temp2 = *StandarHexCell->_u0w0->_nextSurface;
		PatchIndx_r = 5;
		break;
	case u1w1:
		temp = *StandarHexCell->_u1w1;
		temp2 = *StandarHexCell->_u1w1->_nextSurface;
		PatchIndx_r = 4;
		break;
	default:
		break;
	}

	//���Ƶ��Ӧ����ģ���ڱ��
	auto findNumber = [](vector<point4d> tt, Cube * cell)
	{
		vector<int> order;
		for (int i = 0; i != tt.size(); ++i)
		{
			for (int j = 0; j != cell->_Vol_Data._ControlPts->size(); ++j)
			{
				if (abs(tt[i].x - cell->_Vol_Data._ControlPts->at(j).x) < 0.001
					&&abs(tt[i].y - cell->_Vol_Data._ControlPts->at(j).y) < 0.001
					&&abs(tt[i].z - cell->_Vol_Data._ControlPts->at(j).z) < 0.001)
				{
					order.push_back(j);
				}
			}
		}
		return order;
	};
	
	//��ȡfinalorder
	auto GetFinalOrder = [&]()
	{
		tt1 = *temp._surfaceData._ControlPts;
		//�ڴ����������Ѱ�ҵ���Щ���Ƶ�,�õ�һ���ɿ��Ƶ��±�
		IndxOrder_old = findNumber(tt1, WaitHexCell);
		

		tt2 = *temp2._surfaceData._ControlPts;
		IndxOrder_orient = findNumber(tt2, WaitHexCell);

		if (tt1.size() != IndxOrder_old.size() || tt2.size() != IndxOrder_orient.size())
		{
			int a = 0;
			a++;
		}

		//finalorder�ļ���
		for (int i = 0; i != IndxOrder_orient.size(); ++i)
		{
			for (int j = 0; j != IndxOrder_old.size(); ++j)
			{
				if (IndxOrder_old[j] == IndxOrder_orient[i])
				{
					FinalOrder.push_back(j);
					break;
				}
			}
		}
	};

	//Ƭ���������п��Ƶ���ȡ��xΪʶ�����
	auto MySwitch = [&](int x, Cube * b)
	{
		vector<vector<int>> c;
		switch (x)
		{

		case u0v0://uv0����
			for (int i = b->_Vol_Data._w_Num-1; i != -1; i--)
			{
				for (int j = 0; j != b->_Vol_Data._v_Num; ++j)
				{
					for (int k = 0; k != b->_Vol_Data._u_Num; ++k)
					{
						tt3.push_back(b->_Vol_Data.getControlPoint(k, j, i));
					}
				}
				c.push_back(findNumber(tt3, b));
				tt3.clear();
			};
			break;
		case u1v1://uv1����
			for (int i = 0; i != b->_Vol_Data._w_Num; i++)
			{
				for (int j = 0; j != b->_Vol_Data._v_Num; ++j)
				{
					for (int k = 0; k != b->_Vol_Data._u_Num; ++k)
					{
						tt3.push_back(b->_Vol_Data.getControlPoint(k, j, i));
					}
				}
				c.push_back(findNumber(tt3, b));
				tt3.clear();
			};
			break;
		case v0w0://vw0����
			for (int i = b->_Vol_Data._u_Num - 1; i != -1; i--)
			{
				for (int j = 0; j != b->_Vol_Data._w_Num; ++j)
				{
					for (int k = 0; k != b->_Vol_Data._v_Num; ++k)
					{
						tt3.push_back(b->_Vol_Data.getControlPoint(i, k, j));
					}
				}
				c.push_back(findNumber(tt3, b));
				tt3.clear();
			};
			break;
		case v1w1://vw1����
			for (int i = 0; i != b->_Vol_Data._u_Num; i++)
			{
				for (int j = 0; j != b->_Vol_Data._w_Num; ++j)
				{
					for (int k = 0; k != b->_Vol_Data._v_Num; ++k)
					{
						tt3.push_back(b->_Vol_Data.getControlPoint(i, k, j));
					}
				}
				c.push_back(findNumber(tt3, b));
				tt3.clear();
			};
			break;
		case u0w0://uw0����
			for (int i = b->_Vol_Data._v_Num - 1; i != -1; i--)
			{
				for (int j = 0; j != b->_Vol_Data._w_Num; ++j)
				{
					for (int k = 0; k != b->_Vol_Data._u_Num; ++k)
					{
						tt3.push_back(b->_Vol_Data.getControlPoint(k, i, j));
					}
				}
				c.push_back(findNumber(tt3, b));
				tt3.clear();
			};
			break;
		case u1w1://uw1����
			for (int i = 0; i != b->_Vol_Data._v_Num; i++)
			{
				for (int j = 0; j != b->_Vol_Data._w_Num; ++j)
				{
					for (int k = 0; k != b->_Vol_Data._u_Num; ++k)
					{
						tt3.push_back(b->_Vol_Data.getControlPoint(k, i, j));
					}
				}
				c.push_back(findNumber(tt3, b));
				tt3.clear();
			};
			break;
		default:
			break;
		}
		return c;
	};

	//��ȡ������������fianlOreder
	GetFinalOrder();
	//��ȡ������������
	continor = MySwitch(PatchIndx_r, StandarHexCell);
	continor2 = MySwitch(temp2._pointThisNumber, WaitHexCell);

	//Ȼ�����finalOrder����һ������
	for (int i = 0; i != continor2.size(); ++i)
	{
		vector<int> c;
		c.resize(continor2[i].size());
		for (int j = 0; j != continor2[i].size(); ++j)
		{
			c[FinalOrder[j]] = continor2[i][j];
			//c.push_back(continor2[i][FinalOrder[j]]);
		}
		continor2[i] = c;
	}

	//���������ģ�ͽ�������
	vector<point4d> replace;
	replace.resize(WaitHexCell->_Vol_Data._ControlPts->size());
	for (int i = 0; i != continor.size(); ++i)
	{
		for (int j = 0; j != continor[i].size(); ++j)
		{
			replace[continor[i][j]] = (*WaitHexCell->_Vol_Data._ControlPts)[continor2[i][j]];
		}
	}
	//����滻
	WaitHexCell->_Vol_Data._ControlPts->clear();
	for (auto iter = replace.begin(); iter != replace.end(); ++iter)
	{
		WaitHexCell->_Vol_Data._ControlPts->push_back(*iter);
	}

	//����滻֮����Ҫ��ԭ�ȵ����ݽ��и���
	//���ȸ�����

	upDateCubeSurface(WaitHexCell);

	//�������֮�����pointOtherָ����Ҫ����

	vector<Surface *> tempSurface;
	vector<Surface *> cubeSurface;
	auto getCubeSurfaces = [&](Surface *surface)
	{
		cubeSurface.push_back(surface);
		if (surface->_nextSurface != NULL)
		{
			tempSurface.push_back(surface->_nextSurface);
			surface->_haveNextSuf = false;
			surface->_pointOther = NULL;
			surface->_nextSurface = NULL;
		}
	};

	getCubeSurfaces(WaitHexCell->_u0v0);
	getCubeSurfaces(WaitHexCell->_u1v1);
	getCubeSurfaces(WaitHexCell->_u0w0);
	getCubeSurfaces(WaitHexCell->_u1w1);
	getCubeSurfaces(WaitHexCell->_v0w0);
	getCubeSurfaces(WaitHexCell->_v1w1);

	for (auto iter = tempSurface.begin(); iter != tempSurface.end(); ++iter)
	{
		for (auto iter1 = cubeSurface.begin(); iter1 != cubeSurface.end(); ++iter1)
		{
			if ((*iter1)->_surfaceData.isTwoSurfaceSame((*iter)->_surfaceData))
			{
				(*iter1)->_haveNextSuf = true;
				(*iter1)->_nextSurface = (*iter);
				(*iter1)->_pointOther = (*iter)->_pointThis;
				break;
			}
		}
	}

	return true;
}

void YN::WingStruct::setGloblControlPointId()
{
	//����������ȫ�ֿ��Ƶ���ţ�֮�󷵻�
	int idx = 0;

	//��ȡ���Ƶ��������
	auto getIndex = [&](int u,int v,int w,int _u_Num,int _v_Num)
	{
		return u + v * _u_Num + w * (_u_Num*_v_Num);
	};
	//����һ������ȫ�ֿ��Ƶ��Ƭ���������ȫ�ֿ��Ƶ����
	auto setSurfaceGolblControlPointID = [&](Cube * cube)
	{
		int u_Num = cube->_Vol_Data._u_Num;
		int v_Num = cube->_Vol_Data._v_Num;
		int w_Num = cube->_Vol_Data._w_Num;
		int index;
		//u0v0�����
		for (int j = 0; j != cube->_Vol_Data._v_Num; ++j)
		{
			for (int k = 0; k != cube->_Vol_Data._u_Num; ++k)
			{
				index=cube->_volGloblControlPointID.at(getIndex(k,j,0,u_Num, v_Num));
				cube->_u0v0->_surfaceGloblControlPointID.push_back(index);
			}
		}
		//u1v1�����
		for (int j = 0; j != cube->_Vol_Data._v_Num; ++j)
		{
			for (int k = 0; k != cube->_Vol_Data._u_Num; ++k)
			{
				index = cube->_volGloblControlPointID.at(getIndex(k, j, w_Num - 1, u_Num, v_Num));
				cube->_u1v1->_surfaceGloblControlPointID.push_back(index);
			}
		}
		//v0w0�����
		for (int j = 0; j != cube->_Vol_Data._w_Num; ++j)
		{
			for (int k = 0; k != cube->_Vol_Data._v_Num; ++k)
			{
				index = cube->_volGloblControlPointID.at(getIndex(0, k, j, u_Num, v_Num));
				cube->_v0w0->_surfaceGloblControlPointID.push_back(index);
			}
		}
		//v1w1�����
		for (int j = 0; j != cube->_Vol_Data._w_Num; ++j)
		{
			for (int k = 0; k != cube->_Vol_Data._v_Num; ++k)
			{
				index = cube->_volGloblControlPointID.at(getIndex(u_Num - 1, k, j, u_Num, v_Num));
				cube->_v1w1->_surfaceGloblControlPointID.push_back(index);
			}
		}
		//u0w0�����
		for (int j = 0; j != cube->_Vol_Data._w_Num; ++j)
		{
			for (int k = 0; k != cube->_Vol_Data._u_Num; ++k)
			{
				index = cube->_volGloblControlPointID.at(getIndex(k, 0, j, u_Num, v_Num));
				cube->_u0w0->_surfaceGloblControlPointID.push_back(index);
			}
		}
		//u1w1�����
		for (int j = 0; j != cube->_Vol_Data._w_Num; ++j)
		{
			for (int k = 0; k != cube->_Vol_Data._u_Num; ++k)
			{
				index = cube->_volGloblControlPointID.at(getIndex(k, v_Num - 1, j, u_Num, v_Num));
				cube->_u1w1->_surfaceGloblControlPointID.push_back(index);
			}
		}
	};
	//����������Ƶ�ȫ����ŷ���cube��
	auto setControlPointID = [&](Surface * surface,Cube * cube)
	{
		if (surface->_haveNextSuf &&
			!surface->_nextSurface->_surfaceGloblControlPointID.empty())
		{
			for (int i = 0; i != surface->_nextSurface->_surfaceData._ControlPts->size(); ++i)
			{
				for (int j = 0; j != cube->_Vol_Data._ControlPts->size(); ++j)
				{
					if ((*surface->_nextSurface->_surfaceData._ControlPts)[i] ==
						(*cube->_Vol_Data._ControlPts)[j])
					{
						cube->_volGloblControlPointID[j] = surface->_nextSurface->_surfaceGloblControlPointID.at(i);
					}
				}
			}
		}
	};
	//����һ��Ƭ����Ƭ�ڿ��Ƶ�ȫ�ֱ��
	auto setVolGloblControlPointID = [&](Cube * cube)
	{
		cube->_volGloblControlPointID.resize(cube->_Vol_Data._u_Num*cube->_Vol_Data._v_Num*cube->_Vol_Data._w_Num);
		std::fill(cube->_volGloblControlPointID.begin(), cube->_volGloblControlPointID.end(), -1);
		setControlPointID(cube->_u0v0, cube);
		setControlPointID(cube->_u1v1, cube);
		setControlPointID(cube->_v0w0, cube);
		setControlPointID(cube->_v1w1, cube);
		setControlPointID(cube->_u0w0, cube);
		setControlPointID(cube->_u1w1, cube);
		for (auto iter = cube->_volGloblControlPointID.begin(); iter != cube->_volGloblControlPointID.end(); ++iter)
		{
			if (*iter == -1)
			{
				*iter = idx++;
			}
		}
	};

	//����ȫ�ֿ��Ƶ�
	std::list<Cube *> List;
	List.push_back(_root);
	_root->_isOrdered = !_root->_isOrdered;
	while (!List.empty())
	{
		Cube * a = List.front();
		List.pop_front();
		setVolGloblControlPointID(a);
		setSurfaceGolblControlPointID(a);
		if (a->_u0v0->_pointOther != NULL && a->_u0v0->_pointOther->_isOrdered != a->_isOrdered)
		{
			List.push_back(a->_u0v0->_pointOther);
			a->_u0v0->_pointOther->_isOrdered = !a->_u0v0->_pointOther->_isOrdered;
		}	
		if (a->_u1v1->_pointOther != NULL && a->_u1v1->_pointOther->_isOrdered != a->_isOrdered)
		{
			List.push_back(a->_u1v1->_pointOther);
			a->_u1v1->_pointOther->_isOrdered = !a->_u1v1->_pointOther->_isOrdered;
		}	
		if (a->_v0w0->_pointOther != NULL && a->_v0w0->_pointOther->_isOrdered != a->_isOrdered)
		{
			List.push_back(a->_v0w0->_pointOther);
			a->_v0w0->_pointOther->_isOrdered = !a->_v0w0->_pointOther->_isOrdered;
		}
		if (a->_v1w1->_pointOther != NULL && a->_v1w1->_pointOther->_isOrdered != a->_isOrdered)
		{
			List.push_back(a->_v1w1->_pointOther);
			a->_v1w1->_pointOther->_isOrdered = !a->_v1w1->_pointOther->_isOrdered;
		}
		if (a->_u0w0->_pointOther != NULL && a->_u0w0->_pointOther->_isOrdered != a->_isOrdered)
		{
			List.push_back(a->_u0w0->_pointOther);
			a->_u0w0->_pointOther->_isOrdered = !a->_u0w0->_pointOther->_isOrdered;
		}		
		if (a->_u1w1->_pointOther != NULL && a->_u1w1->_pointOther->_isOrdered != a->_isOrdered)
		{
			List.push_back(a->_u1w1->_pointOther);
			a->_u1w1->_pointOther->_isOrdered = !a->_u1w1->_pointOther->_isOrdered;
		}
	}

}

vector<Cube *> YN::WingStruct::deepFirstSearch(Cube * startCube, int index)
{
	vector<Cube *> temp;
	switch (index)
	{
	case u0v0:
		while (startCube!= NULL)
		{
			temp.push_back(startCube);
			startCube = startCube->_u0v0->_pointOther;
		}
		break;
	case u1v1:
		while (startCube != NULL)
		{
			temp.push_back(startCube);
			startCube = startCube->_u1v1->_pointOther;
		}
		break;
	case v0w0:
		while (startCube != NULL)
		{
			temp.push_back(startCube);
			startCube = startCube->_v0w0->_pointOther;
		}
		break;
	case v1w1:
		while (startCube != NULL)
		{
			temp.push_back(startCube);
			startCube = startCube->_v1w1->_pointOther;
		}
		break;
	case u0w0:
		while (startCube != NULL)
		{
			temp.push_back(startCube);
			startCube = startCube->_u0w0->_pointOther;
		}
		break;
	case u1w1:
		while (startCube != NULL)
		{
			temp.push_back(startCube);
			startCube = startCube->_u1w1->_pointOther;
		}
		break;
	default:
		break;
	}

	return temp;
}

vector<Cube *> YN::WingStruct::bordFirstSearch()
{
	//ÿ�α���֮ǰ����Ƭ��_bordFirstSearch������ͬ��,
	//�ѵ�һƬ_bordFirstSearch����Ϊ����ͬ�ı�ǣ��������ʹ��
	_root->_bordFirstSearch = !_root->_bordFirstSearch;
	vector<Cube *> re;
	std::list<Cube *> List;
	List.push_back(_root);
	
	while (!List.empty())
	{
		Cube * temp = List.front();
		List.pop_front();
		re.push_back(temp);
		if (temp->_u0v0->_pointOther!=NULL&&temp->_u0v0->_pointOther->_bordFirstSearch != temp->_bordFirstSearch)
		{
			temp->_u0v0->_pointOther->_bordFirstSearch = temp->_bordFirstSearch;
			List.push_back(temp->_u0v0->_pointOther);
		}
		if (temp->_u1v1->_pointOther != NULL && temp->_u1v1->_pointOther->_bordFirstSearch != temp->_bordFirstSearch)
		{
			temp->_u1v1->_pointOther->_bordFirstSearch = temp->_bordFirstSearch;
			List.push_back(temp->_u1v1->_pointOther);
		}
		if (temp->_u0w0->_pointOther != NULL && temp->_u0w0->_pointOther->_bordFirstSearch != temp->_bordFirstSearch)
		{
			temp->_u0w0->_pointOther->_bordFirstSearch = temp->_bordFirstSearch;
			List.push_back(temp->_u0w0->_pointOther);
		}
		if (temp->_u1w1->_pointOther != NULL && temp->_u1w1->_pointOther->_bordFirstSearch != temp->_bordFirstSearch)
		{
			temp->_u1w1->_pointOther->_bordFirstSearch = temp->_bordFirstSearch;
			List.push_back(temp->_u1w1->_pointOther);
		}
		if (temp->_v0w0->_pointOther != NULL && temp->_v0w0->_pointOther->_bordFirstSearch != temp->_bordFirstSearch)
		{
			temp->_v0w0->_pointOther->_bordFirstSearch = temp->_bordFirstSearch;
			List.push_back(temp->_v0w0->_pointOther);
		}
		if (temp->_v1w1->_pointOther != NULL && temp->_v1w1->_pointOther->_bordFirstSearch != temp->_bordFirstSearch)
		{
			temp->_v1w1->_pointOther->_bordFirstSearch = temp->_bordFirstSearch;
			List.push_back(temp->_v1w1->_pointOther);
		}
	}
	return re;
}

vector<NurbsSurface> YN::WingStruct::getModelSurface()
{
	vector<Cube *> a;
	vector<NurbsSurface> re;
	a = bordFirstSearch();
	for (auto iter = a.begin(); iter != a.end(); ++iter)
	{
		if ((*iter)->_u0v0->_haveNextSuf == false)
		{
			re.push_back((*iter)->_u0v0->_surfaceData);
		}
		if ((*iter)->_u1v1->_haveNextSuf == false)
		{
			re.push_back((*iter)->_u1v1->_surfaceData);
		}
		if ((*iter)->_v0w0->_haveNextSuf == false)
		{
			re.push_back((*iter)->_v0w0->_surfaceData);
		}
		if ((*iter)->_v1w1->_haveNextSuf == false)
		{
			re.push_back((*iter)->_v1w1->_surfaceData);
		}
		if ((*iter)->_u0w0->_haveNextSuf == false)
		{
			re.push_back((*iter)->_u0w0->_surfaceData);
		}
		if ((*iter)->_u1w1->_haveNextSuf == false)
		{
			re.push_back((*iter)->_u1w1->_surfaceData);
		}
	}
	return re;
}

bool YN::WingStruct::removeVol(point4d target)
{
	//ͨ�����������ҵ�Ŀ��Ƭ
	vector<Cube *> a = quickSearch(target, _root);
	if (a.empty())
	{
		return false;
	}
	Cube * b = a.back();

	//ɾ����
	auto deleteSurface = [&](Surface * surface)
	{
		if (surface != NULL)
		{
			surface->_nextSurface->_pointOther = nullptr;
			surface->_nextSurface->_nextSurface = nullptr;
			surface->_nextSurface->_haveNextSuf = false;
			delete surface;
			surface = nullptr;
		}
	};
	//������ɾ��
	if (b->_u0v0->_haveNextSuf)
	{
		deleteSurface(b->_u0v0);
	}
	if (b->_u1v1->_haveNextSuf)
	{
		deleteSurface(b->_u1v1);
	}
	if (b->_v0w0->_haveNextSuf)
	{
		deleteSurface(b->_v0w0);
	}
	if (b->_v1w1->_haveNextSuf)
	{
		deleteSurface(b->_v1w1);
	}
	if (b->_u0w0->_haveNextSuf)
	{
		deleteSurface(b->_u0w0);
	}
	if (b->_u1w1->_haveNextSuf)
	{
		deleteSurface(b->_u1w1);
	}
	//Ƭɾ��
	delete b;
	b = nullptr;

	return true;
}

bool YN::WingStruct::addVolToSurface(point4d target, NurbsLine line, NurbsSurface surface)
{
	//���������ҵ�Ŀ��Ƭ
	vector<Cube *> a = quickSearch(target, _root);
	NurbsVol vol;
	Cube * b = a.back();
	int surfaceNumber = isPointInSurface(target, b);
	if (surfaceNumber == -1)
		return false;
	//������������nurbsvol
	switch (surfaceNumber)
	{
	case 0:
		vol.LoftingNurbsVol(line, surface, b->_u0v0->_surfaceData);
		break;
	case 1:
		vol.LoftingNurbsVol(line, surface, b->_u1v1->_surfaceData);
		break;
	case 2:
		vol.LoftingNurbsVol(line, surface, b->_v0w0->_surfaceData);
		break;
	case 3:
		vol.LoftingNurbsVol(line, surface, b->_v1w1->_surfaceData);
		break;
	case 4:
		vol.LoftingNurbsVol(line, surface, b->_u0w0->_surfaceData);
		break;
	case 5:
		vol.LoftingNurbsVol(line, surface, b->_u1w1->_surfaceData);
		break;
	default:
		break;
	}

	Cube * newCube = new Cube;
	newCube->_Vol_Data = vol;
	upDateCubeSurface(newCube);

	auto function = [&](Surface * surface)
	{
		surface->_haveNextSuf = true;
		surface->_pointOther = newCube;
		surface->_nextSurface = newCube->_u1v1;
		setNearCubeWithSameUVW(b, newCube, surfaceNumber);
	};

	switch (surfaceNumber)
	{
	case 0:
		function(b->_u0v0);
		break;
	case 1:
		function(b->_u1v1);
		break;
	case 2:
		function(b->_v0w0);
		break;
	case 3:
		function(b->_v1w1);
		break;
	case 4:
		function(b->_u0w0);
		break;
	case 5:
		function(b->_u1w1);
		break;
	default:
		break;
	}

	return true;
}

int YN::WingStruct::isPointInSurface(point4d target, Cube * cube)
{
	auto compare = [&](Surface * surface)
	{
		for (auto i = surface->_surfaceData._ControlPts->begin();
			i != surface->_surfaceData._ControlPts->end(); ++i)
		{
			if (*i == target)
				return surface->_pointThisNumber;
		}
		return -1;
	};
	
	int re;
	re = compare(cube->_u0v0);
	if (re != -1)
	{
		return re;
	}	
	
	re = compare(cube->_u1v1);
	if (re != -1)
	{
		return re;
	}
	
	re = compare(cube->_v0w0);
	if (re != -1)
	{
		return re;
	}

	re = compare(cube->_v1w1);
	if (re != -1)
	{
		return re;
	}

	re = compare(cube->_u0w0);
	if (re != -1)
	{
		return re;
	}

	re = compare(cube->_u1w1);
	if (re != -1)
	{
		return re;
	}

	return -1;
}

vector<Cube *> YN::WingStruct::quickSearch(point4d target, Cube * cubes)
{
	//Ŀǰ�ĳ��򻹲��ܴﵽlogN��ƽ��ʱ�临�Ӷȣ�����Ҫ���иĽ�
	//ʹ��RTT����������ķ��������������������и���һ�����򣬷�ֹ����Ŀ������
	//�߼����������߸������򣬵�������Ŀ��ʱ��ȷ���������߶ȣ�������ָ��̵�·��������϶̵�
	//�߶ȣ��Դ��ҳ���̵�·��

	struct Node;
	vector<Cube *> re;
	std::stack<Node *> Stack;

	//�������ڵ�
	struct Node
	{
		Node * father;
		Cube * Data;
		vector<Node *> childs;
	};

	//�Ӹ��ڵ�root��ʼ
	//���������,������������Ϊ�ӽڵ�
	Node * head = new Node;
	if (/*�ཻ����������Ƿ���Ƭ���ڲ���������ڲ�)*/true)
	{
		//���·���ĳ��ȣ������е�·���������Ա�
		int size;

	}
	//��������ڲ��� �����ܵ��� ����Ӧ��Ƭȫ������Ϊ�ӽڵ�
	else
	{
		vector<Surface *> maybeSurfce;
		for (auto iter = maybeSurfce.begin(); iter != maybeSurfce.end(); ++iter)
		{
			Node * newNode = new Node;
			newNode->Data = (*iter)->_pointOther;
			newNode->father = head;
			head->childs.push_back(newNode);
		}
	}

	


	return re;
}

void YN::WingStruct::Modify(Cube* target, int direct, int u, int v, int w)
{
	//�޸İ����˿��Ƶ��ϸ��
	//UVW���� 0 0 0 �ķ�ʽ��¼
	//���������������Ҫ��������Ϊ UW ���� ������
	//U V W
	//1 0 1
	//�������Ϊ5������֮�󰤸����� �� & ����

	vector<Cube *> uDirect, vDirect, wDirect,tempDirect,totalDirect;

	//���U��Ϊ0
	if (direct & 4 != 0)
	{
		tempDirect = deepFirstSearch(target, 0);
		uDirect.insert(uDirect.end(), tempDirect.begin() + 1, tempDirect.end());
		tempDirect = deepFirstSearch(target, 1);
		uDirect.insert(uDirect.end(), tempDirect.begin() + 1, tempDirect.end());
	}
	//���V��Ϊ0
	if (direct & 2 != 0)
	{
		for (auto iter = uDirect.begin(); iter != uDirect.end(); ++iter)
		{
			tempDirect = deepFirstSearch(*iter, 2);
			vDirect.insert(vDirect.end(), tempDirect.begin() + 1, tempDirect.end());
			tempDirect = deepFirstSearch(*iter, 3);
			vDirect.insert(vDirect.end(), tempDirect.begin() + 1, tempDirect.end());
		}
	}
	//���w��Ϊ0
	if (direct & 4 != 0)
	{
		for (auto iter = vDirect.begin(); iter != vDirect.end(); ++iter)
		{
			tempDirect = deepFirstSearch(*iter, 4);
			wDirect.insert(wDirect.end(), tempDirect.begin() + 1, tempDirect.end());
			tempDirect = deepFirstSearch(*iter, 5);
			wDirect.insert(wDirect.end(), tempDirect.begin() + 1, tempDirect.end());
		}
	}

	for (auto iter = wDirect.begin(); iter != wDirect.end(); ++iter)
	{
		int U = max((*iter)->_Vol_Data._u_Degree, u);
		int V = max((*iter)->_Vol_Data._v_Degree, v);
		int W = max((*iter)->_Vol_Data._w_Degree, w);
		(*iter)->_Vol_Data.DegreeElevate(U, V, W);
		upDateCubeSurface((*iter));
	}
}