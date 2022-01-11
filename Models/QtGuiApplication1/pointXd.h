#ifndef _POINTXD_H_
#define _POINTXD_H_ 

#ifndef PI
#define PI 3.14159265358979323846
#endif // !PI

#include <cmath>

class point2d
{
public:
	double x, y;

	point2d();
	point2d(const double x, const double y);

	//��ȡi����ֵ�ĳ�������
	const double& e(int i) const {
		return (&x)[i];
	}

	//��ȡi����ֵ
	double& e(int i) {
		return (&x)[i];
	}

	//����[]ʵ�ֻ�ȡi����ֵ
	double& operator[] (int i) { return (&x)[i]; }

	//����[]ʵ�ֻ�ȡi����ֵ�ĳ�������
	const double& operator[] (int i) const { return (&x)[i]; }

	//���������
	point2d& operator+=(const point2d &a);
	point2d& operator+=(const double &a)
	{
		x += a; y += a;
		return  (*this);
	};
	point2d& operator-=(const point2d &a);
	point2d& operator*=(const double &a);
	point2d& operator/=(const double &a);

	//����ģ����ƽ��
	double SquareMagnitude() const;
	//����ģ��
	double Magnitude() const;

	//������λ��
	point2d Normalize() const;

	// ��С�������
	double Min() const;
	// ����������
	double Max() const;
	// v���������Ƿ����off
	bool Equals(const point2d& v, double off) const;

	//����ȡ����ֵ
	point2d Abs();

	//���������
	friend point2d operator+(const point2d& a, const point2d& b);
	friend point2d operator-(const point2d& a, const point2d& b);
	friend point2d operator/(const point2d& a, const point2d& b);

	friend point2d operator-(const point2d& a);

	friend point2d operator*(const point2d& a, const double& b);
	friend point2d operator/(const point2d& a, const double& b);
	friend point2d operator*(const double& b, const point2d& a);

	friend int operator == (const point2d& v1, const point2d& v2);

	//���
	double Dot(const point2d& a);

	//�������ȡС
	point2d Minimize(const point2d& v1, const point2d& v2);
	//�������ȡ��
	point2d Maximize(const point2d& v1, const point2d& v2);

	//ƽ��
	point2d Trans(double dx, double dy);

	//�������нǣ��ȣ�
	double Angle(const point2d p);
};

inline
point2d::point2d()
{
	x = 0; y = 0;
}

inline
point2d::point2d(double _x, double _y)
{
	x = _x; y = _y;
}

inline point2d&
point2d::operator+=(const point2d &a)
{
	x = x + a.x; y = y + a.y;
	return  (*this);
}

inline point2d&
point2d::operator-=(const point2d &a)
{
	x = x - a.x; y = y - a.y;
	return  (*this);
}

inline point2d&
point2d::operator*=(const double &a)
{
	x = a * x; y = a * y;
	return  (*this);
}

inline point2d&
point2d::operator/=(const double &a)
{
	x = x / a; y = y / a;
	return  (*this);
}

//����ƽ��
inline double
point2d::SquareMagnitude()  const
{
	return x*x + y*y;
}

//����
inline double
point2d::Magnitude()  const
{
	return sqrt(SquareMagnitude());
}


//��λ��
inline point2d
point2d::Normalize() const
{
	double len = Magnitude();
	point2d nor(0.0, 0.0);
	if (len != 0)
	{
		nor.x = x / len;
		nor.y = y / len;
	}
	return nor;
}

// ��С�������
inline double
point2d::Min()  const
{
	double ret = x;
	if (y < ret) ret = y;
	return ret;
}

// ����������
inline double
point2d::Max()  const
{
	double ret = x;
	if (ret < y) ret = y;
	return ret;
}

// v���������Ƿ����off
inline bool
point2d::Equals(const point2d& ver, double off) const
{
	if ((*this - ver).Magnitude() < off) {
		return true;
	}
	else {
		return false;
	}
}

//����ȡ����ֵ
inline point2d
point2d::Abs()
{
	point2d r;
	r.x = fabs(x);
	r.y = fabs(y);
	return r;
}

inline point2d
operator+(const point2d& a, const point2d& b)
{
	point2d ret = a;
	ret += b;
	return ret;
}

inline point2d
operator-(const point2d& a, const point2d& b)
{
	point2d ret = a;
	ret -= b;
	return ret;
}

inline point2d
operator/(const point2d& a, const point2d& b)
{
	point2d ret = a;
	ret -= b;
	return point2d(a.x / b.x, a.y / b.y);
}

inline point2d
operator-(const point2d& a)
{
	return point2d(-a.x, -a.y);
}

inline point2d
operator*(const point2d& a, const double& b)
{
	point2d ret = a;
	ret *= b;
	return ret;
}

inline point2d
operator/(const point2d& a, const double& b)
{
	point2d ret = a;
	ret /= b;
	return ret;
}

inline point2d
operator*(const double& b, const point2d& a)
{
	point2d ret = a;
	ret *= b;
	return ret;
}

//���
inline double
point2d::Dot(const point2d& a)
{
	return a.x*x + a.y*y;
}

//�������ȡС
inline point2d
point2d::Minimize(const point2d& v1, const point2d& v2)
{
	return point2d(v1.x < v2.x ? v1.x : v2.x,
		v1.y < v2.y ? v1.y : v2.y);
}

//�������ȡ��
inline point2d
point2d::Maximize(const point2d& v1, const point2d& v2)
{
	return point2d(v1.x >= v2.x ? v1.x : v2.x,
		v1.y >= v2.y ? v1.y : v2.y);
}

inline int
operator == (const point2d& v1, const point2d& v2)
{
	return v1.x == v2.x && v1.y == v2.y;
}


// Following is added by D. XXX - Apr. 15, 2006


//ƽ��
inline
point2d	point2d::Trans(double dx, double dy) {
	point2d v;

	v.x = x + dx;
	v.y = y + dy;
	return v;
}

inline double point2d::Angle(const point2d p)
{
	double a, b, c;
	a = this->Dot(p);
	b = this->Magnitude();
	c = p.Magnitude();
	if (b == 0 || c == 0)
		return -1;

	double d = a / (b*c);
	double res = acos(d) * 180 / PI;
	return res;
}


class  point3d {
public:
	double x, y, z;
	double U, V, W, Gray;
	bool boundary;
public:
	//Ĭ�Ϲ���
	point3d();
	//����
	point3d(double _x, double _y, double _z);

	point3d(point2d p2d);

	//��������
	void Set(double _x, double _y, double _z);

	//����
	void Zero();

	//��ȡi����ֵ�ĳ�������
	const double& e(int i) const {
		return (&x)[i];
	}
	//��ȡi����ֵ
	double& e(int i) {
		return (&x)[i];
	}

	//����[]ʵ�ֻ�ȡi����ֵ
	double& operator[] (int i) { return (&x)[i]; }

	//����[]ʵ�ֻ�ȡi����ֵ�ĳ�������
	const double& operator[] (int i) const { return (&x)[i]; }

	//���������
	point3d& operator+=(const point3d &a);
	point3d& operator+=(const double &a)
	{
		x += a; y += a; z += a;
		return  (*this);
	};
	point3d& operator-=(const point3d &a);
	point3d& operator*=(const double &a);
	point3d& operator/=(const double &a);

	//����ģ����ƽ��
	double SquareMagnitude() const;
	//����ģ��
	double Magnitude() const;

	//������λ��
	point3d Normalize() const;

	// ��С�������
	double Min() const;
	// ����������
	double Max() const;
	// v���������Ƿ����off
	bool Equals(const point3d& v, double off) const;

	//����ȡ����ֵ
	point3d Abs();

	//���������
	friend point3d operator+(const point3d& a, const point3d& b);
	friend point3d operator-(const point3d& a, const point3d& b);
	friend point3d operator/(const point3d& a, const point3d& b);

	friend point3d operator-(const point3d& a);

	friend point3d operator*(const point3d& a, const double& b);
	friend point3d operator/(const point3d& a, const double& b);
	friend point3d operator*(const double& b, const point3d& a);

	friend int operator == (const point3d& v1, const point3d& v2);
	friend int operator != (const point3d& v1, const point3d& v2);

	//���
	double Dot(const point3d& a)const;

	//���
	point3d Cross(const point3d& a)const;

	//�������ȡС
	point3d Minimize(const point3d& v1, const point3d& v2);
	//�������ȡ��
	point3d Maximize(const point3d& v1, const point3d& v2);

	// Added by D. XXX - Apr. 15, 2006

	//ƽ��
	point3d Trans(double dx, double dy, double dz);
	//��ת
	point3d RotateX(double ang);
	point3d RotateY(double ang);
	point3d RotateZ(double ang);
	point3d RotAxis(point3d vect, double ang);

	// ������xozƽ�棬yozƽ��ļн�
	void ComputeGlbAng(double *xrot, double *yrot);

	//�������нǣ��ȣ�
	double Angle(const point3d p);

};

inline
point3d::point3d()
{
	x = 0; y = 0; z = 0;
}

inline
point3d::point3d(double _x, double _y, double _z)
{
	x = _x; y = _y; z = _z;
}

inline point3d::point3d(point2d p2d)
{
	x = p2d.x;
	y = p2d.y;
	z = 0;
}

inline void
point3d::Set(double _x, double _y, double _z)
{
	x = _x; y = _y; z = _z;
}

//0����
inline void
point3d::Zero()
{
	x = 0.0; y = 0.0; z = 0.0;
}

inline point3d&
point3d::operator+=(const point3d &a)
{
	x = x + a.x; y = y + a.y; z = z + a.z;
	return  (*this);
}

inline point3d&
point3d::operator-=(const point3d &a)
{
	x = x - a.x; y = y - a.y; z = z - a.z;
	return  (*this);
}

inline point3d&
point3d::operator*=(const double &a)
{
	x = a * x; y = a * y; z = a * z;
	return  (*this);
}

inline point3d&
point3d::operator/=(const double &a)
{
	x = x / a; y = y / a; z = z / a;
	return  (*this);
}

//����ƽ��
inline double
point3d::SquareMagnitude()  const
{
	return x*x + y*y + z*z;
}

//����
inline double
point3d::Magnitude()  const
{
	return sqrt(SquareMagnitude());
}


//��λ��
inline point3d
point3d::Normalize() const
{
	double len = Magnitude();
	point3d nor(0.0, 0.0, 0.0);
	if (len != 0)
	{
		nor.x = x / len;
		nor.y = y / len;
		nor.z = z / len;
	}
	return nor;
}

// ��С�������
inline double
point3d::Min()  const
{
	double ret = x;
	if (y < ret) ret = y;
	if (z < ret) ret = z;
	return ret;
}

// ����������
inline double
point3d::Max()  const
{
	double ret = x;
	if (ret < y) ret = y;
	if (ret < z) ret = z;
	return ret;
}

// v���������Ƿ����off
inline bool
point3d::Equals(const point3d& ver, double off) const
{
	if ((*this - ver).Magnitude() < off) {
		return true;
	}
	else{
		return false;
	}
}

//����ȡ����ֵ
inline point3d
point3d::Abs()
{
	point3d r;
	r.x = fabs(x);
	r.y = fabs(y);
	r.z = fabs(z);
	return r;
}

inline point3d
operator+(const point3d& a, const point3d& b)
{
	point3d ret = a;
	ret += b;

	return ret;
}

inline point3d
operator-(const point3d& a, const point3d& b)
{
	point3d ret = a;
	ret -= b;

	return ret;
}

inline point3d
operator/(const point3d& a, const point3d& b)
{
	point3d ret = a;
	ret -= b;

	return point3d(a.x / b.x, a.y / b.y, a.z / b.z);
}

inline point3d
operator-(const point3d& a)
{
	return point3d(-a.x, -a.y, -a.z);
}

inline point3d
operator*(const point3d& a, const double& b)
{
	point3d ret = a;
	ret *= b;

	return ret;
}

inline point3d
operator/(const point3d& a, const double& b)
{
	point3d ret = a;
	ret /= b;

	return ret;
}

inline point3d
operator*(const double& b, const point3d& a)
{
	point3d ret = a;
	ret *= b;

	return ret;
}

//���
inline double
point3d::Dot(const point3d& a)const
{
	return a.x*x + a.y*y + a.z*z;
}

//���
inline point3d
point3d::Cross(const point3d& a)const
{
	point3d r;
	r.x = y * a.z - z * a.y;
	r.y = z * a.x - x * a.z;
	r.z = x * a.y - y * a.x;

	return r;
}

//�������ȡС
inline point3d
point3d::Minimize(const point3d& v1, const point3d& v2)
{
	return point3d(v1.x < v2.x ? v1.x : v2.x,
		v1.y < v2.y ? v1.y : v2.y,
		v1.z < v2.z ? v1.z : v2.z);
}

//�������ȡ��
inline point3d
point3d::Maximize(const point3d& v1, const point3d& v2)
{
	return point3d(v1.x >= v2.x ? v1.x : v2.x,
		v1.y >= v2.y ? v1.y : v2.y,
		v1.z >= v2.z ? v1.z : v2.z);
}

inline int
operator == (const point3d& v1, const point3d& v2)
{
	return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

inline int
operator != (const point3d& v1, const point3d& v2)
{
	return v1.x != v2.x || v1.y != v2.y || v1.z != v2.z;
}

// Following is added by D. XXX - Apr. 15, 2006


//ƽ��
inline
point3d	point3d::Trans(double dx, double dy, double dz){
	point3d v;

	v.x = x + dx;
	v.y = y + dy;
	v.z = z + dz;

	return v;
}

//��X����ת
inline
point3d	point3d::RotateX(double ang){
	point3d v;
	double sa, ca;

	sa = (double)sin(ang);
	ca = (double)cos(ang);

	v.x = x;
	v.y = y*ca - z*sa;
	v.z = y*sa + z*ca;
	return v;
}

//��Y����ת
inline
point3d	point3d::RotateY(double ang){
	point3d v;

	double sa, ca;
	sa = (double)sin(ang);
	ca = (double)cos(ang);

	v.x = x*ca + z*sa;
	v.y = y;
	v.z = -x*sa + z*ca;
	return v;
}

//��Z����ת
inline
point3d	point3d::RotateZ(double ang){
	point3d v;

	double sa, ca;
	sa = (double)sin(ang);
	ca = (double)cos(ang);

	v.x = x*ca - y*sa;
	v.y = x*sa + y*ca;
	v.z = z;
	return v;
}

//��������ת
inline
point3d point3d::RotAxis(point3d vect, double ang)
{
	point3d p;
	double xa, ya;

	p = point3d(x, y, z);

	vect.ComputeGlbAng(&xa, &ya);

	p = p.RotateX(xa);
	p = p.RotateY(ya);
	p = p.RotateZ(ang);
	p = p.RotateY(-ya);
	p = p.RotateX(-xa);

	return p;
}

// Compute the global transfomation angle of a vector
// ������xozƽ���yozƽ��ļн�
inline
void point3d::ComputeGlbAng(double *xrot, double *yrot)
{
	double tmp;
	double pi = 3.14159265358979323846;
	double err = 0.000001f;

	// rotate x-axis
	if (fabs(z) < err&&fabs(y) < err) *xrot = 0;
	else if (fabs(z) < err&&y > 0) *xrot = (double)(pi / 2.0);
	else if (fabs(z) < err&&y < 0) *xrot = (double)(-pi / 2.0);
	else{

		if (fabs(z) < err){
			return;
		}

		*xrot = (double)atan(y / z);
	}

	if (z < 0){
		if (y > 0) *xrot = (double)(*xrot + pi);
		else    *xrot = (double)(*xrot - pi);
	}

	// rotate y-axis
	if (fabs(y*sin(*xrot) + z*cos(*xrot)) < err&&fabs(x) < err) *yrot = 0;
	else if (fabs(y*sin(*xrot) + z*cos(*xrot)) < err&&x > 0) *yrot = (double)(-pi / 2.0);
	else if (fabs(y*sin(*xrot) + z*cos(*xrot)) < err&&x < 0) *yrot = (double)(pi / 2.0);
	else{
		tmp = y*(double)sin(*xrot) + z*(double)cos(*xrot);

		if (fabs(tmp) < err){
			return;
		}

		*yrot = (double)atan(-x / tmp);
	}

	if ((y*sin(*xrot) + z*cos(*xrot)) < 0){
		if (x < 0) *yrot = (double)(*yrot + pi);
		else    *yrot = (double)(*yrot - pi);
	}

	return;
}

//�������нǣ��ȣ�
inline double point3d::Angle(const point3d p)
{
	double a, b, c;
	a = this->Dot(p);
	b = this->Magnitude();
	c = p.Magnitude();
	if (b == 0 || c == 0)
		return -1;

	double d = a / (b*c);
	if (d > 1)d = 1;
	else if (d < -1)d = -1;
	double res = acos(d) * 180 / PI;
	return res;
}

//added by Wenbin He in 2018

class point4d : public point3d
{
public:
	double w;

	//Ĭ�Ϲ���
	point4d();
	//����
	point4d(double _x, double _y, double _z, double _w = 1);
	point4d(double _w);

	//��ʼ��
	void Inital();

	//ת��
	point4d(const point3d& p);

	//ת��&����w
	point4d(const point3d& p, double _w);

	point4d RotateX(double ang);
	point4d RotateY(double ang);
	point4d RotateZ(double ang);

	bool operator < (const point4d& p) const
	{
		double erf = 0.001;					//�������ж���ֵ
		if (x - p.x <0 && abs(x - p.x)>erf)
			return true;
		else if (abs(x - p.x) <= erf && y - p.y < 0 && abs(y - p.y)>erf)
			return true;
		else if (abs(x - p.x) <= erf && abs(y - p.y) <= erf && z - p.z <0 && abs(z - p.z)>erf)
			return true;
		return false;
	}
};

inline
point4d::point4d()
{
	w = 1;
}

inline
point4d::point4d(double _x, double _y, double _z, double _w)
{
	x = _x;
	y = _y;
	z = _z;
	w = _w;
}

inline point4d::point4d(double _w)
{
	w = _w;
}

inline
void point4d::Inital()
{
	x = 0;
	y = 0;
	z = 0;
	w = 1;
}

//ת��
inline
point4d::point4d(const point3d& p)
{
	x = p.x;
	y = p.y;
	z = p.z;
	w = 1;
}
inline point4d::point4d(const point3d & p, double _w)
{
	x = p.x;
	y = p.y;
	z = p.z;
	w = _w;
}
inline point4d point4d::RotateX(double ang)
{
	point4d v;
	double sa, ca;

	sa = (double)sin(ang);
	ca = (double)cos(ang);

	v.x = x;
	v.y = y * ca - z * sa;
	v.z = y * sa + z * ca;
	v.w = w;
	return v;
}
inline point4d point4d::RotateY(double ang)
{
	point4d v;

	double sa, ca;
	sa = (double)sin(ang);
	ca = (double)cos(ang);

	v.x = x * ca + z * sa;
	v.y = y;
	v.z = -x * sa + z * ca;
	v.w = w;
	return v;
}
inline point4d point4d::RotateZ(double ang)
{
	point4d v;

	double sa, ca;
	sa = (double)sin(ang);
	ca = (double)cos(ang);

	v.x = x * ca - y * sa;
	v.y = x * sa + y * ca;
	v.z = z;
	v.w = w;
	return v;
}
#endif/* Not def*/