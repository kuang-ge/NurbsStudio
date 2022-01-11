//  Real-time Simulation Library
//	str.h: �X�g�����O�N���X�錾��
//
//  	Copyright(C) 2001,  T. Ishii
//////////////////////////////////////////////////////////////////////

#ifndef __STRING_H_INCLUDED__141
#define __STRING_H_INCLUDED__141
#include <tchar.h>

namespace base 
{

//String
// ��������������߂̃N���X
class  String {
//�I�y���[�V����
public:
	
	String();//�f�t�H���g�R���X�g���N�^
	String(const String& str); //��������R�s�[���ăR���X�g���N�g
	String(int val);//�R���X�g���N�^
	String(float  val);//�R���X�g���N�^	
	String(double val);//�R���X�g���N�^
	String(const _TCHAR* format, ...);//�R���X�g���N�^
	String(_TCHAR c);//�R���X�g���N�^
#ifdef _UNICODE
	String(const std::wstring&);//�R���X�g���N�^
#else
	//String(const std::string&);//�R���X�g���N�^
#endif
	
	
	~String();//�f�X�g���N�^


	//������̃N���A
	void clear();

	//�����o�b�t�@�ɃA�N�Z�X
	_TCHAR* getbuf() const;
	
	//�����̈�̋����m��
	void reserve(int size);
	
	//������̒����𓾂�
	int len() const;
	
	const _TCHAR& at(int i) const;//i�Ԗڂ̕����𓾂�
	_TCHAR& at(int i) {return m_buf[i];}//i�Ԗڂ̕����𓾂�
	
	//i�Ԗڂ̕�����c����
	void set(int i, _TCHAR c);
	
	//format�ɏ]����������
	String& format(const _TCHAR* format, ...);
	
	//format�ɏ]��������ǉ�
	String& add_format(const _TCHAR* format, ...);
	
	//�����̌���
	//������Ȃ������ꍇ��-1��Ԃ�
	//���������ꍇ�͂��̏ꏊ��Ԃ�
	//offset�͌����J�n�ʒu
	int  find(const _TCHAR c, int offset=0) const;

	//������̌���
	//������Ȃ������ꍇ��-1��Ԃ�
	//���������ꍇ�͂��̏ꏊ��Ԃ�
	//offset�͌����J�n�ʒu
	int  find(const _TCHAR* str, int offset=0) const;

	//�ꕶ���u��
	//all=true�̏ꍇ�́A�S�̂�u�� false�̏ꍇ�͐擪�݂̂�u��
	//offset�͌����J�n�ʒu
	String replace(const _TCHAR from, const _TCHAR to, int all=1, int offset=0);

	//����������u��
	//all=true�̏ꍇ�́A�S�̂�u�� false�̏ꍇ�͐擪�݂̂�u��
	//offset�͌����J�n�ʒu
	String replace(const _TCHAR* from, const _TCHAR* to, int all=1, int offset=0);
	
	//c�Ɠ������������ׂč폜
	String delchar(const _TCHAR* c, int offset=0);

	//������̍�����white�ɓo�^����Ă��镶�������ׂč폜
	String delleftwhite(const _TCHAR* white);

	//���ԕ�����𓾂�
	String mid(int from, int len) const;
	//�����ԕ�����𓾂�
	String left(int len) const;
	//�E���ԕ�����𓾂�
	String right(int len) const;
	//�T�u�X�g�����O�̒��o
	String sub(int from, int len) const;
	//���ׂĂ�啶���ɓ���
	String upper() const;
	//���ׂĂ��������ɓ���
	String lower() const;
	
	//�f���~�^��p���ĒP��̐؂�o�����s��
	//�؂�o�����P��͌��̕����񂩂�폜�����
	String getword(const _TCHAR* delm);
	
	//int�^�ɕϊ��\���H
	int isint()   const;

	//float�^�ɕϊ��\���H
	int isfloat() const;

	//const char*�^�ɕϊ�
	operator const _TCHAR*() const;

	//int�^�ɕϊ�
	int    atoi() const;
	//float�^�ɕϊ�
	double atof() const;

	String& operator=(const _TCHAR* a);
	String& operator=(const String& a);
#ifdef _UNICODE
	String& operator=(const std::wstring&);
#else
	//String& operator=(const std::string&);
#endif
	
	String& operator=(int   a);
	String& operator=(float a);
	String& operator=(double a);
	String& operator+=(const _TCHAR* a);
	String& operator+=(const _TCHAR a);
	String& operator+=(int   a);
	String& operator+=(float a);
	String& operator+=(double a);
#ifdef _UNICODE
	String& operator+=(const std::wstring& str);
#else
	//String& operator+=(const std::string& str);
#endif
	

	friend   int operator==(const String& a, const String& b);
	friend   int operator==(const String& a, const _TCHAR* b);
	friend   int operator!=(const String& a, const String& b);
	friend   int operator!=(const String& a, const _TCHAR* b);
	friend   String operator+(const String& a, const String& b);
	friend   String operator+(const String& a, const _TCHAR* b);
	friend   String operator+(const _TCHAR* b, const String& a);
	friend   String operator+(const String& a, const _TCHAR b);
	friend   String operator+(const _TCHAR b, String& a);

	friend   String operator+(const String& a, const int& b);
	friend   String operator+(const String& a, const float& b);
	friend   String operator+(const String& a, const double& b);
	friend   bool   operator>(const String& a, const String& b);
	friend   bool   operator<(const String& a, const String& b);

	bool equals(const String& a);
//�A�g���r���[�g
private:
	_TCHAR* m_buf;
	int   m_allocsize;
};
#if defined(__BASE_STREAM_H_INCLUDED__)
/* �L���X�goperator ������̂Œ�`���Ȃ��ėǂ�
BaseOStream& operator<<(BaseOStream& ostr, const String& str) 
{
	ostr << (const char*)str;
	return ostr;
}
*/
inline
BaseIStream& operator>>(BaseIStream& istr, String& str) 
{
	if (istr.ascii()) {
		_TCHAR tmp[1024];
		istr.get_token(tmp, 1024);
		str = tmp;
	} else {
		DWORD len;
		istr >> len;
		std::vector<_TCHAR> tmp((len+1)*sizeof(_TCHAR), '\0');
		istr.read(&(*tmp.begin()), len);

		str = (_TCHAR*)(&(*tmp.begin()));	// modified by XXX 2004.7.15
	}
	return istr;
}

inline
BaseOStream& operator<<(BaseOStream& ostr, const String& str) 
{
	if (ostr.ascii()) {
		ostr << (const _TCHAR*)str;
	} else {
		DWORD len = str.len();
		ostr << len;
		ostr.write(str.getbuf(), len);
	}
	return ostr;
}

#endif


}

#endif

