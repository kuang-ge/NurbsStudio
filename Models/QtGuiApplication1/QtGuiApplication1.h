#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_QtGuiApplication1.h"


class QtGuiApplication1 : public QMainWindow
{
	Q_OBJECT

public:
	QtGuiApplication1(QWidget *parent = Q_NULLPTR);

public:
	//�����ļ�����
	void test();

protected:
	//����������ʾ����
	void setDockWidget_data(varray<NurbsVol> vol);
	void setDockWidget_data(varray<NurbsSurface> vol);
	void setDockWidget_data(varray<NurbsLine> vol);

	//����Ȩ��ť
	void authorityOpen();
public slots:
	//˫��������ݷ�Ӧ����
	void slotTest(QTreeWidgetItem * in1, int in2);
	void slotChanged(QTreeWidgetItem * in1, int in2);
private:
	varray<SplineSurface> Spf;
	varray<SplineVolume> Spv;
	varray<Spline> Spl;
	Ui::QtGuiApplication1Class ui;
	Option option;
};
