#include "QtGuiApplication1.h"
#include <qfiledialog.h>
#include <qcolordialog.h>
#include <qinputdialog.h>
#include <qprogressdialog.h>
#include <qdialog.h>

#include "MyDoc.h"
#include "sensitiveSettingDialog.h"
#include "modleSettingDialog.h"
#include "fragmentDialog.h"
#include "singlePic.h"
#include "readFileDialog.h"
#include "ChangeFineScaleDialog.h"
#include "camerLocatonDialog.h"
#include "PolyIGA.h"
#include "ReducerDialog.h"
#include "Geardialog.h"

#define nurbsline 3;
#define nurbssurface 4;
#define nurbsvol 5;

QtGuiApplication1::QtGuiApplication1(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	setWindowState(Qt::WindowMaximized);//���������ʾ
	ui.centralWidget->grabKeyboard();//��widget���м��̲�׽������
	MyDoc::Ptr pdoc = MyDoc::getInstance();
	//�ϲ�����dockwidget
	tabifyDockWidget(ui.dockWidget_tools, ui.dockWidget_data);

	//����dockWidget_data�ĳ�ʼ��ʾ
	ui.treeWidget->setHeaderLabels(QStringList() << QString::fromLocal8Bit("��������") << QString::fromLocal8Bit("��ֵ/����"));

	//�����ļ�
	//----------
	connect(ui.actiongear, &QAction::triggered, [=]() {
		//�������ļ�
		ui.centralWidget->releaseKeyboard();
		Geardialog dialog(this);
		dialog.exec();
		ui.centralWidget->grabKeyboard();
		QString path = dialog.getPath();
		varray<SplineVolume> gear = dialog.getGear();
		if (!path.isEmpty()) {
			RWGeometric rwg;
			rwg.WriteSplineVolume(path.toStdString(), gear);
		}
		QOpenGLContext * openglContext = ui.centralWidget->context();//�����Ⱦ����
		option.drawSplineVolume(openglContext, this, gear);
		option.createCP(openglContext);
		pdoc->cleanUp();

	});

	connect(ui.actionReducer, &QAction::triggered, [=]() {
		//�������ļ�
		ui.centralWidget->releaseKeyboard();
		ReducerDialog dialog(this);
		dialog.exec();
		ui.centralWidget->grabKeyboard();
		QString path = dialog.getPath();
		varray<SplineVolume> reducer = dialog.getReducer();
		if (!path.isEmpty()) {
			RWGeometric rwg;
			rwg.WriteSplineVolume(path.toStdString(), reducer);
		}
		QOpenGLContext * openglContext = ui.centralWidget->context();//�����Ⱦ����
		option.drawSplineVolume(openglContext, this, reducer);
		option.createCP(openglContext);
		pdoc->cleanUp();

	});




	connect(ui.actionclear_all, &QAction::triggered, [=]() {
		pdoc->clearScrn();
	});
	
	//��ȡ�ļ�
	//----------
	connect(ui.actionreadFile, &QAction::triggered, [=]() {
		readFileDialog dialog(this);
		dialog.exec();
		QOpenGLContext * openglContext = ui.centralWidget->context();//�����Ⱦ����

		option.OnBnClickedCreateModel(openglContext, this);
		option.createCP(openglContext);//������ʾ����

		//��ʱ�Ĵ�������������������
		//������Ȼ��������⣬�ڴ�������һ˲���ڴ�ռ������Ȼ�����
		pdoc->cleanUp();

		authorityOpen();
	});

	//����3mf�ļ�
	//----------
	connect(ui.actionsaveAs3MF, &QAction::triggered, [=]() {
		//option.OnBnClickedCreateModel();
		option.OnBnClicked3MF();

		//��ʱ�Ĵ�������������������
		//������Ȼ��������⣬�ڴ�������һ˲���ڴ�ռ������Ȼ�����
		MyDoc::Ptr pdoc = MyDoc::getInstance();
		pdoc->cleanUp();
	});

	//����������
	//----------
	connect(ui.actionsenstive, &QAction::triggered, [=]() {
		ui.centralWidget->releaseKeyboard();//�ͷż��̲�׽������
		sensitiveSettingDialog dialog(this);
		dialog.exec();
		ui.centralWidget->grabKeyboard();
	});

	//��ɫ����
	//----------
	connect(ui.actionbordLineColor, &QAction::triggered, [=]() {
		ui.centralWidget->releaseKeyboard();//�ͷż��̲�׽������
		QColor color = QColorDialog::getColor(Qt::white);
		if (color.isValid())
		{
			MyDoc::Ptr pdoc = MyDoc::getInstance();
			//��ɫ��������ת��
			pdoc->_bordLineColor = { (float)color.red() / (float)255,(float)color.green() / (float)255,
				(float)color.blue() / (float)255 };
		}
		ui.centralWidget->grabKeyboard();
	});

	connect(ui.actionsurfaceColor, &QAction::triggered, [=]() {
		ui.centralWidget->releaseKeyboard();//�ͷż��̲�׽������
		QColor color = QColorDialog::getColor(Qt::white);
		if (color.isValid())
		{
			MyDoc::Ptr pdoc = MyDoc::getInstance();
			pdoc->_surfaceColor = { (float)color.red() / 255.0f,(float)color.green() / 255.0f,
				(float)color.blue() / 255.0f };
		}
		ui.centralWidget->grabKeyboard();
	});

	connect(ui.actionControlpointColor, &QAction::triggered, [=]() {
		ui.centralWidget->releaseKeyboard();//�ͷż��̲�׽������
		QColor color = QColorDialog::getColor(Qt::white);
		if (color.isValid())
		{
			MyDoc::Ptr pdoc = MyDoc::getInstance();
			pdoc->_controlPoinrColor = { (float)color.red() / 255.0f,(float)color.green() / 255.0f,
				(float)color.blue() / 255.0f  };
		}
		ui.centralWidget->grabKeyboard();
	});

	connect(ui.actionsliceSurfaceColor, &QAction::triggered, [=]() {
		ui.centralWidget->releaseKeyboard();//�ͷż��̲�׽������
		QColor color = QColorDialog::getColor(Qt::white);
		if (color.isValid())
		{
			MyDoc::Ptr pdoc = MyDoc::getInstance();
			pdoc->_sliceSurfaceColor = { (float)color.red() / 255.0f,(float)color.green() / 255.0f,
				(float)color.blue() / 255.0f  };
		}
		ui.centralWidget->grabKeyboard();
	});

	connect(ui.actionsliceSupportColor, &QAction::triggered, [=]() {
		ui.centralWidget->releaseKeyboard();//�ͷż��̲�׽������
		QColor color = QColorDialog::getColor(Qt::white);
		if (color.isValid())
		{
			MyDoc::Ptr pdoc = MyDoc::getInstance();
			pdoc->_sliceSuppotrColor = { (float)color.red() / 255.0f,(float)color.green() / 255.0f,
				(float)color.blue() / 255.0f  };
		}
		ui.centralWidget->grabKeyboard();
	});

	//�Ƿ���ʾ����ϵ
	//--------
	connect(ui.actioncoordinate, &QAction::triggered, [=]() {
		pdoc->m_showCoordinates = !pdoc->m_showCoordinates;
	});

	//������ʾģʽ���㣬�ߣ���
	//----------
	connect(ui.actionPoints, &QAction::triggered, [=]() {
		pdoc->m_GLMode = 0;
		ui.actionLines->setChecked(false);
		ui.actionsurface->setChecked(false);
	});

	connect(ui.actionLines, &QAction::triggered, [=]() {
		pdoc->m_GLMode = 3;
		ui.actionPoints->setChecked(false);
		ui.actionsurface->setChecked(false);
	});

	connect(ui.actionsurface, &QAction::triggered, [=]() {
		pdoc->m_GLMode = 7;
		ui.actionPoints->setChecked(false);
		ui.actionLines->setChecked(false);
	});

	//�Ƿ���ʾ���Ƶ�
	//----------
	connect(ui.actioncontrolPoints, &QAction::triggered, [=]() {
		pdoc->m_showCtrlPts = !pdoc->m_showCtrlPts;
	});

	//ģ�Ͳ�������
	//----------
	connect(ui.actionmodelSetting, &QAction::triggered, []() {
		modleSettingDialog *dialog = new modleSettingDialog;
		dialog->show();
		dialog->setAttribute(Qt::WA_DeleteOnClose);
	});

	//��Ƭ��ʾ
	//----------
	connect(ui.actionsinglePic, &QAction::triggered, [=]() {
		ui.centralWidget->releaseKeyboard();//�ͷż��̲�׽������
		singlePic *dialog = new singlePic;
		dialog->show();
		dialog->setAttribute(Qt::WA_DeleteOnClose);
		ui.centralWidget->grabKeyboard();
	});

	//��Ƭ
	//----------
	connect(ui.actionFragment, &QAction::triggered, [=]() {
		ui.centralWidget->releaseKeyboard();//�ͷż��̲�׽������
		fragmentDialog dialog(this);
		dialog.exec();
		ui.centralWidget->grabKeyboard();
	});

	//��������num
	//----------
	connect(ui.actionchangeFineScale, &QAction::triggered, [=]() {
		ui.centralWidget->releaseKeyboard();//�ͷż��̲�׽������
		ChangeFineScaleDialog dialog(this);
		dialog.exec();
		QOpenGLContext * openglContext = ui.centralWidget->context();//�����Ⱦ����
		option.OnBnClickedCreateModel(openglContext,this);//�������ݼ���
		option.createCP(openglContext);//������ʾ����

		//��ʱ�Ĵ�������������������
		//������Ȼ��������⣬�ڴ�������һ˲���ڴ�ռ������Ȼ�����
		pdoc->cleanUp();
		ui.centralWidget->grabKeyboard();
	});

	//������������ڵ�λ��
	//----------
	connect(ui.actioncameraLocation, &QAction::triggered, [=]() {
		ui.centralWidget->releaseKeyboard();//�ͷż��̲�׽������
		camerLocatonDialog dialog(this);
		dialog.exec();
		ui.centralWidget->grabKeyboard();
	});

	//����������ʾ��ر�
	//----------
	connect(ui.actiontools, &QAction::triggered, [=]() {
		if (ui.actiontools->isChecked())
		{
			ui.dockWidget_tools->show();
			ui.dockWidget_data->show();
		}
		else
		{
			ui.dockWidget_tools->close();
			ui.dockWidget_data->close();
		}
	});

	//�����ӽǰ�ť
	//connect(ui.initCameraLocationButton, &QPushButton::clicked, [=]() {
	//	pdoc->m_restView = true;
	//});


	//�Ҳ����ݵ����Ӧ
	connect(ui.treeWidget, &QTreeWidget::itemDoubleClicked,this,&QtGuiApplication1::slotTest);

	//�Ҳ����޸�����֮����Ӧ����
	connect(ui.treeWidget, &QTreeWidget::itemChanged, this, &QtGuiApplication1::slotChanged);
}

//void QtGuiApplication1::authorityOpen()
//{
//	//ԭ���ǻ�ɫ�İ�ť��������ʹ�ã������������֮��Ϳ���ʹ��
//	//---------
//	ui.actionsaveAs3MF->setEnabled(true);
//	ui.actionFragment->setEnabled(true);
//	ui.actionsliceSupportColor->setEnabled(true);
//	ui.actionsliceSurfaceColor->setEnabled(true);
//	ui.actionsinglePic->setEnabled(true);
//}

//void QtGuiApplication1::setDockWidget_data(varray<NurbsVol> vol)
//{
//
//	//����������ʾ����������
//	QList<QTreeWidgetItem* > items;
//	for (int i = 0; i != vol.size(); ++i)
//	{
//		//����һ����Ŀ¼
//		QTreeWidgetItem * item = new QTreeWidgetItem(ui.treeWidget);
//		item->setText(0, QString::fromLocal8Bit("��")+QString::number(i)+ QString::fromLocal8Bit("Ƭ"));
//
//		//������Ŀ¼�������Ŀ¼
//		QList<QTreeWidgetItem *> childlist;
//		//�����������Ƶ�����
//		QTreeWidgetItem * m_uNum = new QTreeWidgetItem(QStringList() << "m_uNum" << QString::number(vol[i].m_uNum));
//		QTreeWidgetItem * m_vNum = new QTreeWidgetItem(QStringList() << "m_vNum" << QString::number(vol[i].m_vNum));
//		QTreeWidgetItem * m_wNum = new QTreeWidgetItem(QStringList() << "m_wNum" << QString::number(vol[i].m_wNum));
//		//������������
//		QTreeWidgetItem * m_uDegree = new QTreeWidgetItem(QStringList() << "m_uDegree" << QString::number(vol[i].m_uDegree));
//		QTreeWidgetItem * m_vDegree = new QTreeWidgetItem(QStringList() << "m_vDegree" << QString::number(vol[i].m_vDegree));
//		QTreeWidgetItem * m_wDegree = new QTreeWidgetItem(QStringList() << "m_wDegree" << QString::number(vol[i].m_wDegree));
//
//		//����ڵ�ʸ��
//		QTreeWidgetItem * m_uKnots = new QTreeWidgetItem(item);
//		QTreeWidgetItem * m_vKnots = new QTreeWidgetItem(item);
//		QTreeWidgetItem * m_wKnots = new QTreeWidgetItem(item);
//
//		m_uKnots->setText(0, "m_uKnots");
//		m_uKnots->setText(1, QString::number(vol[i].m_uKnots.size()));
//		m_vKnots->setText(0, "m_vKnots");
//		m_vKnots->setText(1, QString::number(vol[i].m_vKnots.size()));
//		m_wKnots->setText(0, "m_wKnots");
//		m_wKnots->setText(1, QString::number(vol[i].m_wKnots.size()));
//
//		//����u�ڵ�ʸ��
//		for (int j = 0; j != vol[i].m_uKnots.size(); ++j)
//		{
//			QTreeWidgetItem * mm = new QTreeWidgetItem(m_uKnots);
//			mm->setText(0, QString::number(j));
//			mm->setText(1, QString::number(vol[i].m_uKnots[j]));
//		}
//		//����v�ڵ�ʸ��
//		for (int j = 0; j != vol[i].m_vKnots.size(); ++j)
//		{
//			QTreeWidgetItem * mm = new QTreeWidgetItem(m_vKnots);
//			mm->setText(0, QString::number(j));
//			mm->setText(1, QString::number(vol[i].m_vKnots[j]));
//		}
//		//����w�ڵ�ʸ��
//		for (int j = 0; j != vol[i].m_wKnots.size(); ++j)
//		{
//			QTreeWidgetItem * mm = new QTreeWidgetItem(m_wKnots);
//			mm->setText(0, QString::number(j));
//			mm->setText(1, QString::number(vol[i].m_wKnots[j]));
//		}
//
//		//������Ƶ�
//		QTreeWidgetItem * contralPoints = new QTreeWidgetItem(item);
//		contralPoints->setText(0, "m_CtrlPts");
//		contralPoints->setText(1, QString::number(vol[i].m_CtrlPts.size()));
//		for (int j=0; j != vol[i].m_CtrlPts.size(); ++j)
//		{
//			QTreeWidgetItem * aa = new QTreeWidgetItem(contralPoints);
//			aa->setText(0, QString::number(j));
//			//����������xyz����ֵ
//			QTreeWidgetItem *x = new QTreeWidgetItem(aa);
//			x->setText(0, "x");
//			x->setText(1, QString::number(vol[i].m_CtrlPts[j].x));
//
//			QTreeWidgetItem *y = new QTreeWidgetItem(aa);
//			y->setText(0, "y");
//			y->setText(1, QString::number(vol[i].m_CtrlPts[j].y));
//
//			QTreeWidgetItem *z = new QTreeWidgetItem(aa);
//			z->setText(0, "z");
//			z->setText(1, QString::number(vol[i].m_CtrlPts[j].z));
//		}
//
//		childlist.push_back(m_uNum);
//		childlist.push_back(m_vNum);
//		childlist.push_back(m_wNum);
//		childlist.push_back(m_uDegree);
//		childlist.push_back(m_vDegree);
//		childlist.push_back(m_wDegree);
//		childlist.push_back(m_uKnots);
//		childlist.push_back(m_vKnots);
//		childlist.push_back(m_wKnots);
//		childlist.push_back(contralPoints);
//
//		item->addChildren(childlist);
//
//		//�Ѵ����õ�һ����Ŀ¼�Ž�������
//		items.push_back(item);
//	}
//	ui.treeWidget->addTopLevelItems(items);
//}

//void QtGuiApplication1::setDockWidget_data(varray<NurbsSurface> vol)
//{
//	//����������ʾ����������
//	QList<QTreeWidgetItem* > items;
//	for (int i = 0; i != vol.size(); ++i)
//	{
//		//����һ����Ŀ¼
//		QTreeWidgetItem * item = new QTreeWidgetItem(ui.treeWidget);
//		item->setText(0, QString::fromLocal8Bit("��") + QString::number(i) + QString::fromLocal8Bit("Ƭ"));
//
//		//������Ŀ¼�������Ŀ¼
//		QList<QTreeWidgetItem *> childlist;
//		//�����������Ƶ�����
//		QTreeWidgetItem * m_uNum = new QTreeWidgetItem(QStringList() << "m_uNum" << QString::number(vol[i].m_uNum));
//		QTreeWidgetItem * m_vNum = new QTreeWidgetItem(QStringList() << "m_vNum" << QString::number(vol[i].m_vNum));
//		//������������
//		QTreeWidgetItem * m_uDegree = new QTreeWidgetItem(QStringList() << "m_uDegree" << QString::number(vol[i].m_uDegree));
//		QTreeWidgetItem * m_vDegree = new QTreeWidgetItem(QStringList() << "m_vDegree" << QString::number(vol[i].m_vDegree));
//
//		//����ڵ�ʸ��
//		QTreeWidgetItem * m_uKnots = new QTreeWidgetItem(item);
//		QTreeWidgetItem * m_vKnots = new QTreeWidgetItem(item);
//
//		m_uKnots->setText(0, "m_uKnots");
//		m_uKnots->setText(1, QString::number(vol[i].m_uKnots.size()));
//		m_vKnots->setText(0, "m_vKnots");
//		m_vKnots->setText(1, QString::number(vol[i].m_vKnots.size()));
//
//		//����u�ڵ�ʸ��
//		for (int j = 0; j != vol[i].m_uKnots.size(); ++j)
//		{
//			QTreeWidgetItem * mm = new QTreeWidgetItem(m_uKnots);
//			mm->setText(0, QString::number(j));
//			mm->setText(1, QString::number(vol[i].m_uKnots[j]));
//		}
//		//����v�ڵ�ʸ��
//		for (int j = 0; j != vol[i].m_vKnots.size(); ++j)
//		{
//			QTreeWidgetItem * mm = new QTreeWidgetItem(m_vKnots);
//			mm->setText(0, QString::number(j));
//			mm->setText(1, QString::number(vol[i].m_vKnots[j]));
//		}
//
//		//������Ƶ�
//		QTreeWidgetItem * contralPoints = new QTreeWidgetItem(item);
//		contralPoints->setText(0, "m_CtrlPts");
//		contralPoints->setText(1, QString::number(vol[i].m_CtrlPts.size()));
//		for (int j = 0; j != vol[i].m_CtrlPts.size(); ++j)
//		{
//			QTreeWidgetItem * aa = new QTreeWidgetItem(contralPoints);
//			aa->setText(0, QString::number(j));
//			//����������xyz����ֵ
//			QTreeWidgetItem *x = new QTreeWidgetItem(aa);
//			x->setText(0, "x");
//			x->setText(1, QString::number(vol[i].m_CtrlPts[j].x));
//
//			QTreeWidgetItem *y = new QTreeWidgetItem(aa);
//			y->setText(0, "y");
//			y->setText(1, QString::number(vol[i].m_CtrlPts[j].y));
//
//			QTreeWidgetItem *z = new QTreeWidgetItem(aa);
//			z->setText(0, "z");
//			z->setText(1, QString::number(vol[i].m_CtrlPts[j].z));
//		}
//
//		childlist.push_back(m_uNum);
//		childlist.push_back(m_vNum);
//		childlist.push_back(m_uDegree);
//		childlist.push_back(m_vDegree);
//		childlist.push_back(m_uKnots);
//		childlist.push_back(m_vKnots);
//		childlist.push_back(contralPoints);
//
//		item->addChildren(childlist);
//
//		//�Ѵ����õ�һ����Ŀ¼�Ž�������
//		items.push_back(item);
//	}
//	ui.treeWidget->addTopLevelItems(items);
//}

//void QtGuiApplication1::setDockWidget_data(varray<NurbsLine> vol)
//{
//	//����������ʾ����������
//	QList<QTreeWidgetItem* > items;
//	for (int i = 0; i != vol.size(); ++i)
//	{
//		//����һ����Ŀ¼
//		QTreeWidgetItem * item = new QTreeWidgetItem(ui.treeWidget);
//		item->setText(0, QString::fromLocal8Bit("��") + QString::number(i) + QString::fromLocal8Bit("Ƭ"));
//
//		//������Ŀ¼�������Ŀ¼
//		QList<QTreeWidgetItem *> childlist;
//		//������������
//		QTreeWidgetItem * m_uDegree = new QTreeWidgetItem(QStringList() << "m_uDegree" << QString::number(vol[i].m_Degree));
//
//		//����ڵ�ʸ��
//		QTreeWidgetItem * m_uKnots = new QTreeWidgetItem(item);
//
//		m_uKnots->setText(0, "m_uKnots");
//		m_uKnots->setText(1, QString::number(vol[i].m_Knots.size()));
//
//		//����u�ڵ�ʸ��
//		for (int j = 0; j != vol[i].m_Knots.size(); ++j)
//		{
//			QTreeWidgetItem * mm = new QTreeWidgetItem(m_uKnots);
//			mm->setText(0, QString::number(j));
//			mm->setText(1, QString::number(vol[i].m_Knots[j]));
//		}
//
//		//������Ƶ�
//		QTreeWidgetItem * contralPoints = new QTreeWidgetItem(item);
//		contralPoints->setText(0, "m_CtrlPts");
//		contralPoints->setText(1, QString::number(vol[i].m_CtrlPts.size()));
//		for (int j = 0; j != vol[i].m_CtrlPts.size(); ++j)
//		{
//			QTreeWidgetItem * aa = new QTreeWidgetItem(contralPoints);
//			aa->setText(0, QString::number(j));
//			//����������xyz����ֵ
//			QTreeWidgetItem *x = new QTreeWidgetItem(aa);
//			x->setText(0, "x");
//			x->setText(1, QString::number(vol[i].m_CtrlPts[j].x));
//
//			QTreeWidgetItem *y = new QTreeWidgetItem(aa);
//			y->setText(0, "y");
//			y->setText(1, QString::number(vol[i].m_CtrlPts[j].y));
//
//			QTreeWidgetItem *z = new QTreeWidgetItem(aa);
//			z->setText(0, "z");
//			z->setText(1, QString::number(vol[i].m_CtrlPts[j].z));
//		}
//
//		childlist.push_back(m_uDegree);
//		childlist.push_back(m_uKnots);
//		childlist.push_back(contralPoints);
//
//		item->addChildren(childlist);
//
//		//�Ѵ����õ�һ����Ŀ¼�Ž�������
//		items.push_back(item);
//	}
//	ui.treeWidget->addTopLevelItems(items);
//}

//void QtGuiApplication1::slotTest(QTreeWidgetItem * in1, int in2)
//{
//	MyDoc::Ptr pdoc = MyDoc::getInstance();
//	QString text = in1->text(in2);//�����ѡ������
//	QRegExp temp(QString::fromLocal8Bit("��"));
//	if (text.contains(temp))//�鿴�Ƿ�ƥ�䵥����ʾ��ģ��
//	{
//		QRegularExpression exp("\\d{1,2}");
//		QRegularExpressionMatch match = exp.match(text);
//		if (match.hasMatch())
//		{
//			QString matched = match.captured();
//			pdoc->m_CellIdx = matched.toInt();
//		}
//	}
//	if (in1->parent())
//	{
//		//���������ǿ��Ƶ��е�xyz����һ���������ֵ
//		if (in1->parent()->parent())
//		{
//			//�޸�xyz�������ֵ
//			ui.treeWidget->openPersistentEditor(in1, 1);
//		}
//	}
//}

//void QtGuiApplication1::slotChanged(QTreeWidgetItem * in1, int in2)
//{
//	if (option.isLcokOn())
//	{
//		//��ȡ�ı��ֵ
//		QString a = in1->text(in2);
//		ui.treeWidget->closePersistentEditor(in1, in2);
//		MyDoc::Ptr pdoc = MyDoc::getInstance();
//		pdoc->clearScrn();
//		option.m_again = 1;
//		//���Ƶ�ı�Ϊ�ղ��޸ĵ�ֵ
//		//------------
//		//�ҵ�үү�ڵ㣬Ҳ���ǵڼ�����
//		QString matched;
//		QRegularExpression exp("\\d{1,2}");
//		QRegularExpressionMatch match = exp.match(in1->parent()->parent()->parent()->text(0));
//		if (match.hasMatch())
//		{
//			matched = match.captured();
//		}
//		//�ҵ��ǵڼ����ڵ�
//		QString temp = in1->parent()->text(0);
//		//�ҵ�����λ��
//		if (in1->text(0) == "x")
//		{
//			pdoc->m_ShowData4D[matched.toInt()][temp.toInt()].x = a.toFloat();
//		}else		if (in1->text(0) == "y")
//		{
//			pdoc->m_ShowData4D[matched.toInt()][temp.toInt()].y = a.toFloat();
//		}else		if (in1->text(0) == "z")
//		{
//			pdoc->m_ShowData4D[matched.toInt()][temp.toInt()].z = a.toFloat();
//		}
//		//���»���ͼ��
//		QOpenGLContext * openglContext = ui.centralWidget->context();//�����Ⱦ����
//		option.OnBnClickedCreateModel(openglContext, this);
//		option.createCP(openglContext);//������ʾ����
//	}
//}

void QtGuiApplication1::test()
{
	//���Ƶ�ϸ��6*6
	auto ControlPSpecify = [](varray<NurbsVol> &nvs)
	{
		vector<double> knot = { 0.25,0.5,0.75 };
		varray<double> uknots, vknots, wknots;
		for (auto i : knot)
		{
			uknots.push_back(i);
			vknots.push_back(i);
			wknots.push_back(i);
		}
		for (int i = 0; i != nvs.size(); ++i)
		{
			nvs[i].KnotsRefine(uknots, vknots, wknots);
		}
	};

	auto GloblControlP = [](CPolyParaVolume &cppv)
	{
		//����ȫ�ֿ��Ƶ�
		int idx = 0;
		for (idx; idx < cppv.m_HexVolumes[0].m_splVol.m_vAllCtrlPts.size(); idx++) {
			cppv.m_HexVolumes[0].m_controlPtGlobalIDs.push_back(idx);
		}
		cppv.m_HexVolumes[0].m_bGenerateGlobalID = true;
		for (int i = 1; i < cppv.m_HexVolumes.size(); i++) {
			for (int j = 0; j < cppv.m_HexVolumes[i].m_splVol.m_vAllCtrlPts.size(); j++) {
				bool flag = true;
				for (int k = 0; k < i; k++) {
					for (int t = 0; t < cppv.m_HexVolumes[k].m_splVol.m_vAllCtrlPts.size(); t++) {
						if (abs(cppv.m_HexVolumes[i].m_splVol.m_vAllCtrlPts[j].m_pt.x - cppv.m_HexVolumes[k].m_splVol.m_vAllCtrlPts[t].m_pt.x) < ERRF
							&& abs(cppv.m_HexVolumes[i].m_splVol.m_vAllCtrlPts[j].m_pt.y - cppv.m_HexVolumes[k].m_splVol.m_vAllCtrlPts[t].m_pt.y) < ERRF
							&& abs(cppv.m_HexVolumes[i].m_splVol.m_vAllCtrlPts[j].m_pt.z - cppv.m_HexVolumes[k].m_splVol.m_vAllCtrlPts[t].m_pt.z) < ERRF)
						{
							cppv.m_HexVolumes[i].m_controlPtGlobalIDs.push_back(cppv.m_HexVolumes[k].m_controlPtGlobalIDs[t]);
							flag = false;
							break;
						}
					}
					if (flag == false)
						break;
				}
				if (flag == true) {
					cppv.m_HexVolumes[i].m_controlPtGlobalIDs.push_back(idx);
					idx++;
				}
			}
			cppv.m_HexVolumes[i].m_bGenerateGlobalID = true;
		}
	};

	//auto ChangeOrderedVol = [](const varray<NurbsVol> &vol, varray<NurbsVol> &orderVol) {
	//	int patchNum = vol.size();
	//	for (int i = 0; i < patchNum; i++) {
	//		
	//	}
	//};

	QOpenGLContext * openglContext = ui.centralWidget->context();//�����Ⱦ����
	MyDoc::Ptr pdoc = MyDoc::getInstance();


	RWGeometric rwg;
	varray<NurbsVol> nvs;
	rwg.ReadNurbsVol("D:\\quadTest\\��װ�2\\duokongban-vol.txt", nvs);
	for (int i = 0; i != nvs.size(); ++i)
	{
		nvs[i].DegreeElevate(2, 2, 2);
	}
	rwg.WriteNurbsVol("D:\\quadTest\\��װ�2\\duokongban-vol-2degree.txt", nvs);
	//ControlPSpecify(nvs);
	CPolyParaVolume cppv;
	cppv = nvs;
	cppv.SetAdjacentHexIdx();
	cppv.Order();
	GloblControlP(cppv);
	//cppv.OutputParaVolumeDataTxt("test/1.txt", "test/2.txt");
	cppv.OutputParaVolumeDataVTK("D:\\quadTest\\��װ�2\\duokongban-vol.vtk");
	nvs.clear();
	option.TransCnurbsvolToPolyIGA(nvs, cppv);
	rwg.WriteNurbsVol("D:\\quadTest\\��װ�2\\duokongban-vol-ordered.txt", nvs);
	//setDockWidget_data(nvs);

	//�����ļ���ȡ֮��û��ѡȡ��ʾ��ģ�ͣ��ر�ѡȡһ��
	//pdoc->modelType = 5;

	option.drawNurbsVol(openglContext, this, nvs);
	option.OnBnClickedCreateModel(openglContext, this);
	option.createCP(openglContext);//������ʾ����
	//���������ʾ�ĺ��������뱣��
	option.setLockOn();

	////�������ļ�
	//QOpenGLContext * openglContext = ui.centralWidget->context();//�����Ⱦ����
	////�µĲ��Ժ�������д������
	//RWGeometric rwg;
	//varray<NurbsVol> nvs;
	//rwg.ReadNurbsVol("test/new_allvols20_3��������02 - ����.txt", nvs);

	////���Ƶ�ϸ��6*6
	///*vector<double> knot={0.25,0.5,0.75};
	//varray<double> uknots,vknots,wknots;
	//for (auto i : knot)
	//{
	//	uknots.push_back(i);
	//	vknots.push_back(i);
	//	wknots.push_back(i);
	//}
	//for (int i = 0; i != nvs.size(); ++i)
	//{
	//	nvs[i].KnotsRefine(uknots, vknots, wknots);
	//}*/
	//CPolyParaVolume cppv;
	//cppv = nvs;
	//cppv.SetAdjacentHexIdx();
	//cppv.Order();
	//
	////����ȫ�ֿ��Ƶ�
	//int idx = 0;
	//for (idx; idx < cppv.m_HexVolumes[0].m_splVol.m_vAllCtrlPts.size(); idx++) {
	//	cppv.m_HexVolumes[0].m_controlPtGlobalIDs.push_back(idx);
	//}

	//for (int i = 1; i < cppv.m_HexVolumes.size(); i++) {
	//	for (int j = 0; j < cppv.m_HexVolumes[i].m_splVol.m_vAllCtrlPts.size(); j++) {
	//		bool flag = true;
	//		for (int k = 0; k < i; k++) {
	//			for (int t = 0; t < cppv.m_HexVolumes[k].m_splVol.m_vAllCtrlPts.size(); t++) {
	//				if (abs(cppv.m_HexVolumes[i].m_splVol.m_vAllCtrlPts[j].m_pt.x - cppv.m_HexVolumes[k].m_splVol.m_vAllCtrlPts[t].m_pt.x) < ERRF
	//					&& abs(cppv.m_HexVolumes[i].m_splVol.m_vAllCtrlPts[j].m_pt.y - cppv.m_HexVolumes[k].m_splVol.m_vAllCtrlPts[t].m_pt.y) < ERRF
	//					&& abs(cppv.m_HexVolumes[i].m_splVol.m_vAllCtrlPts[j].m_pt.z - cppv.m_HexVolumes[k].m_splVol.m_vAllCtrlPts[t].m_pt.z) < ERRF)
	//				{
	//					cppv.m_HexVolumes[i].m_controlPtGlobalIDs.push_back(cppv.m_HexVolumes[k].m_controlPtGlobalIDs[t]);
	//					flag = false;
	//					break;
	//				}
	//			}
	//			if (flag == false)
	//				break;
	//		}
	//		if (flag == true) {
	//			cppv.m_HexVolumes[i].m_controlPtGlobalIDs.push_back(idx);
	//			idx++;
	//		}
	//	}
	//	cppv.m_HexVolumes[i].m_bGenerateGlobalID = true;
	//}
	//cppv.OutputParaVolumeDataTxt("test/1.txt", "test/2.txt");

	//MyDoc::Ptr pdoc = MyDoc::getInstance();
	////setDockWidget_data(nvs);

	////�����ļ���ȡ֮��û��ѡȡ��ʾ��ģ�ͣ��ر�ѡȡһ��
	//pdoc->modelType = 5;

	//option.drawNurbsVol(openglContext, this, nvs);
	//option.OnBnClickedCreateModel(openglContext, this);
	//option.createCP(openglContext);//������ʾ����
	////���������ʾ�ĺ��������뱣��
	//option.setLockOn();
}

void QtGuiApplication1::authorityOpen()
{
	//ԭ���ǻ�ɫ�İ�ť��������ʹ�ã������������֮��Ϳ���ʹ��
	//---------
	ui.actionsaveAs3MF->setEnabled(true);
	ui.actionFragment->setEnabled(true);
	ui.actionsliceSupportColor->setEnabled(true);
	ui.actionsliceSurfaceColor->setEnabled(true);
	ui.actionsinglePic->setEnabled(true);
}

void QtGuiApplication1::setDockWidget_data(varray<NurbsVol> vol)
{

	//����������ʾ����������
	QList<QTreeWidgetItem* > items;
	for (int i = 0; i != vol.size(); ++i)
	{
		//����һ����Ŀ¼
		QTreeWidgetItem * item = new QTreeWidgetItem(ui.treeWidget);
		item->setText(0, QString::fromLocal8Bit("��") + QString::number(i) + QString::fromLocal8Bit("Ƭ"));

		//������Ŀ¼�������Ŀ¼
		QList<QTreeWidgetItem *> childlist;
		//�����������Ƶ�����
		QTreeWidgetItem * m_uNum = new QTreeWidgetItem(QStringList() << "m_uNum" << QString::number(vol[i].m_uNum));
		QTreeWidgetItem * m_vNum = new QTreeWidgetItem(QStringList() << "m_vNum" << QString::number(vol[i].m_vNum));
		QTreeWidgetItem * m_wNum = new QTreeWidgetItem(QStringList() << "m_wNum" << QString::number(vol[i].m_wNum));
		//������������
		QTreeWidgetItem * m_uDegree = new QTreeWidgetItem(QStringList() << "m_uDegree" << QString::number(vol[i].m_uDegree));
		QTreeWidgetItem * m_vDegree = new QTreeWidgetItem(QStringList() << "m_vDegree" << QString::number(vol[i].m_vDegree));
		QTreeWidgetItem * m_wDegree = new QTreeWidgetItem(QStringList() << "m_wDegree" << QString::number(vol[i].m_wDegree));

		//����ڵ�ʸ��
		QTreeWidgetItem * m_uKnots = new QTreeWidgetItem(item);
		QTreeWidgetItem * m_vKnots = new QTreeWidgetItem(item);
		QTreeWidgetItem * m_wKnots = new QTreeWidgetItem(item);

		m_uKnots->setText(0, "m_uKnots");
		m_uKnots->setText(1, QString::number(vol[i].m_uKnots.size()));
		m_vKnots->setText(0, "m_vKnots");
		m_vKnots->setText(1, QString::number(vol[i].m_vKnots.size()));
		m_wKnots->setText(0, "m_wKnots");
		m_wKnots->setText(1, QString::number(vol[i].m_wKnots.size()));

		//����u�ڵ�ʸ��
		for (int j = 0; j != vol[i].m_uKnots.size(); ++j)
		{
			QTreeWidgetItem * mm = new QTreeWidgetItem(m_uKnots);
			mm->setText(0, QString::number(j));
			mm->setText(1, QString::number(vol[i].m_uKnots[j]));
		}
		//����v�ڵ�ʸ��
		for (int j = 0; j != vol[i].m_vKnots.size(); ++j)
		{
			QTreeWidgetItem * mm = new QTreeWidgetItem(m_vKnots);
			mm->setText(0, QString::number(j));
			mm->setText(1, QString::number(vol[i].m_vKnots[j]));
		}
		//����w�ڵ�ʸ��
		for (int j = 0; j != vol[i].m_wKnots.size(); ++j)
		{
			QTreeWidgetItem * mm = new QTreeWidgetItem(m_wKnots);
			mm->setText(0, QString::number(j));
			mm->setText(1, QString::number(vol[i].m_wKnots[j]));
		}

		//������Ƶ�
		QTreeWidgetItem * contralPoints = new QTreeWidgetItem(item);
		contralPoints->setText(0, "m_CtrlPts");
		contralPoints->setText(1, QString::number(vol[i].m_CtrlPts.size()));
		for (int j = 0; j != vol[i].m_CtrlPts.size(); ++j)
		{
			QTreeWidgetItem * aa = new QTreeWidgetItem(contralPoints);
			aa->setText(0, QString::number(j));
			//����������xyz����ֵ
			QTreeWidgetItem *x = new QTreeWidgetItem(aa);
			x->setText(0, "x");
			x->setText(1, QString::number(vol[i].m_CtrlPts[j].x));

			QTreeWidgetItem *y = new QTreeWidgetItem(aa);
			y->setText(0, "y");
			y->setText(1, QString::number(vol[i].m_CtrlPts[j].y));

			QTreeWidgetItem *z = new QTreeWidgetItem(aa);
			z->setText(0, "z");
			z->setText(1, QString::number(vol[i].m_CtrlPts[j].z));
		}

		childlist.push_back(m_uNum);
		childlist.push_back(m_vNum);
		childlist.push_back(m_wNum);
		childlist.push_back(m_uDegree);
		childlist.push_back(m_vDegree);
		childlist.push_back(m_wDegree);
		childlist.push_back(m_uKnots);
		childlist.push_back(m_vKnots);
		childlist.push_back(m_wKnots);
		childlist.push_back(contralPoints);

		item->addChildren(childlist);

		//�Ѵ����õ�һ����Ŀ¼�Ž�������
		items.push_back(item);
	}
	ui.treeWidget->addTopLevelItems(items);
}

void QtGuiApplication1::setDockWidget_data(varray<NurbsSurface> vol)
{
	//����������ʾ����������
	QList<QTreeWidgetItem* > items;
	for (int i = 0; i != vol.size(); ++i)
	{
		//����һ����Ŀ¼
		QTreeWidgetItem * item = new QTreeWidgetItem(ui.treeWidget);
		item->setText(0, QString::fromLocal8Bit("��") + QString::number(i) + QString::fromLocal8Bit("Ƭ"));

		//������Ŀ¼�������Ŀ¼
		QList<QTreeWidgetItem *> childlist;
		//�����������Ƶ�����
		QTreeWidgetItem * m_uNum = new QTreeWidgetItem(QStringList() << "m_uNum" << QString::number(vol[i].m_uNum));
		QTreeWidgetItem * m_vNum = new QTreeWidgetItem(QStringList() << "m_vNum" << QString::number(vol[i].m_vNum));
		//������������
		QTreeWidgetItem * m_uDegree = new QTreeWidgetItem(QStringList() << "m_uDegree" << QString::number(vol[i].m_uDegree));
		QTreeWidgetItem * m_vDegree = new QTreeWidgetItem(QStringList() << "m_vDegree" << QString::number(vol[i].m_vDegree));

		//����ڵ�ʸ��
		QTreeWidgetItem * m_uKnots = new QTreeWidgetItem(item);
		QTreeWidgetItem * m_vKnots = new QTreeWidgetItem(item);

		m_uKnots->setText(0, "m_uKnots");
		m_uKnots->setText(1, QString::number(vol[i].m_uKnots.size()));
		m_vKnots->setText(0, "m_vKnots");
		m_vKnots->setText(1, QString::number(vol[i].m_vKnots.size()));

		//����u�ڵ�ʸ��
		for (int j = 0; j != vol[i].m_uKnots.size(); ++j)
		{
			QTreeWidgetItem * mm = new QTreeWidgetItem(m_uKnots);
			mm->setText(0, QString::number(j));
			mm->setText(1, QString::number(vol[i].m_uKnots[j]));
		}
		//����v�ڵ�ʸ��
		for (int j = 0; j != vol[i].m_vKnots.size(); ++j)
		{
			QTreeWidgetItem * mm = new QTreeWidgetItem(m_vKnots);
			mm->setText(0, QString::number(j));
			mm->setText(1, QString::number(vol[i].m_vKnots[j]));
		}

		//������Ƶ�
		QTreeWidgetItem * contralPoints = new QTreeWidgetItem(item);
		contralPoints->setText(0, "m_CtrlPts");
		contralPoints->setText(1, QString::number(vol[i].m_CtrlPts.size()));
		for (int j = 0; j != vol[i].m_CtrlPts.size(); ++j)
		{
			QTreeWidgetItem * aa = new QTreeWidgetItem(contralPoints);
			aa->setText(0, QString::number(j));
			//����������xyz����ֵ
			QTreeWidgetItem *x = new QTreeWidgetItem(aa);
			x->setText(0, "x");
			x->setText(1, QString::number(vol[i].m_CtrlPts[j].x));

			QTreeWidgetItem *y = new QTreeWidgetItem(aa);
			y->setText(0, "y");
			y->setText(1, QString::number(vol[i].m_CtrlPts[j].y));

			QTreeWidgetItem *z = new QTreeWidgetItem(aa);
			z->setText(0, "z");
			z->setText(1, QString::number(vol[i].m_CtrlPts[j].z));
		}

		childlist.push_back(m_uNum);
		childlist.push_back(m_vNum);
		childlist.push_back(m_uDegree);
		childlist.push_back(m_vDegree);
		childlist.push_back(m_uKnots);
		childlist.push_back(m_vKnots);
		childlist.push_back(contralPoints);

		item->addChildren(childlist);

		//�Ѵ����õ�һ����Ŀ¼�Ž�������
		items.push_back(item);
	}
	ui.treeWidget->addTopLevelItems(items);
}

void QtGuiApplication1::setDockWidget_data(varray<NurbsLine> vol)
{
	//����������ʾ����������
	QList<QTreeWidgetItem* > items;
	for (int i = 0; i != vol.size(); ++i)
	{
		//����һ����Ŀ¼
		QTreeWidgetItem * item = new QTreeWidgetItem(ui.treeWidget);
		item->setText(0, QString::fromLocal8Bit("��") + QString::number(i) + QString::fromLocal8Bit("Ƭ"));

		//������Ŀ¼�������Ŀ¼
		QList<QTreeWidgetItem *> childlist;
		//������������
		QTreeWidgetItem * m_uDegree = new QTreeWidgetItem(QStringList() << "m_uDegree" << QString::number(vol[i].m_Degree));

		//����ڵ�ʸ��
		QTreeWidgetItem * m_uKnots = new QTreeWidgetItem(item);

		m_uKnots->setText(0, "m_uKnots");
		m_uKnots->setText(1, QString::number(vol[i].m_Knots.size()));

		//����u�ڵ�ʸ��
		for (int j = 0; j != vol[i].m_Knots.size(); ++j)
		{
			QTreeWidgetItem * mm = new QTreeWidgetItem(m_uKnots);
			mm->setText(0, QString::number(j));
			mm->setText(1, QString::number(vol[i].m_Knots[j]));
		}

		//������Ƶ�
		QTreeWidgetItem * contralPoints = new QTreeWidgetItem(item);
		contralPoints->setText(0, "m_CtrlPts");
		contralPoints->setText(1, QString::number(vol[i].m_CtrlPts.size()));
		for (int j = 0; j != vol[i].m_CtrlPts.size(); ++j)
		{
			QTreeWidgetItem * aa = new QTreeWidgetItem(contralPoints);
			aa->setText(0, QString::number(j));
			//����������xyz����ֵ
			QTreeWidgetItem *x = new QTreeWidgetItem(aa);
			x->setText(0, "x");
			x->setText(1, QString::number(vol[i].m_CtrlPts[j].x));

			QTreeWidgetItem *y = new QTreeWidgetItem(aa);
			y->setText(0, "y");
			y->setText(1, QString::number(vol[i].m_CtrlPts[j].y));

			QTreeWidgetItem *z = new QTreeWidgetItem(aa);
			z->setText(0, "z");
			z->setText(1, QString::number(vol[i].m_CtrlPts[j].z));
		}

		childlist.push_back(m_uDegree);
		childlist.push_back(m_uKnots);
		childlist.push_back(contralPoints);

		item->addChildren(childlist);

		//�Ѵ����õ�һ����Ŀ¼�Ž�������
		items.push_back(item);
	}
	ui.treeWidget->addTopLevelItems(items);
}

void QtGuiApplication1::slotTest(QTreeWidgetItem * in1, int in2)
{
	MyDoc::Ptr pdoc = MyDoc::getInstance();
	QString text = in1->text(in2);//�����ѡ������
	QRegExp temp(QString::fromLocal8Bit("��"));
	if (text.contains(temp))//�鿴�Ƿ�ƥ�䵥����ʾ��ģ��
	{
		QRegularExpression exp("\\d{1,2}");
		QRegularExpressionMatch match = exp.match(text);
		if (match.hasMatch())
		{
			QString matched = match.captured();
			pdoc->m_CellIdx = matched.toInt();
		}
	}
	if (in1->parent())
	{
		//���������ǿ��Ƶ��е�xyz����һ���������ֵ
		if (in1->parent()->parent())
		{
			//�޸�xyz�������ֵ
			ui.treeWidget->openPersistentEditor(in1, 1);
		}
	}
}

void QtGuiApplication1::slotChanged(QTreeWidgetItem * in1, int in2)
{
	if (option.isLcokOn())
	{
		//��ȡ�ı��ֵ
		QString a = in1->text(in2);
		ui.treeWidget->closePersistentEditor(in1, in2);
		MyDoc::Ptr pdoc = MyDoc::getInstance();
		pdoc->clearScrn();
		option.m_again = 1;
		//���Ƶ�ı�Ϊ�ղ��޸ĵ�ֵ
		//------------
		//�ҵ�үү�ڵ㣬Ҳ���ǵڼ�����
		QString matched;
		QRegularExpression exp("\\d{1,2}");
		QRegularExpressionMatch match = exp.match(in1->parent()->parent()->parent()->text(0));
		if (match.hasMatch())
		{
			matched = match.captured();
		}
		//�ҵ��ǵڼ����ڵ�
		QString temp = in1->parent()->text(0);
		//�ҵ�����λ��
		if (in1->text(0) == "x")
		{
			pdoc->m_ShowData4D[matched.toInt()][temp.toInt()].x = a.toFloat();
		}
		else		if (in1->text(0) == "y")
		{
			pdoc->m_ShowData4D[matched.toInt()][temp.toInt()].y = a.toFloat();
		}
		else		if (in1->text(0) == "z")
		{
			pdoc->m_ShowData4D[matched.toInt()][temp.toInt()].z = a.toFloat();
		}
		//���»���ͼ��
		QOpenGLContext * openglContext = ui.centralWidget->context();//�����Ⱦ����
		option.OnBnClickedCreateModel(openglContext, this);
		option.createCP(openglContext);//������ʾ����
	}
}

