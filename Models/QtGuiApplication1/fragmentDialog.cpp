#include "fragmentDialog.h"
#include <qmessagebox.h>

#include "MyDoc.h"

fragmentDialog::fragmentDialog(QWidget *parent)
	: QDialog(parent)
{
	ui.setupUi(this);
	
	MyDoc::Ptr pdoc = MyDoc::getInstance();
	
	connect(ui.cancle, &QPushButton::clicked, this, &QWidget::close);

	ui.total->setValidator(new QIntValidator(0, 1000, this));//����ֻ����������
	ui.begin->setValidator(new QIntValidator(0, 1000, this));//����ֻ����������
	ui.end->setValidator(new QIntValidator(0, 1000, this));//����ֻ����������
	ui.a->setValidator(new QIntValidator(0, 1000, this));//����ֻ����������
	ui.b->setValidator(new QIntValidator(0, 1000, this));//����ֻ����������
	ui.c->setValidator(new QIntValidator(0, 1000, this));//����ֻ����������

	connect(ui.confirm, &QPushButton::clicked, [=]() {
		if (ui.total->text() == NULL|| ui.begin->text()==NULL|| ui.a->text()==NULL|| ui.b->text()==NULL|| ui.c->text()==NULL)
		{
			QMessageBox box(QMessageBox::Warning, QString::fromLocal8Bit("����"), QString::fromLocal8Bit("���ݲ���Ϊ��"));
			box.setStandardButtons(QMessageBox::Ok);
			box.setButtonText(QMessageBox::Ok, QString::fromLocal8Bit("ȷ ��"));
			box.exec();
		}
		else
		{
			pdoc->m_fragmentDegree = ui.total->text().toInt();
			pdoc->m_fragmentBegin = ui.begin->text().toInt();
			pdoc->m_fragmentEnd = ui.end->text().toInt();
			pdoc->m_fragmentA = ui.a->text().toInt();
			pdoc->m_fragmentB = ui.b->text().toInt();
			pdoc->m_fragmentC = ui.c->text().toInt();
			pdoc->m_Fragment = true;
			this->close();
		}
	});
}

fragmentDialog::~fragmentDialog()
{
}
