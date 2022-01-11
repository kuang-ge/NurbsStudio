#include "camerLocatonDialog.h"
#include "MyDoc.h"
#include <qmessagebox.h>
#include <QKeyEvent>

camerLocatonDialog::camerLocatonDialog(QWidget *parent)
	: QDialog(parent)
{
	ui.setupUi(this);
	MyDoc::Ptr pdoc = MyDoc::getInstance();

	ui.X->setText(QString::number(pdoc->cameraX));//��ȡ��ǰ����Ϣ
	ui.Y->setText(QString::number(pdoc->cameraY));
	ui.Z->setText(QString::number(pdoc->cameraZ));

	ui.X->setValidator(new QIntValidator(0, 10000, this));//����ֻ���������֣��ҷ�Χ��0-10000
	ui.Y->setValidator(new QIntValidator(0, 10000, this));
	ui.Z->setValidator(new QIntValidator(0, 10000, this));

	connect(ui.conform, &QPushButton::clicked, [=]() {
		if (ui.X->text() == NULL || ui.X->text() == NULL || ui.X->text() == NULL)
		{
			QMessageBox box(QMessageBox::Warning, QString::fromLocal8Bit("����"), QString::fromLocal8Bit("���ݲ���Ϊ��"));
			box.setStandardButtons(QMessageBox::Ok);
			box.setButtonText(QMessageBox::Ok, QString::fromLocal8Bit("ȷ ��"));
			box.exec();
		}
		else
		{
			pdoc->cameraX = ui.X->text().toFloat();
			pdoc->cameraY = ui.Y->text().toFloat();
			pdoc->cameraZ = ui.Z->text().toFloat();
			pdoc->m_restView = true;
			this->close();
		}
	});

	connect(ui.cancle, &QPushButton::clicked, [=]() {
		this->close();
	});
}

camerLocatonDialog::~camerLocatonDialog()
{
}
