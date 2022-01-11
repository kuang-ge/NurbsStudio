#include "readFileDialog.h"
#include "qfiledialog.h"
#include "MyDoc.h"
#include "qmessagebox.h"

readFileDialog::readFileDialog(QWidget *parent)
	: QDialog(parent)
{
	ui.setupUi(this);

	connect(ui.cancle, &QPushButton::clicked, this, &QWidget::close);

	connect(ui.selectFile, &QPushButton::clicked, [=]() {
		QString filePath = QFileDialog::getOpenFileName(this, QString::fromLocal8Bit("ѡ���ļ�"), "", "*.txt");
		ui.filePath->setText(filePath);
	});

	connect(ui.conform, &QPushButton::clicked, [=]() {
		if (ui.filePath->text() == NULL)
		{
			QMessageBox box(QMessageBox::Warning, QString::fromLocal8Bit("����"), QString::fromLocal8Bit("·������Ϊ��"));
			box.setStandardButtons(QMessageBox::Ok);
			box.setButtonText(QMessageBox::Ok, QString::fromLocal8Bit("ȷ ��"));
			box.exec();
		}
		else {
			MyDoc::Ptr pdoc = MyDoc::getInstance();
			pdoc->modelType = ui.modelType->currentIndex();
			pdoc->path = ui.filePath->text();
			this->close();
		}
	});
}

readFileDialog::~readFileDialog()
{
}
