#include "FeatureNtkDialog.h"

FeatureNtkDialog::FeatureNtkDialog(QWidget *parent)
	: QDialog(parent)
{
	ui.setupUi(this);
	connect(ui.button_select, &QPushButton::clicked, [=]() {
		QString text = ui.select_shapes->currentText();

		QTextCursor textCursor = ui.text_input->textCursor();//�õ���ǰ�Ĺ��

		textCursor.movePosition(QTextCursor::End);

		if (textCursor.hasSelection())//�����ѡ�У���ȡ������������Ӱ��
			textCursor.clearSelection();
		textCursor.deletePreviousChar();//ɾ��ǰһ���ַ�
		int cur_mode = -1;
		ui.text_input->setTextCursor(textCursor);
		if (text == QString::fromLocal8Bit("������")) {
			cur_mode = 1;
			QString context;
			context = "RecTangle (double L0,double H0,double x,double y,double z,double alph)";
			ui.text_output->setText(context);
		}
		
		connect(ui.button_generate, &QPushButton::clicked, [=]() {
			mode = cur_mode;
			this->close();
		});



	});


}

FeatureNtkDialog::~FeatureNtkDialog()
{
}
