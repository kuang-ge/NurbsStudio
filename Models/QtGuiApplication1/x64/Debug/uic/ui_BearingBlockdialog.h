/********************************************************************************
** Form generated from reading UI file 'BearingBlockdialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_BEARINGBLOCKDIALOG_H
#define UI_BEARINGBLOCKDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>

QT_BEGIN_NAMESPACE

class Ui_BearingBlockdialog
{
public:
    QCheckBox *checkBox;
    QCheckBox *checkBox_2;
    QPushButton *pushButton;
    QLineEdit *lineEdit;
    QLabel *label;
    QLineEdit *lineEdit_2;

    void setupUi(QDialog *BearingBlockdialog)
    {
        if (BearingBlockdialog->objectName().isEmpty())
            BearingBlockdialog->setObjectName(QStringLiteral("BearingBlockdialog"));
        BearingBlockdialog->resize(466, 423);
        checkBox = new QCheckBox(BearingBlockdialog);
        checkBox->setObjectName(QStringLiteral("checkBox"));
        checkBox->setGeometry(QRect(50, 290, 71, 21));
        checkBox_2 = new QCheckBox(BearingBlockdialog);
        checkBox_2->setObjectName(QStringLiteral("checkBox_2"));
        checkBox_2->setGeometry(QRect(160, 330, 71, 21));
        pushButton = new QPushButton(BearingBlockdialog);
        pushButton->setObjectName(QStringLiteral("pushButton"));
        pushButton->setGeometry(QRect(160, 360, 71, 31));
        lineEdit = new QLineEdit(BearingBlockdialog);
        lineEdit->setObjectName(QStringLiteral("lineEdit"));
        lineEdit->setGeometry(QRect(260, 290, 161, 21));
        label = new QLabel(BearingBlockdialog);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(230, 290, 51, 21));
        lineEdit_2 = new QLineEdit(BearingBlockdialog);
        lineEdit_2->setObjectName(QStringLiteral("lineEdit_2"));
        lineEdit_2->setGeometry(QRect(60, 40, 113, 20));

        retranslateUi(BearingBlockdialog);

        QMetaObject::connectSlotsByName(BearingBlockdialog);
    } // setupUi

    void retranslateUi(QDialog *BearingBlockdialog)
    {
        BearingBlockdialog->setWindowTitle(QApplication::translate("BearingBlockdialog", "BearingBlockdialog", Q_NULLPTR));
        checkBox->setText(QApplication::translate("BearingBlockdialog", "\350\276\223\345\207\272\346\250\241\345\236\213", Q_NULLPTR));
        checkBox_2->setText(QApplication::translate("BearingBlockdialog", "\351\273\230\350\256\244\345\217\202\346\225\260", Q_NULLPTR));
        pushButton->setText(QApplication::translate("BearingBlockdialog", "PushButton", Q_NULLPTR));
        label->setText(QApplication::translate("BearingBlockdialog", "\350\267\257\345\276\204", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class BearingBlockdialog: public Ui_BearingBlockdialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_BEARINGBLOCKDIALOG_H
