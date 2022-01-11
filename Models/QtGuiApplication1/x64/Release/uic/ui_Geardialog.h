/********************************************************************************
** Form generated from reading UI file 'Geardialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GEARDIALOG_H
#define UI_GEARDIALOG_H

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

class Ui_Geardialog
{
public:
    QLineEdit *lineEdit_Dk;
    QLineEdit *lineEdit_B;
    QLineEdit *lineEdit_x;
    QLineEdit *lineEdit_cx;
    QLineEdit *lineEdit_alph;
    QLineEdit *lineEdit_m;
    QLineEdit *lineEdit_hax;
    QLineEdit *lineEdit_z;
    QLabel *label;
    QLabel *label_2;
    QLabel *label_3;
    QLabel *label_4;
    QLabel *label_5;
    QLabel *label_6;
    QLabel *label_7;
    QLabel *label_8;
    QPushButton *Button_ok;
    QLineEdit *lineEdit_path;
    QLabel *label_9;
    QCheckBox *zeroparams;
    QCheckBox *outputmodel;

    void setupUi(QDialog *Geardialog)
    {
        if (Geardialog->objectName().isEmpty())
            Geardialog->setObjectName(QStringLiteral("Geardialog"));
        Geardialog->resize(423, 424);
        lineEdit_Dk = new QLineEdit(Geardialog);
        lineEdit_Dk->setObjectName(QStringLiteral("lineEdit_Dk"));
        lineEdit_Dk->setGeometry(QRect(239, 235, 97, 20));
        lineEdit_Dk->setMaximumSize(QSize(100, 20));
        lineEdit_B = new QLineEdit(Geardialog);
        lineEdit_B->setObjectName(QStringLiteral("lineEdit_B"));
        lineEdit_B->setGeometry(QRect(239, 171, 97, 20));
        lineEdit_B->setMaximumSize(QSize(100, 20));
        lineEdit_x = new QLineEdit(Geardialog);
        lineEdit_x->setObjectName(QStringLiteral("lineEdit_x"));
        lineEdit_x->setGeometry(QRect(47, 235, 97, 20));
        lineEdit_x->setMaximumSize(QSize(100, 20));
        lineEdit_cx = new QLineEdit(Geardialog);
        lineEdit_cx->setObjectName(QStringLiteral("lineEdit_cx"));
        lineEdit_cx->setGeometry(QRect(47, 171, 97, 20));
        lineEdit_cx->setMaximumSize(QSize(100, 20));
        lineEdit_alph = new QLineEdit(Geardialog);
        lineEdit_alph->setObjectName(QStringLiteral("lineEdit_alph"));
        lineEdit_alph->setGeometry(QRect(47, 108, 97, 20));
        lineEdit_alph->setMaximumSize(QSize(100, 20));
        lineEdit_m = new QLineEdit(Geardialog);
        lineEdit_m->setObjectName(QStringLiteral("lineEdit_m"));
        lineEdit_m->setGeometry(QRect(47, 44, 97, 20));
        lineEdit_m->setMaximumSize(QSize(100, 20));
        lineEdit_hax = new QLineEdit(Geardialog);
        lineEdit_hax->setObjectName(QStringLiteral("lineEdit_hax"));
        lineEdit_hax->setGeometry(QRect(239, 108, 97, 20));
        lineEdit_hax->setMaximumSize(QSize(100, 20));
        lineEdit_z = new QLineEdit(Geardialog);
        lineEdit_z->setObjectName(QStringLiteral("lineEdit_z"));
        lineEdit_z->setGeometry(QRect(239, 44, 97, 20));
        lineEdit_z->setMaximumSize(QSize(100, 20));
        label = new QLabel(Geardialog);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(70, 20, 30, 16));
        label_2 = new QLabel(Geardialog);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(270, 22, 51, 20));
        label_3 = new QLabel(Geardialog);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(250, 85, 81, 21));
        label_4 = new QLabel(Geardialog);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setGeometry(QRect(260, 150, 51, 16));
        label_5 = new QLabel(Geardialog);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setGeometry(QRect(250, 215, 81, 21));
        label_6 = new QLabel(Geardialog);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setGeometry(QRect(60, 210, 71, 16));
        label_7 = new QLabel(Geardialog);
        label_7->setObjectName(QStringLiteral("label_7"));
        label_7->setGeometry(QRect(61, 150, 71, 16));
        label_8 = new QLabel(Geardialog);
        label_8->setObjectName(QStringLiteral("label_8"));
        label_8->setGeometry(QRect(70, 90, 45, 16));
        Button_ok = new QPushButton(Geardialog);
        Button_ok->setObjectName(QStringLiteral("Button_ok"));
        Button_ok->setGeometry(QRect(160, 360, 75, 23));
        lineEdit_path = new QLineEdit(Geardialog);
        lineEdit_path->setObjectName(QStringLiteral("lineEdit_path"));
        lineEdit_path->setGeometry(QRect(210, 300, 171, 20));
        label_9 = new QLabel(Geardialog);
        label_9->setObjectName(QStringLiteral("label_9"));
        label_9->setGeometry(QRect(180, 300, 51, 21));
        zeroparams = new QCheckBox(Geardialog);
        zeroparams->setObjectName(QStringLiteral("zeroparams"));
        zeroparams->setGeometry(QRect(140, 330, 101, 21));
        outputmodel = new QCheckBox(Geardialog);
        outputmodel->setObjectName(QStringLiteral("outputmodel"));
        outputmodel->setGeometry(QRect(60, 300, 91, 21));

        retranslateUi(Geardialog);

        QMetaObject::connectSlotsByName(Geardialog);
    } // setupUi

    void retranslateUi(QDialog *Geardialog)
    {
        Geardialog->setWindowTitle(QApplication::translate("Geardialog", "Geardialog", Q_NULLPTR));
        label->setText(QApplication::translate("Geardialog", "\346\250\241\346\225\260", Q_NULLPTR));
        label_2->setText(QApplication::translate("Geardialog", "\351\275\277\346\225\260", Q_NULLPTR));
        label_3->setText(QApplication::translate("Geardialog", "\351\275\277\351\241\266\351\253\230\347\263\273\346\225\260", Q_NULLPTR));
        label_4->setText(QApplication::translate("Geardialog", "\351\275\277\350\275\256\345\216\232", Q_NULLPTR));
        label_5->setText(QApplication::translate("Geardialog", "\351\275\277\350\275\256\350\275\264\347\233\264\345\276\204", Q_NULLPTR));
        label_6->setText(QApplication::translate("Geardialog", "\345\217\230\344\275\215\347\263\273\346\225\260", Q_NULLPTR));
        label_7->setText(QApplication::translate("Geardialog", "\351\241\266\351\232\231\347\263\273\346\225\260", Q_NULLPTR));
        label_8->setText(QApplication::translate("Geardialog", "\345\216\213\345\212\233\350\247\222", Q_NULLPTR));
        Button_ok->setText(QApplication::translate("Geardialog", "\347\241\256\345\256\232", Q_NULLPTR));
        label_9->setText(QApplication::translate("Geardialog", "\350\267\257\345\276\204", Q_NULLPTR));
        zeroparams->setText(QApplication::translate("Geardialog", "\351\273\230\350\256\244\345\217\202\346\225\260", Q_NULLPTR));
        outputmodel->setText(QApplication::translate("Geardialog", "\350\276\223\345\207\272\346\250\241\345\236\213", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Geardialog: public Ui_Geardialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GEARDIALOG_H
