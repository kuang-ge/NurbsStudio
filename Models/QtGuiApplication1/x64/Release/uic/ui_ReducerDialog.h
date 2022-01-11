/********************************************************************************
** Form generated from reading UI file 'ReducerDialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_REDUCERDIALOG_H
#define UI_REDUCERDIALOG_H

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

class Ui_ReducerDialog
{
public:
    QLineEdit *lineEdit_gear_r1;
    QLineEdit *lineEdit_gear_r2;
    QLineEdit *lineEdit_shaft_r1;
    QLineEdit *lineEdit_shaft_r2;
    QLineEdit *lineEdit_shaft_r3;
    QLineEdit *lineEdit_shaft_l1;
    QLineEdit *lineEdit_shaft_l2;
    QLineEdit *lineEdit_shaft_h;
    QLineEdit *lineEdit_shaft_t;
    QLineEdit *lineEdit_path;
    QCheckBox *outputmodel;
    QCheckBox *zeroparams;
    QLabel *label;
    QPushButton *button_ok;
    QLabel *label_2;
    QLabel *label_3;
    QLabel *label_4;
    QLabel *label_5;
    QLabel *label_6;
    QLabel *label_7;
    QLabel *label_8;
    QLabel *label_9;
    QLineEdit *lineEdit_box_l;
    QLineEdit *lineEdit_box_t;
    QLabel *label_10;
    QLabel *label_11;
    QLabel *label_12;

    void setupUi(QDialog *ReducerDialog)
    {
        if (ReducerDialog->objectName().isEmpty())
            ReducerDialog->setObjectName(QStringLiteral("ReducerDialog"));
        ReducerDialog->resize(463, 513);
        lineEdit_gear_r1 = new QLineEdit(ReducerDialog);
        lineEdit_gear_r1->setObjectName(QStringLiteral("lineEdit_gear_r1"));
        lineEdit_gear_r1->setGeometry(QRect(50, 40, 113, 20));
        lineEdit_gear_r2 = new QLineEdit(ReducerDialog);
        lineEdit_gear_r2->setObjectName(QStringLiteral("lineEdit_gear_r2"));
        lineEdit_gear_r2->setGeometry(QRect(240, 40, 113, 20));
        lineEdit_shaft_r1 = new QLineEdit(ReducerDialog);
        lineEdit_shaft_r1->setObjectName(QStringLiteral("lineEdit_shaft_r1"));
        lineEdit_shaft_r1->setGeometry(QRect(50, 90, 113, 20));
        lineEdit_shaft_r2 = new QLineEdit(ReducerDialog);
        lineEdit_shaft_r2->setObjectName(QStringLiteral("lineEdit_shaft_r2"));
        lineEdit_shaft_r2->setGeometry(QRect(240, 90, 113, 20));
        lineEdit_shaft_r3 = new QLineEdit(ReducerDialog);
        lineEdit_shaft_r3->setObjectName(QStringLiteral("lineEdit_shaft_r3"));
        lineEdit_shaft_r3->setGeometry(QRect(50, 140, 113, 20));
        lineEdit_shaft_l1 = new QLineEdit(ReducerDialog);
        lineEdit_shaft_l1->setObjectName(QStringLiteral("lineEdit_shaft_l1"));
        lineEdit_shaft_l1->setGeometry(QRect(240, 140, 113, 20));
        lineEdit_shaft_l2 = new QLineEdit(ReducerDialog);
        lineEdit_shaft_l2->setObjectName(QStringLiteral("lineEdit_shaft_l2"));
        lineEdit_shaft_l2->setGeometry(QRect(50, 190, 113, 20));
        lineEdit_shaft_h = new QLineEdit(ReducerDialog);
        lineEdit_shaft_h->setObjectName(QStringLiteral("lineEdit_shaft_h"));
        lineEdit_shaft_h->setGeometry(QRect(240, 190, 113, 20));
        lineEdit_shaft_t = new QLineEdit(ReducerDialog);
        lineEdit_shaft_t->setObjectName(QStringLiteral("lineEdit_shaft_t"));
        lineEdit_shaft_t->setGeometry(QRect(50, 240, 113, 20));
        lineEdit_path = new QLineEdit(ReducerDialog);
        lineEdit_path->setObjectName(QStringLiteral("lineEdit_path"));
        lineEdit_path->setGeometry(QRect(220, 340, 171, 20));
        outputmodel = new QCheckBox(ReducerDialog);
        outputmodel->setObjectName(QStringLiteral("outputmodel"));
        outputmodel->setGeometry(QRect(80, 340, 81, 21));
        zeroparams = new QCheckBox(ReducerDialog);
        zeroparams->setObjectName(QStringLiteral("zeroparams"));
        zeroparams->setGeometry(QRect(150, 390, 71, 21));
        label = new QLabel(ReducerDialog);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(180, 340, 51, 21));
        button_ok = new QPushButton(ReducerDialog);
        button_ok->setObjectName(QStringLiteral("button_ok"));
        button_ok->setGeometry(QRect(150, 420, 81, 31));
        label_2 = new QLabel(ReducerDialog);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(70, 10, 61, 21));
        label_3 = new QLabel(ReducerDialog);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(260, 10, 61, 21));
        label_4 = new QLabel(ReducerDialog);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setGeometry(QRect(70, 70, 91, 21));
        label_5 = new QLabel(ReducerDialog);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setGeometry(QRect(260, 70, 81, 21));
        label_6 = new QLabel(ReducerDialog);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setGeometry(QRect(70, 120, 91, 21));
        label_7 = new QLabel(ReducerDialog);
        label_7->setObjectName(QStringLiteral("label_7"));
        label_7->setGeometry(QRect(260, 120, 91, 21));
        label_8 = new QLabel(ReducerDialog);
        label_8->setObjectName(QStringLiteral("label_8"));
        label_8->setGeometry(QRect(70, 170, 91, 21));
        label_9 = new QLabel(ReducerDialog);
        label_9->setObjectName(QStringLiteral("label_9"));
        label_9->setGeometry(QRect(260, 170, 91, 21));
        lineEdit_box_l = new QLineEdit(ReducerDialog);
        lineEdit_box_l->setObjectName(QStringLiteral("lineEdit_box_l"));
        lineEdit_box_l->setGeometry(QRect(240, 240, 113, 20));
        lineEdit_box_t = new QLineEdit(ReducerDialog);
        lineEdit_box_t->setObjectName(QStringLiteral("lineEdit_box_t"));
        lineEdit_box_t->setGeometry(QRect(50, 290, 113, 20));
        label_10 = new QLabel(ReducerDialog);
        label_10->setObjectName(QStringLiteral("label_10"));
        label_10->setGeometry(QRect(80, 220, 91, 21));
        label_11 = new QLabel(ReducerDialog);
        label_11->setObjectName(QStringLiteral("label_11"));
        label_11->setGeometry(QRect(270, 220, 91, 21));
        label_12 = new QLabel(ReducerDialog);
        label_12->setObjectName(QStringLiteral("label_12"));
        label_12->setGeometry(QRect(80, 270, 91, 21));

        retranslateUi(ReducerDialog);

        QMetaObject::connectSlotsByName(ReducerDialog);
    } // setupUi

    void retranslateUi(QDialog *ReducerDialog)
    {
        ReducerDialog->setWindowTitle(QApplication::translate("ReducerDialog", "ReducerDialog", Q_NULLPTR));
        outputmodel->setText(QApplication::translate("ReducerDialog", "\350\276\223\345\207\272\346\250\241\345\236\213", Q_NULLPTR));
        zeroparams->setText(QApplication::translate("ReducerDialog", "\351\273\230\350\256\244\345\217\202\346\225\260", Q_NULLPTR));
        label->setText(QApplication::translate("ReducerDialog", "\350\267\257\345\276\204", Q_NULLPTR));
        button_ok->setText(QApplication::translate("ReducerDialog", "\347\241\256\345\256\232", Q_NULLPTR));
        label_2->setText(QApplication::translate("ReducerDialog", "\345\244\247\351\275\277\350\275\256\345\215\212\345\276\204", Q_NULLPTR));
        label_3->setText(QApplication::translate("ReducerDialog", "\345\260\217\351\275\277\350\275\256\345\215\212\345\276\204", Q_NULLPTR));
        label_4->setText(QApplication::translate("ReducerDialog", "\345\244\247\351\275\277\350\275\256\350\275\264\345\215\212\345\276\204", Q_NULLPTR));
        label_5->setText(QApplication::translate("ReducerDialog", "\345\260\217\351\275\277\350\275\256\350\275\2641\345\215\212\345\276\204", Q_NULLPTR));
        label_6->setText(QApplication::translate("ReducerDialog", "\345\260\217\351\275\277\350\275\256\350\275\2642\345\215\212\345\276\204", Q_NULLPTR));
        label_7->setText(QApplication::translate("ReducerDialog", "\345\244\247\345\260\217\351\275\277\350\275\256\350\275\264\351\227\264\350\267\235", Q_NULLPTR));
        label_8->setText(QApplication::translate("ReducerDialog", "\345\260\217\351\275\277\350\275\256\350\275\264\351\227\264\350\267\235", Q_NULLPTR));
        label_9->setText(QApplication::translate("ReducerDialog", "\351\275\277\350\275\256\350\275\264\351\253\230\345\272\246", Q_NULLPTR));
        label_10->setText(QApplication::translate("ReducerDialog", "\351\275\277\350\275\256\347\255\222\345\216\232", Q_NULLPTR));
        label_11->setText(QApplication::translate("ReducerDialog", "\347\256\261\344\275\223\345\256\275", Q_NULLPTR));
        label_12->setText(QApplication::translate("ReducerDialog", "\347\256\261\344\275\223\345\216\232", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class ReducerDialog: public Ui_ReducerDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_REDUCERDIALOG_H
