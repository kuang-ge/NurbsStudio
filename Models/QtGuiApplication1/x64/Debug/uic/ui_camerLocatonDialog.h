/********************************************************************************
** Form generated from reading UI file 'camerLocatonDialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CAMERLOCATONDIALOG_H
#define UI_CAMERLOCATONDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_camerLocatonDialog
{
public:
    QVBoxLayout *verticalLayout;
    QWidget *widget;
    QHBoxLayout *horizontalLayout;
    QSpacerItem *horizontalSpacer_3;
    QLabel *label;
    QLineEdit *X;
    QSpacerItem *horizontalSpacer_4;
    QWidget *widget_2;
    QHBoxLayout *horizontalLayout_2;
    QSpacerItem *horizontalSpacer_5;
    QLabel *label_2;
    QLineEdit *Y;
    QSpacerItem *horizontalSpacer_6;
    QWidget *widget_3;
    QHBoxLayout *horizontalLayout_3;
    QSpacerItem *horizontalSpacer_8;
    QLabel *label_3;
    QLineEdit *Z;
    QSpacerItem *horizontalSpacer_7;
    QWidget *widget_4;
    QHBoxLayout *horizontalLayout_4;
    QSpacerItem *horizontalSpacer;
    QPushButton *conform;
    QPushButton *cancle;
    QSpacerItem *horizontalSpacer_2;

    void setupUi(QDialog *camerLocatonDialog)
    {
        if (camerLocatonDialog->objectName().isEmpty())
            camerLocatonDialog->setObjectName(QStringLiteral("camerLocatonDialog"));
        camerLocatonDialog->resize(250, 250);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(camerLocatonDialog->sizePolicy().hasHeightForWidth());
        camerLocatonDialog->setSizePolicy(sizePolicy);
        camerLocatonDialog->setMinimumSize(QSize(250, 250));
        camerLocatonDialog->setMaximumSize(QSize(250, 250));
        verticalLayout = new QVBoxLayout(camerLocatonDialog);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        widget = new QWidget(camerLocatonDialog);
        widget->setObjectName(QStringLiteral("widget"));
        horizontalLayout = new QHBoxLayout(widget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalSpacer_3 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_3);

        label = new QLabel(widget);
        label->setObjectName(QStringLiteral("label"));

        horizontalLayout->addWidget(label);

        X = new QLineEdit(widget);
        X->setObjectName(QStringLiteral("X"));
        X->setClearButtonEnabled(false);

        horizontalLayout->addWidget(X);

        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_4);


        verticalLayout->addWidget(widget);

        widget_2 = new QWidget(camerLocatonDialog);
        widget_2->setObjectName(QStringLiteral("widget_2"));
        horizontalLayout_2 = new QHBoxLayout(widget_2);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        horizontalSpacer_5 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_5);

        label_2 = new QLabel(widget_2);
        label_2->setObjectName(QStringLiteral("label_2"));

        horizontalLayout_2->addWidget(label_2);

        Y = new QLineEdit(widget_2);
        Y->setObjectName(QStringLiteral("Y"));

        horizontalLayout_2->addWidget(Y);

        horizontalSpacer_6 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_6);


        verticalLayout->addWidget(widget_2);

        widget_3 = new QWidget(camerLocatonDialog);
        widget_3->setObjectName(QStringLiteral("widget_3"));
        horizontalLayout_3 = new QHBoxLayout(widget_3);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalSpacer_8 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_8);

        label_3 = new QLabel(widget_3);
        label_3->setObjectName(QStringLiteral("label_3"));

        horizontalLayout_3->addWidget(label_3);

        Z = new QLineEdit(widget_3);
        Z->setObjectName(QStringLiteral("Z"));

        horizontalLayout_3->addWidget(Z);

        horizontalSpacer_7 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_7);


        verticalLayout->addWidget(widget_3);

        widget_4 = new QWidget(camerLocatonDialog);
        widget_4->setObjectName(QStringLiteral("widget_4"));
        horizontalLayout_4 = new QHBoxLayout(widget_4);
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer);

        conform = new QPushButton(widget_4);
        conform->setObjectName(QStringLiteral("conform"));

        horizontalLayout_4->addWidget(conform);

        cancle = new QPushButton(widget_4);
        cancle->setObjectName(QStringLiteral("cancle"));

        horizontalLayout_4->addWidget(cancle);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_2);


        verticalLayout->addWidget(widget_4);


        retranslateUi(camerLocatonDialog);

        QMetaObject::connectSlotsByName(camerLocatonDialog);
    } // setupUi

    void retranslateUi(QDialog *camerLocatonDialog)
    {
        camerLocatonDialog->setWindowTitle(QApplication::translate("camerLocatonDialog", "camerLocatonDialog", Q_NULLPTR));
        label->setText(QApplication::translate("camerLocatonDialog", "X\350\275\264\344\275\215\347\275\256", Q_NULLPTR));
        label_2->setText(QApplication::translate("camerLocatonDialog", "Y\350\275\264\344\275\215\347\275\256", Q_NULLPTR));
        label_3->setText(QApplication::translate("camerLocatonDialog", "Z\350\275\264\344\275\215\347\275\256", Q_NULLPTR));
        conform->setText(QApplication::translate("camerLocatonDialog", "\347\241\256\345\256\232", Q_NULLPTR));
        cancle->setText(QApplication::translate("camerLocatonDialog", "\345\217\226\346\266\210", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class camerLocatonDialog: public Ui_camerLocatonDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CAMERLOCATONDIALOG_H
