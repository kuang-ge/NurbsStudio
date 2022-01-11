/********************************************************************************
** Form generated from reading UI file 'ChangeFineScaleDialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CHANGEFINESCALEDIALOG_H
#define UI_CHANGEFINESCALEDIALOG_H

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

class Ui_ChangeFineScaleDialog
{
public:
    QVBoxLayout *verticalLayout;
    QWidget *widget;
    QHBoxLayout *horizontalLayout;
    QSpacerItem *horizontalSpacer_6;
    QLabel *label;
    QLineEdit *U_NUM;
    QSpacerItem *horizontalSpacer_7;
    QWidget *widget_2;
    QHBoxLayout *horizontalLayout_2;
    QSpacerItem *horizontalSpacer_8;
    QLabel *label_2;
    QLineEdit *V_NUM;
    QSpacerItem *horizontalSpacer_9;
    QWidget *widget_3;
    QHBoxLayout *horizontalLayout_3;
    QSpacerItem *horizontalSpacer_4;
    QLabel *label_3;
    QLineEdit *W_NUM;
    QSpacerItem *horizontalSpacer_5;
    QWidget *widget_4;
    QHBoxLayout *horizontalLayout_4;
    QSpacerItem *horizontalSpacer;
    QPushButton *conform;
    QSpacerItem *horizontalSpacer_3;
    QPushButton *cancle;
    QSpacerItem *horizontalSpacer_2;

    void setupUi(QDialog *ChangeFineScaleDialog)
    {
        if (ChangeFineScaleDialog->objectName().isEmpty())
            ChangeFineScaleDialog->setObjectName(QStringLiteral("ChangeFineScaleDialog"));
        ChangeFineScaleDialog->resize(250, 250);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(ChangeFineScaleDialog->sizePolicy().hasHeightForWidth());
        ChangeFineScaleDialog->setSizePolicy(sizePolicy);
        ChangeFineScaleDialog->setMinimumSize(QSize(250, 250));
        ChangeFineScaleDialog->setMaximumSize(QSize(250, 250));
        verticalLayout = new QVBoxLayout(ChangeFineScaleDialog);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        widget = new QWidget(ChangeFineScaleDialog);
        widget->setObjectName(QStringLiteral("widget"));
        horizontalLayout = new QHBoxLayout(widget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalSpacer_6 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_6);

        label = new QLabel(widget);
        label->setObjectName(QStringLiteral("label"));

        horizontalLayout->addWidget(label);

        U_NUM = new QLineEdit(widget);
        U_NUM->setObjectName(QStringLiteral("U_NUM"));

        horizontalLayout->addWidget(U_NUM);

        horizontalSpacer_7 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_7);


        verticalLayout->addWidget(widget);

        widget_2 = new QWidget(ChangeFineScaleDialog);
        widget_2->setObjectName(QStringLiteral("widget_2"));
        horizontalLayout_2 = new QHBoxLayout(widget_2);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        horizontalSpacer_8 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_8);

        label_2 = new QLabel(widget_2);
        label_2->setObjectName(QStringLiteral("label_2"));

        horizontalLayout_2->addWidget(label_2);

        V_NUM = new QLineEdit(widget_2);
        V_NUM->setObjectName(QStringLiteral("V_NUM"));

        horizontalLayout_2->addWidget(V_NUM);

        horizontalSpacer_9 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_9);


        verticalLayout->addWidget(widget_2);

        widget_3 = new QWidget(ChangeFineScaleDialog);
        widget_3->setObjectName(QStringLiteral("widget_3"));
        horizontalLayout_3 = new QHBoxLayout(widget_3);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_4);

        label_3 = new QLabel(widget_3);
        label_3->setObjectName(QStringLiteral("label_3"));

        horizontalLayout_3->addWidget(label_3);

        W_NUM = new QLineEdit(widget_3);
        W_NUM->setObjectName(QStringLiteral("W_NUM"));

        horizontalLayout_3->addWidget(W_NUM);

        horizontalSpacer_5 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_5);


        verticalLayout->addWidget(widget_3);

        widget_4 = new QWidget(ChangeFineScaleDialog);
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

        horizontalSpacer_3 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_3);

        cancle = new QPushButton(widget_4);
        cancle->setObjectName(QStringLiteral("cancle"));

        horizontalLayout_4->addWidget(cancle);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_2);


        verticalLayout->addWidget(widget_4);


        retranslateUi(ChangeFineScaleDialog);

        QMetaObject::connectSlotsByName(ChangeFineScaleDialog);
    } // setupUi

    void retranslateUi(QDialog *ChangeFineScaleDialog)
    {
        ChangeFineScaleDialog->setWindowTitle(QApplication::translate("ChangeFineScaleDialog", "ChangeFineScaleDialog", Q_NULLPTR));
        label->setText(QApplication::translate("ChangeFineScaleDialog", "U\346\226\271\345\220\221\347\273\206\345\210\206\345\272\246", Q_NULLPTR));
        label_2->setText(QApplication::translate("ChangeFineScaleDialog", "V\346\226\271\345\220\221\347\273\206\345\210\206\345\272\246", Q_NULLPTR));
        label_3->setText(QApplication::translate("ChangeFineScaleDialog", "W\346\226\271\345\220\221\347\273\206\345\210\206\345\272\246", Q_NULLPTR));
        conform->setText(QApplication::translate("ChangeFineScaleDialog", "\347\241\256\345\256\232", Q_NULLPTR));
        cancle->setText(QApplication::translate("ChangeFineScaleDialog", "\345\217\226\346\266\210", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class ChangeFineScaleDialog: public Ui_ChangeFineScaleDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CHANGEFINESCALEDIALOG_H
