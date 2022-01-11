/********************************************************************************
** Form generated from reading UI file 'fragmentDialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_FRAGMENTDIALOG_H
#define UI_FRAGMENTDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_fragmentDialog
{
public:
    QVBoxLayout *verticalLayout;
    QWidget *widget;
    QHBoxLayout *horizontalLayout;
    QSpacerItem *horizontalSpacer_3;
    QLabel *label;
    QSpacerItem *horizontalSpacer_17;
    QLineEdit *total;
    QSpacerItem *horizontalSpacer_4;
    QWidget *widget_2;
    QHBoxLayout *horizontalLayout_2;
    QSpacerItem *horizontalSpacer_5;
    QLineEdit *begin;
    QLabel *label_2;
    QSpacerItem *horizontalSpacer_6;
    QWidget *widget_3;
    QHBoxLayout *horizontalLayout_3;
    QSpacerItem *horizontalSpacer_7;
    QLineEdit *end;
    QLabel *label_3;
    QSpacerItem *horizontalSpacer_8;
    QWidget *widget_4;
    QHBoxLayout *horizontalLayout_5;
    QSpacerItem *horizontalSpacer_9;
    QLabel *label_4;
    QLineEdit *a;
    QSpacerItem *horizontalSpacer_10;
    QWidget *widget_5;
    QHBoxLayout *horizontalLayout_7;
    QSpacerItem *horizontalSpacer_12;
    QLabel *label_5;
    QLineEdit *b;
    QSpacerItem *horizontalSpacer_11;
    QWidget *widget_6;
    QHBoxLayout *horizontalLayout_4;
    QSpacerItem *horizontalSpacer_13;
    QLabel *label_6;
    QLineEdit *c;
    QSpacerItem *horizontalSpacer_15;
    QWidget *widget_7;
    QHBoxLayout *horizontalLayout_6;
    QSpacerItem *horizontalSpacer_14;
    QLabel *label_7;
    QLineEdit *d;
    QSpacerItem *horizontalSpacer_16;
    QWidget *widget_8;
    QHBoxLayout *horizontalLayout_8;
    QSpacerItem *horizontalSpacer;
    QPushButton *confirm;
    QPushButton *cancle;
    QSpacerItem *horizontalSpacer_2;

    void setupUi(QWidget *fragmentDialog)
    {
        if (fragmentDialog->objectName().isEmpty())
            fragmentDialog->setObjectName(QStringLiteral("fragmentDialog"));
        fragmentDialog->resize(400, 422);
        fragmentDialog->setMinimumSize(QSize(400, 400));
        fragmentDialog->setMaximumSize(QSize(400, 422));
        verticalLayout = new QVBoxLayout(fragmentDialog);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        widget = new QWidget(fragmentDialog);
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

        horizontalSpacer_17 = new QSpacerItem(10, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_17);

        total = new QLineEdit(widget);
        total->setObjectName(QStringLiteral("total"));

        horizontalLayout->addWidget(total);

        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_4);


        verticalLayout->addWidget(widget);

        widget_2 = new QWidget(fragmentDialog);
        widget_2->setObjectName(QStringLiteral("widget_2"));
        horizontalLayout_2 = new QHBoxLayout(widget_2);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        horizontalSpacer_5 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_5);

        begin = new QLineEdit(widget_2);
        begin->setObjectName(QStringLiteral("begin"));

        horizontalLayout_2->addWidget(begin);

        label_2 = new QLabel(widget_2);
        label_2->setObjectName(QStringLiteral("label_2"));

        horizontalLayout_2->addWidget(label_2);

        horizontalSpacer_6 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_6);


        verticalLayout->addWidget(widget_2);

        widget_3 = new QWidget(fragmentDialog);
        widget_3->setObjectName(QStringLiteral("widget_3"));
        horizontalLayout_3 = new QHBoxLayout(widget_3);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalSpacer_7 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_7);

        end = new QLineEdit(widget_3);
        end->setObjectName(QStringLiteral("end"));

        horizontalLayout_3->addWidget(end);

        label_3 = new QLabel(widget_3);
        label_3->setObjectName(QStringLiteral("label_3"));

        horizontalLayout_3->addWidget(label_3);

        horizontalSpacer_8 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_8);


        verticalLayout->addWidget(widget_3);

        widget_4 = new QWidget(fragmentDialog);
        widget_4->setObjectName(QStringLiteral("widget_4"));
        horizontalLayout_5 = new QHBoxLayout(widget_4);
        horizontalLayout_5->setSpacing(6);
        horizontalLayout_5->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        horizontalSpacer_9 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_5->addItem(horizontalSpacer_9);

        label_4 = new QLabel(widget_4);
        label_4->setObjectName(QStringLiteral("label_4"));

        horizontalLayout_5->addWidget(label_4);

        a = new QLineEdit(widget_4);
        a->setObjectName(QStringLiteral("a"));

        horizontalLayout_5->addWidget(a);

        horizontalSpacer_10 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_5->addItem(horizontalSpacer_10);


        verticalLayout->addWidget(widget_4);

        widget_5 = new QWidget(fragmentDialog);
        widget_5->setObjectName(QStringLiteral("widget_5"));
        horizontalLayout_7 = new QHBoxLayout(widget_5);
        horizontalLayout_7->setSpacing(6);
        horizontalLayout_7->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        horizontalSpacer_12 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_7->addItem(horizontalSpacer_12);

        label_5 = new QLabel(widget_5);
        label_5->setObjectName(QStringLiteral("label_5"));

        horizontalLayout_7->addWidget(label_5);

        b = new QLineEdit(widget_5);
        b->setObjectName(QStringLiteral("b"));

        horizontalLayout_7->addWidget(b);

        horizontalSpacer_11 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_7->addItem(horizontalSpacer_11);


        verticalLayout->addWidget(widget_5);

        widget_6 = new QWidget(fragmentDialog);
        widget_6->setObjectName(QStringLiteral("widget_6"));
        horizontalLayout_4 = new QHBoxLayout(widget_6);
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        horizontalSpacer_13 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_13);

        label_6 = new QLabel(widget_6);
        label_6->setObjectName(QStringLiteral("label_6"));

        horizontalLayout_4->addWidget(label_6);

        c = new QLineEdit(widget_6);
        c->setObjectName(QStringLiteral("c"));

        horizontalLayout_4->addWidget(c);

        horizontalSpacer_15 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_15);


        verticalLayout->addWidget(widget_6);

        widget_7 = new QWidget(fragmentDialog);
        widget_7->setObjectName(QStringLiteral("widget_7"));
        horizontalLayout_6 = new QHBoxLayout(widget_7);
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        horizontalSpacer_14 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_6->addItem(horizontalSpacer_14);

        label_7 = new QLabel(widget_7);
        label_7->setObjectName(QStringLiteral("label_7"));

        horizontalLayout_6->addWidget(label_7);

        d = new QLineEdit(widget_7);
        d->setObjectName(QStringLiteral("d"));

        horizontalLayout_6->addWidget(d);

        horizontalSpacer_16 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_6->addItem(horizontalSpacer_16);


        verticalLayout->addWidget(widget_7);

        widget_8 = new QWidget(fragmentDialog);
        widget_8->setObjectName(QStringLiteral("widget_8"));
        horizontalLayout_8 = new QHBoxLayout(widget_8);
        horizontalLayout_8->setSpacing(6);
        horizontalLayout_8->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));
        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_8->addItem(horizontalSpacer);

        confirm = new QPushButton(widget_8);
        confirm->setObjectName(QStringLiteral("confirm"));

        horizontalLayout_8->addWidget(confirm);

        cancle = new QPushButton(widget_8);
        cancle->setObjectName(QStringLiteral("cancle"));

        horizontalLayout_8->addWidget(cancle);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_8->addItem(horizontalSpacer_2);


        verticalLayout->addWidget(widget_8);


        retranslateUi(fragmentDialog);

        QMetaObject::connectSlotsByName(fragmentDialog);
    } // setupUi

    void retranslateUi(QWidget *fragmentDialog)
    {
        fragmentDialog->setWindowTitle(QApplication::translate("fragmentDialog", "fragmentDialog", Q_NULLPTR));
        label->setText(QApplication::translate("fragmentDialog", "\346\200\273\347\211\207\346\225\260", Q_NULLPTR));
        label_2->setText(QApplication::translate("fragmentDialog", "\350\265\267\345\247\213\344\275\215\347\275\256", Q_NULLPTR));
        label_3->setText(QApplication::translate("fragmentDialog", "\347\273\210\347\202\271\344\275\215\347\275\256", Q_NULLPTR));
        label_4->setText(QApplication::translate("fragmentDialog", "A", Q_NULLPTR));
        label_5->setText(QApplication::translate("fragmentDialog", "B", Q_NULLPTR));
        label_6->setText(QApplication::translate("fragmentDialog", "C", Q_NULLPTR));
        label_7->setText(QApplication::translate("fragmentDialog", "D", Q_NULLPTR));
        confirm->setText(QApplication::translate("fragmentDialog", "\347\241\256\345\256\232", Q_NULLPTR));
        cancle->setText(QApplication::translate("fragmentDialog", "\345\217\226\346\266\210", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class fragmentDialog: public Ui_fragmentDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_FRAGMENTDIALOG_H
