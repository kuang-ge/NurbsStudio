/********************************************************************************
** Form generated from reading UI file 'singlePic.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SINGLEPIC_H
#define UI_SINGLEPIC_H

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

class Ui_singlePic
{
public:
    QVBoxLayout *verticalLayout;
    QWidget *widget;
    QHBoxLayout *horizontalLayout;
    QSpacerItem *horizontalSpacer_3;
    QLabel *totalPic;
    QSpacerItem *horizontalSpacer_7;
    QLineEdit *totalPicEdit;
    QSpacerItem *horizontalSpacer_4;
    QWidget *widget_4;
    QHBoxLayout *horizontalLayout_4;
    QSpacerItem *horizontalSpacer_8;
    QLabel *label_2;
    QLineEdit *lineEdit;
    QSpacerItem *horizontalSpacer_9;
    QWidget *widget_2;
    QHBoxLayout *horizontalLayout_2;
    QSpacerItem *horizontalSpacer_5;
    QLabel *label;
    QLineEdit *currentPicEdit;
    QSpacerItem *horizontalSpacer_6;
    QWidget *widget_3;
    QHBoxLayout *horizontalLayout_3;
    QSpacerItem *horizontalSpacer;
    QPushButton *lastPicButton;
    QPushButton *nextPicButton;
    QPushButton *closeButton;
    QSpacerItem *horizontalSpacer_2;

    void setupUi(QWidget *singlePic)
    {
        if (singlePic->objectName().isEmpty())
            singlePic->setObjectName(QStringLiteral("singlePic"));
        singlePic->resize(300, 250);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(singlePic->sizePolicy().hasHeightForWidth());
        singlePic->setSizePolicy(sizePolicy);
        singlePic->setMinimumSize(QSize(300, 250));
        singlePic->setMaximumSize(QSize(300, 250));
        verticalLayout = new QVBoxLayout(singlePic);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        widget = new QWidget(singlePic);
        widget->setObjectName(QStringLiteral("widget"));
        horizontalLayout = new QHBoxLayout(widget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalSpacer_3 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_3);

        totalPic = new QLabel(widget);
        totalPic->setObjectName(QStringLiteral("totalPic"));

        horizontalLayout->addWidget(totalPic);

        horizontalSpacer_7 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_7);

        totalPicEdit = new QLineEdit(widget);
        totalPicEdit->setObjectName(QStringLiteral("totalPicEdit"));
        totalPicEdit->setReadOnly(true);

        horizontalLayout->addWidget(totalPicEdit);

        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_4);


        verticalLayout->addWidget(widget);

        widget_4 = new QWidget(singlePic);
        widget_4->setObjectName(QStringLiteral("widget_4"));
        horizontalLayout_4 = new QHBoxLayout(widget_4);
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        horizontalSpacer_8 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_8);

        label_2 = new QLabel(widget_4);
        label_2->setObjectName(QStringLiteral("label_2"));

        horizontalLayout_4->addWidget(label_2);

        lineEdit = new QLineEdit(widget_4);
        lineEdit->setObjectName(QStringLiteral("lineEdit"));
        lineEdit->setReadOnly(false);

        horizontalLayout_4->addWidget(lineEdit);

        horizontalSpacer_9 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_9);


        verticalLayout->addWidget(widget_4);

        widget_2 = new QWidget(singlePic);
        widget_2->setObjectName(QStringLiteral("widget_2"));
        horizontalLayout_2 = new QHBoxLayout(widget_2);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        horizontalSpacer_5 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_5);

        label = new QLabel(widget_2);
        label->setObjectName(QStringLiteral("label"));

        horizontalLayout_2->addWidget(label);

        currentPicEdit = new QLineEdit(widget_2);
        currentPicEdit->setObjectName(QStringLiteral("currentPicEdit"));
        currentPicEdit->setReadOnly(false);

        horizontalLayout_2->addWidget(currentPicEdit);

        horizontalSpacer_6 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_6);


        verticalLayout->addWidget(widget_2);

        widget_3 = new QWidget(singlePic);
        widget_3->setObjectName(QStringLiteral("widget_3"));
        horizontalLayout_3 = new QHBoxLayout(widget_3);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer);

        lastPicButton = new QPushButton(widget_3);
        lastPicButton->setObjectName(QStringLiteral("lastPicButton"));

        horizontalLayout_3->addWidget(lastPicButton);

        nextPicButton = new QPushButton(widget_3);
        nextPicButton->setObjectName(QStringLiteral("nextPicButton"));

        horizontalLayout_3->addWidget(nextPicButton);

        closeButton = new QPushButton(widget_3);
        closeButton->setObjectName(QStringLiteral("closeButton"));

        horizontalLayout_3->addWidget(closeButton);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_2);


        verticalLayout->addWidget(widget_3);


        retranslateUi(singlePic);

        QMetaObject::connectSlotsByName(singlePic);
    } // setupUi

    void retranslateUi(QWidget *singlePic)
    {
        singlePic->setWindowTitle(QApplication::translate("singlePic", "singlePic", Q_NULLPTR));
        totalPic->setText(QApplication::translate("singlePic", "\346\200\273\347\211\207\346\225\260", Q_NULLPTR));
        label_2->setText(QApplication::translate("singlePic", "\345\215\225\351\235\242\347\211\207\346\225\260", Q_NULLPTR));
        label->setText(QApplication::translate("singlePic", "\345\275\223\345\211\215\347\211\207\346\225\260", Q_NULLPTR));
        lastPicButton->setText(QApplication::translate("singlePic", "\344\270\212\344\270\200\347\211\207", Q_NULLPTR));
        nextPicButton->setText(QApplication::translate("singlePic", "\344\270\213\344\270\200\347\211\207", Q_NULLPTR));
        closeButton->setText(QApplication::translate("singlePic", "\345\205\263\351\227\255", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class singlePic: public Ui_singlePic {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SINGLEPIC_H
