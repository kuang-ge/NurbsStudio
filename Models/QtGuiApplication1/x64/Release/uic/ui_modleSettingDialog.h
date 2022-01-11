/********************************************************************************
** Form generated from reading UI file 'modleSettingDialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MODLESETTINGDIALOG_H
#define UI_MODLESETTINGDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QScrollBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_modleSettingDialog
{
public:
    QVBoxLayout *verticalLayout;
    QWidget *widget;
    QHBoxLayout *horizontalLayout;
    QLabel *lineWidth;
    QScrollBar *linWidthScrollBar;
    QFrame *frame;
    QHBoxLayout *horizontalLayout_3;
    QLabel *controlPts;
    QScrollBar *controlPtsScrollBar;
    QWidget *widget_2;
    QHBoxLayout *horizontalLayout_2;
    QLabel *Pts;
    QScrollBar *PtsScrollBar;

    void setupUi(QWidget *modleSettingDialog)
    {
        if (modleSettingDialog->objectName().isEmpty())
            modleSettingDialog->setObjectName(QStringLiteral("modleSettingDialog"));
        modleSettingDialog->resize(400, 300);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(modleSettingDialog->sizePolicy().hasHeightForWidth());
        modleSettingDialog->setSizePolicy(sizePolicy);
        modleSettingDialog->setMinimumSize(QSize(400, 300));
        modleSettingDialog->setMaximumSize(QSize(400, 300));
        verticalLayout = new QVBoxLayout(modleSettingDialog);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        widget = new QWidget(modleSettingDialog);
        widget->setObjectName(QStringLiteral("widget"));
        horizontalLayout = new QHBoxLayout(widget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        lineWidth = new QLabel(widget);
        lineWidth->setObjectName(QStringLiteral("lineWidth"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(100);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(lineWidth->sizePolicy().hasHeightForWidth());
        lineWidth->setSizePolicy(sizePolicy1);
        lineWidth->setMinimumSize(QSize(100, 0));

        horizontalLayout->addWidget(lineWidth);

        linWidthScrollBar = new QScrollBar(widget);
        linWidthScrollBar->setObjectName(QStringLiteral("linWidthScrollBar"));
        linWidthScrollBar->setMinimum(1);
        linWidthScrollBar->setMaximum(10);
        linWidthScrollBar->setOrientation(Qt::Horizontal);

        horizontalLayout->addWidget(linWidthScrollBar);


        verticalLayout->addWidget(widget);

        frame = new QFrame(modleSettingDialog);
        frame->setObjectName(QStringLiteral("frame"));
        frame->setFrameShape(QFrame::StyledPanel);
        frame->setFrameShadow(QFrame::Raised);
        horizontalLayout_3 = new QHBoxLayout(frame);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        controlPts = new QLabel(frame);
        controlPts->setObjectName(QStringLiteral("controlPts"));
        sizePolicy1.setHeightForWidth(controlPts->sizePolicy().hasHeightForWidth());
        controlPts->setSizePolicy(sizePolicy1);
        controlPts->setMinimumSize(QSize(100, 0));

        horizontalLayout_3->addWidget(controlPts);

        controlPtsScrollBar = new QScrollBar(frame);
        controlPtsScrollBar->setObjectName(QStringLiteral("controlPtsScrollBar"));
        controlPtsScrollBar->setMinimum(1);
        controlPtsScrollBar->setMaximum(20);
        controlPtsScrollBar->setOrientation(Qt::Horizontal);

        horizontalLayout_3->addWidget(controlPtsScrollBar);


        verticalLayout->addWidget(frame);

        widget_2 = new QWidget(modleSettingDialog);
        widget_2->setObjectName(QStringLiteral("widget_2"));
        horizontalLayout_2 = new QHBoxLayout(widget_2);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        Pts = new QLabel(widget_2);
        Pts->setObjectName(QStringLiteral("Pts"));
        sizePolicy1.setHeightForWidth(Pts->sizePolicy().hasHeightForWidth());
        Pts->setSizePolicy(sizePolicy1);
        Pts->setMinimumSize(QSize(100, 0));

        horizontalLayout_2->addWidget(Pts);

        PtsScrollBar = new QScrollBar(widget_2);
        PtsScrollBar->setObjectName(QStringLiteral("PtsScrollBar"));
        PtsScrollBar->setMinimum(1);
        PtsScrollBar->setMaximum(20);
        PtsScrollBar->setOrientation(Qt::Horizontal);

        horizontalLayout_2->addWidget(PtsScrollBar);


        verticalLayout->addWidget(widget_2);


        retranslateUi(modleSettingDialog);

        QMetaObject::connectSlotsByName(modleSettingDialog);
    } // setupUi

    void retranslateUi(QWidget *modleSettingDialog)
    {
        modleSettingDialog->setWindowTitle(QApplication::translate("modleSettingDialog", "modleSettingDialog", Q_NULLPTR));
        lineWidth->setText(QApplication::translate("modleSettingDialog", "\347\272\277\345\256\275", Q_NULLPTR));
        controlPts->setText(QApplication::translate("modleSettingDialog", "\346\216\247\345\210\266\347\202\271\345\244\247\345\260\217", Q_NULLPTR));
        Pts->setText(QApplication::translate("modleSettingDialog", "\346\230\276\347\244\272\347\202\271\345\244\247\345\260\217", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class modleSettingDialog: public Ui_modleSettingDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MODLESETTINGDIALOG_H
