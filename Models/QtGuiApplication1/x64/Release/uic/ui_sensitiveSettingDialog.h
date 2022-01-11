/********************************************************************************
** Form generated from reading UI file 'sensitiveSettingDialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SENSITIVESETTINGDIALOG_H
#define UI_SENSITIVESETTINGDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QScrollBar>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_sensitiveSettingDialog
{
public:
    QVBoxLayout *verticalLayout;
    QWidget *widget_2;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_2;
    QSpinBox *keyBordSpinBox;
    QScrollBar *keyBordScrollBar;
    QWidget *widget;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QSpinBox *mouseSpinBox;
    QScrollBar *mouseScrollBar;
    QWidget *widget_3;
    QHBoxLayout *horizontalLayout_3;
    QSpacerItem *horizontalSpacer;
    QPushButton *confirmButton;
    QPushButton *cancleButton;
    QSpacerItem *horizontalSpacer_2;

    void setupUi(QWidget *sensitiveSettingDialog)
    {
        if (sensitiveSettingDialog->objectName().isEmpty())
            sensitiveSettingDialog->setObjectName(QStringLiteral("sensitiveSettingDialog"));
        sensitiveSettingDialog->resize(400, 300);
        QSizePolicy sizePolicy(QSizePolicy::Maximum, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(sensitiveSettingDialog->sizePolicy().hasHeightForWidth());
        sensitiveSettingDialog->setSizePolicy(sizePolicy);
        sensitiveSettingDialog->setMinimumSize(QSize(400, 300));
        sensitiveSettingDialog->setMaximumSize(QSize(400, 300));
        verticalLayout = new QVBoxLayout(sensitiveSettingDialog);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        widget_2 = new QWidget(sensitiveSettingDialog);
        widget_2->setObjectName(QStringLiteral("widget_2"));
        horizontalLayout_2 = new QHBoxLayout(widget_2);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label_2 = new QLabel(widget_2);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setMinimumSize(QSize(80, 20));
        label_2->setMaximumSize(QSize(80, 20));

        horizontalLayout_2->addWidget(label_2);

        keyBordSpinBox = new QSpinBox(widget_2);
        keyBordSpinBox->setObjectName(QStringLiteral("keyBordSpinBox"));
        keyBordSpinBox->setMinimumSize(QSize(50, 0));
        keyBordSpinBox->setMaximumSize(QSize(50, 16777215));
        keyBordSpinBox->setMinimum(1);
        keyBordSpinBox->setMaximum(1000);

        horizontalLayout_2->addWidget(keyBordSpinBox);

        keyBordScrollBar = new QScrollBar(widget_2);
        keyBordScrollBar->setObjectName(QStringLiteral("keyBordScrollBar"));
        keyBordScrollBar->setMinimum(1);
        keyBordScrollBar->setMaximum(1000);
        keyBordScrollBar->setOrientation(Qt::Horizontal);

        horizontalLayout_2->addWidget(keyBordScrollBar);


        verticalLayout->addWidget(widget_2);

        widget = new QWidget(sensitiveSettingDialog);
        widget->setObjectName(QStringLiteral("widget"));
        horizontalLayout = new QHBoxLayout(widget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        label = new QLabel(widget);
        label->setObjectName(QStringLiteral("label"));
        label->setMinimumSize(QSize(80, 0));
        label->setMaximumSize(QSize(80, 16777215));

        horizontalLayout->addWidget(label);

        mouseSpinBox = new QSpinBox(widget);
        mouseSpinBox->setObjectName(QStringLiteral("mouseSpinBox"));
        mouseSpinBox->setMinimumSize(QSize(50, 0));
        mouseSpinBox->setMaximumSize(QSize(50, 16777215));
        mouseSpinBox->setMinimum(1);
        mouseSpinBox->setMaximum(10);

        horizontalLayout->addWidget(mouseSpinBox);

        mouseScrollBar = new QScrollBar(widget);
        mouseScrollBar->setObjectName(QStringLiteral("mouseScrollBar"));
        mouseScrollBar->setMinimum(1);
        mouseScrollBar->setMaximum(10);
        mouseScrollBar->setOrientation(Qt::Horizontal);

        horizontalLayout->addWidget(mouseScrollBar);


        verticalLayout->addWidget(widget);

        widget_3 = new QWidget(sensitiveSettingDialog);
        widget_3->setObjectName(QStringLiteral("widget_3"));
        horizontalLayout_3 = new QHBoxLayout(widget_3);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer);

        confirmButton = new QPushButton(widget_3);
        confirmButton->setObjectName(QStringLiteral("confirmButton"));

        horizontalLayout_3->addWidget(confirmButton);

        cancleButton = new QPushButton(widget_3);
        cancleButton->setObjectName(QStringLiteral("cancleButton"));

        horizontalLayout_3->addWidget(cancleButton);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_2);


        verticalLayout->addWidget(widget_3);


        retranslateUi(sensitiveSettingDialog);

        QMetaObject::connectSlotsByName(sensitiveSettingDialog);
    } // setupUi

    void retranslateUi(QWidget *sensitiveSettingDialog)
    {
        sensitiveSettingDialog->setWindowTitle(QApplication::translate("sensitiveSettingDialog", "\347\201\265\346\225\217\345\272\246\350\256\276\347\275\256", Q_NULLPTR));
        label_2->setText(QApplication::translate("sensitiveSettingDialog", "\351\224\256\347\233\230\347\201\265\346\225\217\345\272\246", Q_NULLPTR));
        label->setText(QApplication::translate("sensitiveSettingDialog", "\351\274\240\346\240\207\347\201\265\346\225\217\345\272\246", Q_NULLPTR));
        confirmButton->setText(QApplication::translate("sensitiveSettingDialog", "\347\241\256\345\256\232", Q_NULLPTR));
        cancleButton->setText(QApplication::translate("sensitiveSettingDialog", "\345\217\226\346\266\210", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class sensitiveSettingDialog: public Ui_sensitiveSettingDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SENSITIVESETTINGDIALOG_H
