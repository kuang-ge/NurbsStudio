/********************************************************************************
** Form generated from reading UI file 'readFileDialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_READFILEDIALOG_H
#define UI_READFILEDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_readFileDialog
{
public:
    QVBoxLayout *verticalLayout;
    QWidget *widget_2;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QLineEdit *filePath;
    QPushButton *selectFile;
    QWidget *widget_3;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_2;
    QComboBox *modelType;
    QWidget *widget;
    QGridLayout *gridLayout;
    QPushButton *conform;
    QPushButton *cancle;
    QSpacerItem *horizontalSpacer;
    QSpacerItem *horizontalSpacer_2;

    void setupUi(QWidget *readFileDialog)
    {
        if (readFileDialog->objectName().isEmpty())
            readFileDialog->setObjectName(QStringLiteral("readFileDialog"));
        readFileDialog->resize(400, 200);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(readFileDialog->sizePolicy().hasHeightForWidth());
        readFileDialog->setSizePolicy(sizePolicy);
        readFileDialog->setMinimumSize(QSize(400, 200));
        readFileDialog->setMaximumSize(QSize(400, 200));
        verticalLayout = new QVBoxLayout(readFileDialog);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        widget_2 = new QWidget(readFileDialog);
        widget_2->setObjectName(QStringLiteral("widget_2"));
        sizePolicy.setHeightForWidth(widget_2->sizePolicy().hasHeightForWidth());
        widget_2->setSizePolicy(sizePolicy);
        horizontalLayout = new QHBoxLayout(widget_2);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        label = new QLabel(widget_2);
        label->setObjectName(QStringLiteral("label"));

        horizontalLayout->addWidget(label);

        filePath = new QLineEdit(widget_2);
        filePath->setObjectName(QStringLiteral("filePath"));

        horizontalLayout->addWidget(filePath);

        selectFile = new QPushButton(widget_2);
        selectFile->setObjectName(QStringLiteral("selectFile"));

        horizontalLayout->addWidget(selectFile);


        verticalLayout->addWidget(widget_2);

        widget_3 = new QWidget(readFileDialog);
        widget_3->setObjectName(QStringLiteral("widget_3"));
        horizontalLayout_2 = new QHBoxLayout(widget_3);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label_2 = new QLabel(widget_3);
        label_2->setObjectName(QStringLiteral("label_2"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(label_2->sizePolicy().hasHeightForWidth());
        label_2->setSizePolicy(sizePolicy1);

        horizontalLayout_2->addWidget(label_2);

        modelType = new QComboBox(widget_3);
        modelType->setObjectName(QStringLiteral("modelType"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(modelType->sizePolicy().hasHeightForWidth());
        modelType->setSizePolicy(sizePolicy2);

        horizontalLayout_2->addWidget(modelType);


        verticalLayout->addWidget(widget_3);

        widget = new QWidget(readFileDialog);
        widget->setObjectName(QStringLiteral("widget"));
        gridLayout = new QGridLayout(widget);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        conform = new QPushButton(widget);
        conform->setObjectName(QStringLiteral("conform"));

        gridLayout->addWidget(conform, 0, 1, 1, 1);

        cancle = new QPushButton(widget);
        cancle->setObjectName(QStringLiteral("cancle"));

        gridLayout->addWidget(cancle, 0, 2, 1, 1);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer, 0, 0, 1, 1);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_2, 0, 3, 1, 1);


        verticalLayout->addWidget(widget);


        retranslateUi(readFileDialog);

        QMetaObject::connectSlotsByName(readFileDialog);
    } // setupUi

    void retranslateUi(QWidget *readFileDialog)
    {
        readFileDialog->setWindowTitle(QApplication::translate("readFileDialog", "readFileDialog", Q_NULLPTR));
        label->setText(QApplication::translate("readFileDialog", "\346\226\207\344\273\266\350\267\257\345\276\204", Q_NULLPTR));
        selectFile->setText(QApplication::translate("readFileDialog", "\351\200\211\346\213\251\346\226\207\344\273\266", Q_NULLPTR));
        label_2->setText(QApplication::translate("readFileDialog", "\346\250\241\345\236\213\345\261\236\346\200\247", Q_NULLPTR));
        modelType->clear();
        modelType->insertItems(0, QStringList()
         << QApplication::translate("readFileDialog", "point3d", Q_NULLPTR)
         << QApplication::translate("readFileDialog", "point4d", Q_NULLPTR)
         << QApplication::translate("readFileDialog", "bezierLine", Q_NULLPTR)
         << QApplication::translate("readFileDialog", "nurbusLine", Q_NULLPTR)
         << QApplication::translate("readFileDialog", "nurbusSurface", Q_NULLPTR)
         << QApplication::translate("readFileDialog", "nurbusVol", Q_NULLPTR)
        );
        conform->setText(QApplication::translate("readFileDialog", "\347\241\256\345\256\232", Q_NULLPTR));
        cancle->setText(QApplication::translate("readFileDialog", "\345\217\226\346\266\210", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class readFileDialog: public Ui_readFileDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_READFILEDIALOG_H
