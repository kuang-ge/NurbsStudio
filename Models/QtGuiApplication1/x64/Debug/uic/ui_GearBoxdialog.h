/********************************************************************************
** Form generated from reading UI file 'GearBoxdialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GEARBOXDIALOG_H
#define UI_GEARBOXDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QHeaderView>

QT_BEGIN_NAMESPACE

class Ui_GearBoxdialog
{
public:

    void setupUi(QDialog *GearBoxdialog)
    {
        if (GearBoxdialog->objectName().isEmpty())
            GearBoxdialog->setObjectName(QStringLiteral("GearBoxdialog"));
        GearBoxdialog->resize(400, 300);

        retranslateUi(GearBoxdialog);

        QMetaObject::connectSlotsByName(GearBoxdialog);
    } // setupUi

    void retranslateUi(QDialog *GearBoxdialog)
    {
        GearBoxdialog->setWindowTitle(QApplication::translate("GearBoxdialog", "GearBoxdialog", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class GearBoxdialog: public Ui_GearBoxdialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GEARBOXDIALOG_H
