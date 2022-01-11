/********************************************************************************
** Form generated from reading UI file 'QtGuiApplication1.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QTGUIAPPLICATION1_H
#define UI_QTGUIAPPLICATION1_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QTreeWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include "win.h"

QT_BEGIN_NAMESPACE

class Ui_QtGuiApplication1Class
{
public:
    QAction *actionreadFile;
    QAction *actionnewFile;
    QAction *actionsenstive;
    QAction *actionsaveAs3MF;
    QAction *actionPoints;
    QAction *actionLines;
    QAction *actionsurface;
    QAction *actionmodelConfigure;
    QAction *actionmodelSetting;
    QAction *actionFragment;
    QAction *actioncoordinate;
    QAction *actionbordLineColor;
    QAction *actionsurfaceColor;
    QAction *actionControlpointColor;
    QAction *actionsliceSurfaceColor;
    QAction *actionsliceSupportColor;
    QAction *actioncontrolPoints;
    QAction *actionsinglePic;
    QAction *actionchangeFineScale;
    QAction *actioncameraLocation;
    QAction *actiongear;
    QAction *actionReducer;
    QAction *actionGearbox;
    QAction *actionBearingBlock;
    QAction *actionclear_all;
    QAction *actiontools;
    win *centralWidget;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuOption;
    QMenu *menumodelType;
    QMenu *menuSet;
    QMenu *menu_3;
    QMenu *menu;
    QMenu *menu_4;
    QMenu *menu_2;
    QToolBar *toolBar;
    QDockWidget *dockWidget_tools;
    QWidget *dockWidgetContents_tools;
    QFormLayout *formLayout;
    QPushButton *initCameraLocationButton;
    QDockWidget *dockWidget_data;
    QWidget *dockWidgetContents_data;
    QVBoxLayout *verticalLayout;
    QTreeWidget *treeWidget;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *QtGuiApplication1Class)
    {
        if (QtGuiApplication1Class->objectName().isEmpty())
            QtGuiApplication1Class->setObjectName(QStringLiteral("QtGuiApplication1Class"));
        QtGuiApplication1Class->setWindowModality(Qt::NonModal);
        QtGuiApplication1Class->resize(1075, 1013);
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(QtGuiApplication1Class->sizePolicy().hasHeightForWidth());
        QtGuiApplication1Class->setSizePolicy(sizePolicy);
        QtGuiApplication1Class->setMinimumSize(QSize(0, 28));
        QtGuiApplication1Class->setMaximumSize(QSize(2000, 2000));
        actionreadFile = new QAction(QtGuiApplication1Class);
        actionreadFile->setObjectName(QStringLiteral("actionreadFile"));
        actionnewFile = new QAction(QtGuiApplication1Class);
        actionnewFile->setObjectName(QStringLiteral("actionnewFile"));
        actionsenstive = new QAction(QtGuiApplication1Class);
        actionsenstive->setObjectName(QStringLiteral("actionsenstive"));
        actionsenstive->setEnabled(true);
        actionsaveAs3MF = new QAction(QtGuiApplication1Class);
        actionsaveAs3MF->setObjectName(QStringLiteral("actionsaveAs3MF"));
        actionsaveAs3MF->setCheckable(false);
        actionsaveAs3MF->setEnabled(false);
        actionPoints = new QAction(QtGuiApplication1Class);
        actionPoints->setObjectName(QStringLiteral("actionPoints"));
        actionPoints->setCheckable(true);
        actionLines = new QAction(QtGuiApplication1Class);
        actionLines->setObjectName(QStringLiteral("actionLines"));
        actionLines->setCheckable(true);
        actionsurface = new QAction(QtGuiApplication1Class);
        actionsurface->setObjectName(QStringLiteral("actionsurface"));
        actionsurface->setCheckable(true);
        actionsurface->setChecked(true);
        actionmodelConfigure = new QAction(QtGuiApplication1Class);
        actionmodelConfigure->setObjectName(QStringLiteral("actionmodelConfigure"));
        actionmodelSetting = new QAction(QtGuiApplication1Class);
        actionmodelSetting->setObjectName(QStringLiteral("actionmodelSetting"));
        actionFragment = new QAction(QtGuiApplication1Class);
        actionFragment->setObjectName(QStringLiteral("actionFragment"));
        actionFragment->setEnabled(false);
        actioncoordinate = new QAction(QtGuiApplication1Class);
        actioncoordinate->setObjectName(QStringLiteral("actioncoordinate"));
        actioncoordinate->setCheckable(true);
        actioncoordinate->setChecked(true);
        actioncoordinate->setEnabled(true);
        actionbordLineColor = new QAction(QtGuiApplication1Class);
        actionbordLineColor->setObjectName(QStringLiteral("actionbordLineColor"));
        actionsurfaceColor = new QAction(QtGuiApplication1Class);
        actionsurfaceColor->setObjectName(QStringLiteral("actionsurfaceColor"));
        actionControlpointColor = new QAction(QtGuiApplication1Class);
        actionControlpointColor->setObjectName(QStringLiteral("actionControlpointColor"));
        actionsliceSurfaceColor = new QAction(QtGuiApplication1Class);
        actionsliceSurfaceColor->setObjectName(QStringLiteral("actionsliceSurfaceColor"));
        actionsliceSurfaceColor->setEnabled(false);
        actionsliceSupportColor = new QAction(QtGuiApplication1Class);
        actionsliceSupportColor->setObjectName(QStringLiteral("actionsliceSupportColor"));
        actionsliceSupportColor->setEnabled(false);
        actioncontrolPoints = new QAction(QtGuiApplication1Class);
        actioncontrolPoints->setObjectName(QStringLiteral("actioncontrolPoints"));
        actioncontrolPoints->setCheckable(true);
        actionsinglePic = new QAction(QtGuiApplication1Class);
        actionsinglePic->setObjectName(QStringLiteral("actionsinglePic"));
        actionsinglePic->setEnabled(false);
        actionchangeFineScale = new QAction(QtGuiApplication1Class);
        actionchangeFineScale->setObjectName(QStringLiteral("actionchangeFineScale"));
        actioncameraLocation = new QAction(QtGuiApplication1Class);
        actioncameraLocation->setObjectName(QStringLiteral("actioncameraLocation"));
        actiongear = new QAction(QtGuiApplication1Class);
        actiongear->setObjectName(QStringLiteral("actiongear"));
        actionReducer = new QAction(QtGuiApplication1Class);
        actionReducer->setObjectName(QStringLiteral("actionReducer"));
        actionGearbox = new QAction(QtGuiApplication1Class);
        actionGearbox->setObjectName(QStringLiteral("actionGearbox"));
        actionBearingBlock = new QAction(QtGuiApplication1Class);
        actionBearingBlock->setObjectName(QStringLiteral("actionBearingBlock"));
        actionclear_all = new QAction(QtGuiApplication1Class);
        actionclear_all->setObjectName(QStringLiteral("actionclear_all"));
        actiontools = new QAction(QtGuiApplication1Class);
        actiontools->setObjectName(QStringLiteral("actiontools"));
        actiontools->setCheckable(true);
        centralWidget = new win(QtGuiApplication1Class);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        centralWidget->setMaximumSize(QSize(2000, 2000));
        QtGuiApplication1Class->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(QtGuiApplication1Class);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1075, 30));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QStringLiteral("menuFile"));
        menuOption = new QMenu(menuBar);
        menuOption->setObjectName(QStringLiteral("menuOption"));
        menumodelType = new QMenu(menuOption);
        menumodelType->setObjectName(QStringLiteral("menumodelType"));
        menuSet = new QMenu(menuBar);
        menuSet->setObjectName(QStringLiteral("menuSet"));
        menu_3 = new QMenu(menuSet);
        menu_3->setObjectName(QStringLiteral("menu_3"));
        menu = new QMenu(menuBar);
        menu->setObjectName(QStringLiteral("menu"));
        menu_4 = new QMenu(menu);
        menu_4->setObjectName(QStringLiteral("menu_4"));
        menu_2 = new QMenu(menuBar);
        menu_2->setObjectName(QStringLiteral("menu_2"));
        QtGuiApplication1Class->setMenuBar(menuBar);
        toolBar = new QToolBar(QtGuiApplication1Class);
        toolBar->setObjectName(QStringLiteral("toolBar"));
        QtGuiApplication1Class->addToolBar(Qt::TopToolBarArea, toolBar);
        dockWidget_tools = new QDockWidget(QtGuiApplication1Class);
        dockWidget_tools->setObjectName(QStringLiteral("dockWidget_tools"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(dockWidget_tools->sizePolicy().hasHeightForWidth());
        dockWidget_tools->setSizePolicy(sizePolicy1);
        dockWidget_tools->setMinimumSize(QSize(147, 90));
        dockWidgetContents_tools = new QWidget();
        dockWidgetContents_tools->setObjectName(QStringLiteral("dockWidgetContents_tools"));
        formLayout = new QFormLayout(dockWidgetContents_tools);
        formLayout->setSpacing(6);
        formLayout->setContentsMargins(11, 11, 11, 11);
        formLayout->setObjectName(QStringLiteral("formLayout"));
        initCameraLocationButton = new QPushButton(dockWidgetContents_tools);
        initCameraLocationButton->setObjectName(QStringLiteral("initCameraLocationButton"));

        formLayout->setWidget(0, QFormLayout::LabelRole, initCameraLocationButton);

        dockWidget_tools->setWidget(dockWidgetContents_tools);
        QtGuiApplication1Class->addDockWidget(static_cast<Qt::DockWidgetArea>(1), dockWidget_tools);
        dockWidget_data = new QDockWidget(QtGuiApplication1Class);
        dockWidget_data->setObjectName(QStringLiteral("dockWidget_data"));
        dockWidgetContents_data = new QWidget();
        dockWidgetContents_data->setObjectName(QStringLiteral("dockWidgetContents_data"));
        verticalLayout = new QVBoxLayout(dockWidgetContents_data);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        treeWidget = new QTreeWidget(dockWidgetContents_data);
        QTreeWidgetItem *__qtreewidgetitem = new QTreeWidgetItem();
        __qtreewidgetitem->setText(0, QStringLiteral("1"));
        treeWidget->setHeaderItem(__qtreewidgetitem);
        treeWidget->setObjectName(QStringLiteral("treeWidget"));

        verticalLayout->addWidget(treeWidget);

        dockWidget_data->setWidget(dockWidgetContents_data);
        QtGuiApplication1Class->addDockWidget(static_cast<Qt::DockWidgetArea>(1), dockWidget_data);
        statusBar = new QStatusBar(QtGuiApplication1Class);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        QtGuiApplication1Class->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuOption->menuAction());
        menuBar->addAction(menuSet->menuAction());
        menuBar->addAction(menu->menuAction());
        menuBar->addAction(menu_2->menuAction());
        menuFile->addAction(actionreadFile);
        menuFile->addAction(actionnewFile);
        menuFile->addSeparator();
        menuFile->addAction(actionsaveAs3MF);
        menuOption->addAction(menumodelType->menuAction());
        menuOption->addAction(actionmodelSetting);
        menuOption->addAction(actionsinglePic);
        menuOption->addSeparator();
        menuOption->addAction(actioncontrolPoints);
        menuOption->addSeparator();
        menuOption->addAction(actionFragment);
        menumodelType->addAction(actionPoints);
        menumodelType->addAction(actionLines);
        menumodelType->addAction(actionsurface);
        menuSet->addAction(actionsenstive);
        menuSet->addSeparator();
        menuSet->addAction(menu_3->menuAction());
        menuSet->addAction(actioncameraLocation);
        menuSet->addSeparator();
        menuSet->addAction(actioncoordinate);
        menuSet->addAction(actionchangeFineScale);
        menu_3->addAction(actionbordLineColor);
        menu_3->addAction(actionsurfaceColor);
        menu_3->addAction(actionControlpointColor);
        menu_3->addAction(actionsliceSurfaceColor);
        menu_3->addAction(actionsliceSupportColor);
        menu->addAction(menu_4->menuAction());
        menu_4->addAction(actiongear);
        menu_4->addAction(actionReducer);
        menu_4->addAction(actionGearbox);
        menu_4->addAction(actionBearingBlock);
        toolBar->addSeparator();
        toolBar->addAction(actionmodelSetting);
        toolBar->addAction(actionPoints);
        toolBar->addAction(actionLines);
        toolBar->addAction(actionsurface);
        toolBar->addAction(actioncontrolPoints);
        toolBar->addAction(actionsinglePic);
        toolBar->addAction(actionchangeFineScale);
        toolBar->addAction(actioncoordinate);
        toolBar->addAction(actionclear_all);

        retranslateUi(QtGuiApplication1Class);

        QMetaObject::connectSlotsByName(QtGuiApplication1Class);
    } // setupUi

    void retranslateUi(QMainWindow *QtGuiApplication1Class)
    {
        QtGuiApplication1Class->setWindowTitle(QApplication::translate("QtGuiApplication1Class", "\344\275\223\345\217\202\346\225\260\345\214\226\345\273\272\346\250\241", Q_NULLPTR));
        actionreadFile->setText(QApplication::translate("QtGuiApplication1Class", "\350\257\273\345\217\226\346\226\207\344\273\266", Q_NULLPTR));
        actionnewFile->setText(QApplication::translate("QtGuiApplication1Class", "\346\226\260\346\226\207\344\273\266", Q_NULLPTR));
        actionsenstive->setText(QApplication::translate("QtGuiApplication1Class", "\347\201\265\346\225\217\345\272\246", Q_NULLPTR));
        actionsaveAs3MF->setText(QApplication::translate("QtGuiApplication1Class", "\345\217\246\345\255\230\344\270\2723MF", Q_NULLPTR));
        actionPoints->setText(QApplication::translate("QtGuiApplication1Class", "\347\246\273\346\225\243\347\202\271\346\230\276\347\244\272", Q_NULLPTR));
        actionLines->setText(QApplication::translate("QtGuiApplication1Class", "\347\275\221\346\240\274\346\230\276\347\244\272", Q_NULLPTR));
        actionsurface->setText(QApplication::translate("QtGuiApplication1Class", "\351\235\242\347\211\207\346\230\276\347\244\272", Q_NULLPTR));
        actionmodelConfigure->setText(QApplication::translate("QtGuiApplication1Class", "\346\250\241\345\236\213\345\217\202\346\225\260", Q_NULLPTR));
        actionmodelSetting->setText(QApplication::translate("QtGuiApplication1Class", "\346\250\241\345\236\213\345\217\202\346\225\260", Q_NULLPTR));
        actionFragment->setText(QApplication::translate("QtGuiApplication1Class", "\345\210\207\347\211\207", Q_NULLPTR));
        actioncoordinate->setText(QApplication::translate("QtGuiApplication1Class", "\345\235\220\346\240\207\347\263\273", Q_NULLPTR));
        actionbordLineColor->setText(QApplication::translate("QtGuiApplication1Class", "\350\276\271\347\225\214\347\272\277\351\242\234\350\211\262", Q_NULLPTR));
        actionsurfaceColor->setText(QApplication::translate("QtGuiApplication1Class", "\351\235\242\347\211\207\351\242\234\350\211\262", Q_NULLPTR));
        actionControlpointColor->setText(QApplication::translate("QtGuiApplication1Class", "\346\216\247\345\210\266\347\202\271\351\242\234\350\211\262", Q_NULLPTR));
        actionsliceSurfaceColor->setText(QApplication::translate("QtGuiApplication1Class", "\345\210\207\347\211\207\351\235\242\351\242\234\350\211\262", Q_NULLPTR));
        actionsliceSupportColor->setText(QApplication::translate("QtGuiApplication1Class", "\345\210\207\347\211\207\346\224\257\346\222\221\351\242\234\350\211\262", Q_NULLPTR));
        actioncontrolPoints->setText(QApplication::translate("QtGuiApplication1Class", "\346\216\247\345\210\266\347\202\271", Q_NULLPTR));
        actionsinglePic->setText(QApplication::translate("QtGuiApplication1Class", "\345\215\225\347\211\207\346\230\276\347\244\272", Q_NULLPTR));
        actionchangeFineScale->setText(QApplication::translate("QtGuiApplication1Class", "\346\233\264\346\224\271\347\273\206\345\210\206\345\272\246", Q_NULLPTR));
        actioncameraLocation->setText(QApplication::translate("QtGuiApplication1Class", "\346\221\204\345\203\217\346\234\272\344\275\215\347\275\256", Q_NULLPTR));
        actiongear->setText(QApplication::translate("QtGuiApplication1Class", "Spur gear", Q_NULLPTR));
        actionReducer->setText(QApplication::translate("QtGuiApplication1Class", "Reducer", Q_NULLPTR));
        actionGearbox->setText(QApplication::translate("QtGuiApplication1Class", "Gearbox", Q_NULLPTR));
        actionBearingBlock->setText(QApplication::translate("QtGuiApplication1Class", "Bearing Block", Q_NULLPTR));
        actionclear_all->setText(QApplication::translate("QtGuiApplication1Class", "clear all", Q_NULLPTR));
        actiontools->setText(QApplication::translate("QtGuiApplication1Class", "\345\274\200\345\217\221\345\267\245\345\205\267", Q_NULLPTR));
        menuFile->setTitle(QApplication::translate("QtGuiApplication1Class", "\346\226\207\344\273\266", Q_NULLPTR));
        menuOption->setTitle(QApplication::translate("QtGuiApplication1Class", "\351\200\211\351\241\271", Q_NULLPTR));
        menumodelType->setTitle(QApplication::translate("QtGuiApplication1Class", "\346\250\241\345\236\213\346\230\276\347\244\272\345\275\242\345\274\217", Q_NULLPTR));
        menuSet->setTitle(QApplication::translate("QtGuiApplication1Class", "\350\256\276\347\275\256", Q_NULLPTR));
        menu_3->setTitle(QApplication::translate("QtGuiApplication1Class", "\351\242\234\350\211\262", Q_NULLPTR));
        menu->setTitle(QApplication::translate("QtGuiApplication1Class", "\345\267\245\345\205\267\346\240\217", Q_NULLPTR));
        menu_4->setTitle(QApplication::translate("QtGuiApplication1Class", "FeatureNetwork", Q_NULLPTR));
        menu_2->setTitle(QApplication::translate("QtGuiApplication1Class", "\345\270\256\345\212\251", Q_NULLPTR));
        toolBar->setWindowTitle(QApplication::translate("QtGuiApplication1Class", "toolBar", Q_NULLPTR));
        dockWidget_tools->setWindowTitle(QApplication::translate("QtGuiApplication1Class", "\345\267\245\345\205\267", Q_NULLPTR));
        initCameraLocationButton->setText(QApplication::translate("QtGuiApplication1Class", "\350\247\206\350\247\222\350\277\230\345\216\237", Q_NULLPTR));
        dockWidget_data->setWindowTitle(QApplication::translate("QtGuiApplication1Class", "\346\225\260\346\215\256", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class QtGuiApplication1Class: public Ui_QtGuiApplication1Class {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QTGUIAPPLICATION1_H
