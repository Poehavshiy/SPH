#-------------------------------------------------
#
# Project created by QtCreator 2016-06-08T14:35:57
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = SPH_QUI
TEMPLATE = app
CONFIG += c++14


SOURCES = \
    MyQGraphicsView.cpp \
    MainWindow.cpp \
    main.cpp \
    SpaceParsing.cpp \
    Particile.cpp \
    Help.cpp \
    Flow.cpp \
    Calculator.cpp

HEADERS = \
    myqgraphicsview.h \
    MainWindow.h \
    SpaceParsing.h \
    Particile.h \
    Help.h \
    Flow.h \
    Calculator.h


FORMS    += mainwindow.ui

        #main.cpp
       # Particile.cpp Particile.h
       # Flow.cpp Flow.h
       # Help.cpp Help.h
       # SpaceParsing.cpp SpaceParsing.h
       # Calculator.cpp Calculator.h

DISTFILES += \
    testI_1.txt
