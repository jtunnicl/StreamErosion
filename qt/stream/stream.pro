#-------------------------------------------------
#
# Project created by QtCreator 2015-11-17T15:11:37
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = stream
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    ../../flowrouting.cpp \
    ../../streamerosion.cpp \
    ../../topo.cpp

HEADERS  += mainwindow.h \
    ../../flowrouting.h \
    ../../streamerosion.h \
    ../../topo.h

FORMS    += mainwindow.ui
