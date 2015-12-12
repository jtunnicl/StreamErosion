#
#  QCustomPlot Plot Examples
#

QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = plot-examples
TEMPLATE = app

SOURCES += main.cpp\
           qcustomplot.cpp \
    topo.cpp \
    streamerosion.cpp \
    flowrouting.cpp \
    illumin.cpp \
    mainwindow.cpp

HEADERS  += streamerosion.h \
            qcustomplot.h \
    sort.h \
    nr3.h \
    topo.h \
    flowrouting.h \
    illumin.h \
    mainwindow.h

FORMS    += \
    mainwindow.ui

