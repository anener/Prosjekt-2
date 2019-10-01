TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    bisection.cpp \
        jacobi.cpp \
    lanczosalgoritme.cpp \
        main.cpp \

INCLUDEPATH += C:\armadillo-9.700.2\include
DEPENDPATH += C:\armadillo-9.700.2\include

LIBS += \
    -LC:\armadillo-9.700.2\examples\lib_win64 \
    -llapack_win64_MT \
    -lblas_win64_MT

HEADERS += \ \
    bisectionmetode.h \
    catch.hpp \
    jacobi.h \
    lanczo.h

