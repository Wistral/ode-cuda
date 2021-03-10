SUBDIRS = UnitTest++

AM_CPPFLAGS = -I $(srcdir)/UnitTest++/src \
              -I $(top_srcdir)/include \
              -I $(top_srcdir)/ode/src

if GIMPACT
    AM_CPPFLAGS += -DdTRIMESH_ENABLED -DdTRIMESH_GIMPACT
endif
if OPCODE
    AM_CPPFLAGS += -DdTRIMESH_ENABLED -DdTRIMESH_OPCODE
endif

LDADD = $(builddir)/UnitTest++/src/libunittestpp.la \
        $(top_builddir)/ode/src/libode.la \
        $(top_builddir)/ode/src/cuda_helper.lo \
        $(top_builddir)/ode/src/cuda_matrix.lo \
        -L/usr/local/cuda/lib64 -lcudart -lcuda -ldl -lrt

check_PROGRAMS = tests

TESTS = tests

tests_SOURCES = main.cpp joint.cpp odemath.cpp collision.cpp \
                joints/ball.cpp \
                joints/fixed.cpp \
                joints/hinge.cpp \
                joints/hinge2.cpp \
                joints/piston.cpp \
                joints/pr.cpp \
                joints/pu.cpp \
                joints/slider.cpp \
                joints/universal.cpp

