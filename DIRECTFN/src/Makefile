ROOT_DIR=..

MF_PREFIX=$(ROOT_DIR)/settings
include $(MF_PREFIX)/Makefile.in


.PHONY: all clean dependencies

LOCAL_LIB_NAME = lib$(DIRECTFN_LIB_NAME).$(SO)
LOCAL_INC_PATH += -I$(ROOT_DIR)/include


SOURCES = directfn_common.cpp \
          directfn_contour.cpp \
          directfn_gl_1d.cpp \
          directfn_greenfunc.cpp \
          directfn_greenfunc_voxhenry.cpp \
          directfn_kernel_base.cpp \
          directfn_kernel_tri.cpp \
          directfn_kernel_quad_geom.cpp \
          directfn_kernel_quad_scal.cpp \
          directfn_kernel_quad_vect.cpp \
          directfn_kernel_quad_voxhenry.cpp \
          directfn_kernel_arrsum.cpp \
          directfn_interface.cpp \
          directfn_algorithm_st.cpp \
          directfn_algorithm_ea.cpp \
          directfn_algorithm_va.cpp \
          directfn_algorithm_voxhenry.cpp \
          directfn_quad.cpp \
          directfn_quad_voxhenry.cpp \
          directfn_triag.cpp


HIP=$(ROOT_DIR)/include
HEADERS = $(HIP)/directfn_common.h \
          $(HIP)/directfn_contour.h \
          $(HIP)/directfn_defs.h \
          $(HIP)/directfn_greenfunc.h \
          $(HIP)/directfn_greenfunc_voxhenry.h \
          $(HIP)/directfn_kernel_base.h \
          $(HIP)/directfn_kernel_tri.h \
          $(HIP)/directfn_kernel_quad_geom.h \
          $(HIP)/directfn_kernel_quad_scal.h \
          $(HIP)/directfn_kernel_quad_vect.h \
          $(HIP)/directfn_kernel_quad_voxhenry.h \
          $(HIP)/directfn_kernel_arrsum.h \
          $(HIP)/directfn_interface.h \
          $(HIP)/directfn_algorithm_st.h \
          $(HIP)/directfn_algorithm_ea.h \
          $(HIP)/directfn_algorithm_va.h \
          $(HIP)/directfn_algorithm_voxhenry.h \
          $(HIP)/directfn_quad.h \
          $(HIP)/directfn_quad_voxhenry.h \
          $(HIP)/directfn_triag.h


OBJECTS = $(patsubst %.cpp,%.$(obj),$(SOURCES))

all: $(OBJECTS)
	$(AR) $(FlagsAR) $(LOCAL_LIB_NAME) $^
	$(MOVE) $(LOCAL_LIB_NAME) $(ROOT_DIR)/lib/unix


clean:
	$(RM) *.$(obj) dependencies

%.$(obj): %.cpp
	$(CompCXX) $(FlagsCXX) $(LOCAL_INC_PATH) -o $@ -c $< 

dependencies: $(SOURCES) $(HEADERS)
	$(CompCXX) $(FlagsCXX) -MM $(LOCAL_INC_PATH) $(SOURCES) >$@

include dependencies

# End of file
