ROOT_DIR = ../..
MF_PREFIX=$(ROOT_DIR)/settings
include  $(MF_PREFIX)/Makefile.in

LOC_INCL_PATH = -I$(ROOT_DIR)/include
LINK_FLAGS += -L$(ROOT_DIR)/lib/unix
LINK_LIBS  += -l$(DIRECTFN_LIB_NAME)

TRIAG_TEST_ST = test_triag_st
TRIAG_TEST_EA = test_triag_ea
TRIAG_TEST_VA = test_triag_va
QUAD_TEST_ST  = test_quad_st
QUAD_TEST_EA  = test_quad_ea
QUAD_TEST_VA  = test_quad_va

.PHONY: all clean 

all: $(TRIAG_TEST_ST) \
     $(TRIAG_TEST_EA) \
     $(TRIAG_TEST_VA) \
     $(QUAD_TEST_ST) \
     $(QUAD_TEST_EA) \
     $(QUAD_TEST_VA) 

OBJ_TRG_ST   = test_directfn_triag_st.$(obj)
OBJ_TRG_EA   = test_directfn_triag_ea.$(obj)
OBJ_TRG_VA   = test_directfn_triag_va.$(obj)
OBJ_QUAD_ST  = test_directfn_quad_st.$(obj)
OBJ_QUAD_EA  = test_directfn_quad_ea.$(obj)
OBJ_QUAD_VA  = test_directfn_quad_va.$(obj)


$(TRIAG_TEST_ST): $(OBJ_TRG_ST)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)

$(TRIAG_TEST_EA): $(OBJ_TRG_EA)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)

$(TRIAG_TEST_VA): $(OBJ_TRG_VA)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)

$(QUAD_TEST_ST): $(OBJ_QUAD_ST)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)

$(QUAD_TEST_EA): $(OBJ_QUAD_EA) 
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)

$(QUAD_TEST_VA): $(OBJ_QUAD_VA)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)


%.$(obj): %.cpp 
	$(CompCXX) $(FlagsCXX) -c $(LOC_INCL_PATH) $< 

ALL_RM = $(TRIAG_TEST_ST) $(QUAD_TEST_ST) \
         $(TRIAG_TEST_EA) $(QUAD_TEST_EA) \
         $(TRIAG_TEST_VA) $(QUAD_TEST_VA) 

clean:
	rm *.$(obj) $(ALL_RM)

# End of file

