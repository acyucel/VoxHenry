ROOT_DIR = ../../..
MF_PREFIX=$(ROOT_DIR)/settings
include  $(MF_PREFIX)/Makefile.in

LOC_INCL_PATH = -I$(ROOT_DIR)/include
LINK_FLAGS += -L$(ROOT_DIR)/lib/unix
LINK_LIBS  += -l$(DIRECTFN_LIB_NAME)

TRI_VS_QUAD = tri_vs_quad
ST_ONE_PARAM_CONST = st_one_param_const
ST_CURV = st_curv


.PHONY: all clean 

all: $(TRI_VS_QUAD) \
     $(ST_ONE_PARAM_CONST) \
     $(ST_CURV) 
     #$(QUAD_TEST_ST) \
     #$(QUAD_TEST_EA) \
     #$(QUAD_TEST_VA) 

OBJ_TRI_VS_QUAD   = tri_vs_quad.$(obj)
OBJ_ST_ONE_PARAM_CONST   = st_one_param_const.$(obj)
OBJ_ST_CURV = st_curv.$(obj)



$(TRI_VS_QUAD): $(OBJ_TRI_VS_QUAD)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)

$(ST_ONE_PARAM_CONST): $(OBJ_ST_ONE_PARAM_CONST)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)
	
$(ST_CURV): $(OBJ_ST_CURV)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)
	



%.$(obj): %.cpp 
	$(CompCXX) $(FlagsCXX) -c $(LOC_INCL_PATH) $< 

ALL_RM = $(TRI_VS_QUAD) $(ST_ONE_PARAM_CONST) $(ST_CURV)

clean:
	rm *.$(obj)

# End of file

