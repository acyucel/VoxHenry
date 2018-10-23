ROOT_DIR = ../../..
MF_PREFIX=$(ROOT_DIR)/settings
include  $(MF_PREFIX)/Makefile.in

LOC_INCL_PATH = -I$(ROOT_DIR)/include
LINK_FLAGS += -L$(ROOT_DIR)/lib/unix
LINK_LIBS  += -l$(DIRECTFN_LIB_NAME)

TRI_VS_QUAD = tri_vs_quad
WS_ST_ONE_N = ws_st_one_N
SS_EA_ONE_N = ss_ea_one_N
SS_VA_ONE_N = ss_va_one_N
CURV = curvilinear

.PHONY: all clean 

all: $(TRI_VS_QUAD) \
     $(WS_ST_ONE_N) \
     $(SS_EA_ONE_N) \
     $(SS_VA_ONE_N) \
     $(CURV)  

OBJ_TRI_VS_QUAD   = tri_vs_quad.$(obj)
OBJ_WS_ST_ONE_N   = ws_st_one_N.$(obj)
OBJ_SS_EA_ONE_N   = ss_ea_one_N.$(obj)
OBJ_SS_VA_ONE_N  = ss_va_one_N.$(obj)
OBJ_CURV = curvilinear.$(obj)


$(TRI_VS_QUAD): $(OBJ_TRI_VS_QUAD)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)

$(WS_ST_ONE_N): $(OBJ_WS_ST_ONE_N)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)

$(SS_EA_ONE_N): $(OBJ_SS_EA_ONE_N)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)
	
$(SS_VA_ONE_N): $(OBJ_SS_VA_ONE_N)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)

$(CURV): $(OBJ_CURV)
	$(CompCXX) $(FlagsCXX) -o $@  $^ $(LINK_FLAGS) $(LINK_LIBS)

%.$(obj): %.cpp 
	$(CompCXX) $(FlagsCXX) -c $(LOC_INCL_PATH) $< 

ALL_RM = $(TRI_VS_QUAD) $(WS_ST_ONE_N) \
         $(SS_EA_ONE_N) $(SS_VA_ONE_N) \
		 $(CURV)

clean:
	rm *.$(obj) 

# End of file

