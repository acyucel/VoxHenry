include ../../make.inc
INCL_PATH = -I../../include

allobj: a_functions_ws_st.($obj) lambda_limits_ws_st.($obj) n_functions_ws_st.($obj) psi_limits_ws_st.($obj) quadric_ws_st.($obj) subtriangles_ws_st.($obj) u_limits_ws_st.($obj)

a_functions_ws_st.($obj): a_functions_ws_st.cpp
	$(CompCXX) -c $(INCL_PATH) a_functions_ws_st.cpp

lambda_limits_ws_st.($obj): lambda_limits_ws_st.cpp
	$(CompCXX) -c $(INCL_PATH) lambda_limits_ws_st.cpp

n_functions_ws_st.($obj): n_functions_ws_st.cpp
	$(CompCXX) -c $(INCL_PATH) n_functions_ws_st.cpp

psi_limits_ws_st.($obj): psi_limits_ws_st.cpp
	$(CompCXX) -c $(INCL_PATH) psi_limits_ws_st.cpp

quadric_ws_st.($obj): quadric_ws_st.cpp
	$(CompCXX) -c $(INCL_PATH) quadric_ws_st.cpp

subtriangles_ws_st.($obj): subtriangles_ws_st.cpp
	$(CompCXX) -c $(INCL_PATH) subtriangles_ws_st.cpp

u_limits_ws_st.($obj): u_limits_ws_st.cpp
	$(CompCXX) -c $(INCL_PATH) u_limits_ws_st.cpp

clean:
	rm *.$(obj)