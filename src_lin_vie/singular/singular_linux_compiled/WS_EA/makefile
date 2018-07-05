include ../../make.inc
INCL_PATH = -I../../include

allobj: A_functions_ws_ea.($obj) N_functions_ws_ea.($obj) PSI_limits_ws_ea.($obj) quadric_ws_ea.($obj) Simplex_ws_ea.($obj) THETA_limits_ws_ea.($obj)

A_functions_ws_ea.($obj): A_functions_ws_ea.cpp
	$(CompCXX) -c $(INCL_PATH) A_functions_ws_ea.cpp

N_functions_ws_ea.($obj): N_functions_ws_ea.cpp
	$(CompCXX) -c $(INCL_PATH) N_functions_ws_ea.cpp

PSI_limits_ws_ea.($obj): PSI_limits_ws_ea.cpp
	$(CompCXX) -c $(INCL_PATH) PSI_limits_ws_ea.cpp

quadric_ws_ea.($obj): quadric_ws_ea.cpp
	$(CompCXX) -c $(INCL_PATH) quadric_ws_ea.cpp

Simplex_ws_ea.($obj): Simplex_ws_ea.cpp
	$(CompCXX) -c $(INCL_PATH) Simplex_ws_ea.cpp

THETA_limits_ws_ea.($obj): THETA_limits_ws_ea.cpp
	$(CompCXX) -c $(INCL_PATH) THETA_limits_ws_ea.cpp

clean:
	rm *.$(obj)