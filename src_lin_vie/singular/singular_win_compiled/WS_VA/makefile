include ../../make.inc
INCL_PATH = -I../../include

allobj: quadric_ws_va.($obj) rho_limits_ws_va.($obj) theta_p_limits_ws_va.($obj) theta_q_limits_ws_va.($obj)

quadric_ws_va.($obj): quadric_ws_va.cpp
	$(CompCXX) -c $(INCL_PATH) quadric_ws_va.cpp

rho_limits_ws_va.($obj): rho_limits_ws_va.cpp
	$(CompCXX) -c $(INCL_PATH) rho_limits_ws_va.cpp

theta_p_limits_ws_va.($obj): theta_p_limits_ws_va.cpp
	$(CompCXX) -c $(INCL_PATH) theta_p_limits_ws_va.cpp

theta_q_limits_ws_va.($obj): theta_q_limits_ws_va.cpp
	$(CompCXX) -c $(INCL_PATH) theta_q_limits_ws_va.cpp

clean:
	rm *.$(obj)