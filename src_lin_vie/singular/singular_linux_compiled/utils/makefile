include ../make.inc
INCL_PATH = -I../include

allobj: GL_1D.($obj)

GL_1D.($obj): GL_1D.cpp
	$(CompCXX) -c $(INCL_PATH) GL_1D.cpp

clean:
	rm *.$(obj)