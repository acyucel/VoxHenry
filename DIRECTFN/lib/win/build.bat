@ECHO OFF
cd ..\..\src
python create_makefile.py
nmake -f Makefile_win
DEL *.obj
cd ..\examples\Full_Paper_Examples\c++
python create_makefile.py
nmake -f Makefile_win
cd ..\..\EuCap2017_Examples\c++
python create_makefile.py
nmake -f Makefile_win
cd ..\..\..\lib\win
