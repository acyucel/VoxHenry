
import os
makefile = open('Makefile_win','w')
makefile.write('ROOT_DIR = ../../\n\n')
makefile.write('include ../../settings/Makefile.nmake_win\n\n')
makefile.write('INCL_PATH = -I../../include\n')


makefile.write('LIB_PATH = ../../lib/win\n')
makefile.write('LOCAL_LIB_NAME = $(DIRECTFN_LIB_NAME).$(SO)\n\n\n')
makefile.write('all: allobj allexe\n\n')

makefile.write('allobj: ')

names = [scr.split('.')[0] for scr in os.listdir('./') if scr.endswith('.cpp')]
all_obj = [name+'.$(obj) \\' + '\n\t\t' for name in names]
all_obj[-1] = all_obj[-1].split('.')[0]+'.$(obj)\n\t'
for object in all_obj:
	makefile.write(object)

makefile.write('\n')
makefile.write('allexe: ')
all_exe = [name+'.$(exe) \\' + '\n\t\t' for name in names]
all_exe[-1] = all_exe[-1].split('.')[0]+'.$(exe)\n\t'
for exe in all_exe:
	makefile.write(exe)

makefile.write('\n')
for name in names:
    line = '%s.$(obj): %s.cpp\n\t$(CompCXX) -c $(FlagsCXX) $(INCL_PATH) %s.cpp\n\n' % (name,name,name)
    makefile.write(line)

makefile.write('\n')	
	
for name in names:
    line = '%s.$(exe): %s.obj\n\t link $(FlagsAR) /LIBPATH:$(LIB_PATH) %s.obj $(LOCAL_LIB_NAME)\n\n' % (name,name,name)
    makefile.write(line)



makefile.write('clean:\n\t del *.obj *.exe\n')

makefile.close()






