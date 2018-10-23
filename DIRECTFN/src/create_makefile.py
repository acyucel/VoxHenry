
# coding: utf-8

# In[1]:

import os
makefile = open('Makefile_win','w')


# In[2]:

line1 = 'include ../settings/Makefile.nmake_win'
line2 = 'INCL_PATH = -I../include'
line3 = 'LOCAL_LIB_NAME = $(DIRECTFN_LIB_NAME).$(SO)\n\n'
line4 = 'all: allobj stalib'

# In[3]:

makefile.write("%s\n%s\n%s\n\n%s\n\n" % (line1, line2, line3, line4))


names = [scr.split('.')[0] for scr in os.listdir('./') if scr.endswith('.cpp')]


objects = [name+'.($obj)' for name in names]
string_objects =  ' '.join(objects)


makefile.write('allobj: %s\n\n' %(string_objects))


for name in names:
    line = '%s.($obj): %s.cpp\n\t$(CompCXX)  -c $(FlagsCXX) $(INCL_PATH) %s.cpp\n\n' % (name,name,name)
    makefile.write(line)

command_lib = "stalib:\n\tlib $(FlagsAR) /OUT:../lib/win/$(LOCAL_LIB_NAME) \\\n\t\t"
makefile.write("%s" % (command_lib))

objects = [src.split('.')[0]+'.$(obj) \\' +'\n\t\t' for src in os.listdir('./') if src.endswith('.cpp')]  
objects[-1]  = objects[-1].split('.')[0]+'.$(obj)\n\t\t'
for object in objects:
    makefile.write('%s' % (object))
makefile.write('\n\t')
makefile.write('cd ../lib/win\n\n')

makefile.close()




