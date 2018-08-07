
# MF_PREFIX is used in all subdirs Makefiles with src
MF_PREFIX=../../settings
include  $(MF_PREFIX)/Makefile.in

.PHONY: directfn clean 

directfn:  
	(cd ../../src && make all && make clean)

clean:
	rm lib*.$(SO)

# End of the make-lib file



