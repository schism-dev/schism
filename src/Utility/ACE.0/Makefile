all:
	
	make -C elio
	make -C lib
	cp lib/fli/libfpp.a lib/
	cp elio/libelio.a lib/ 
	mkdir -p include
	cp lib/fli/fpfli.h include/
	cp lib/graphics/patterns.h include/
	cp elio/elio.h include/
	cp lib/graphics/hersh.h include/
	cp lib/graphics/special.h include/

	make -C xmgredit5/ 
	cp xmgredit5/xmgredit5 $(HOME)/bin

	make -C  xmvis6/ 
	cp xmvis6/xmvis6 $(HOME)/bin

	cd xmgr5
	rm -f pars.c pars.h 
	cd ..
	make -C xmgr5 pars.c 
#	make -C xmgr5 clean
	make -C xmgr5
	cp xmgr5/xmgr5 $(HOME)/bin

