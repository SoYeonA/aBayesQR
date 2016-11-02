
CFLAGS  =	-O3 -DHAVE_TRUNC -fPIC -std=gnu++0x

EXOBJS	=	aBayesQR.o


all:	aBayesQR

clean:
	rm $(EXOBJS) aBayesQR

aBayesQR:	$(EXOBJS)
	g++ $(CFLAGS)  -o $@  $(EXOBJS) 
aBayesQR.o:	aBayesQR.cpp 
	g++   $(CFLAGS)  -c -o $@  $<
