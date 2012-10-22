
CC = mpixlc_r -qarch=qp -qtune=qp # -qsmp=omp:auto:opt

#CC = /bgsys/drivers/ppcfloor/comm/gcc/bin/mpicxx -fopenmp -DNDEBUG -Wall -O6

LIBS =	-lm

OBJSM = model.o dist.o gauss.o
OBJSD =  dysfunction.o gauss.o

all:	model dysfunction

model:	$(OBJSM)
		$(CC) $(CFLAGS) -o model $(OBJSM) $(LIBS)

dysfunction: $(OBJSD)
		$(CC) $(CFLAGS) -o dysfunction $(OBJSD) $(LIBS)

clean:
		/bin/rm -f $(OBJSM) $(OBJSD) *.bak

dist.o: dist.h Makefile
model.o: model.h dist.h Makefile gauss.h

dysfunction.o: model.h dysfunction.h Makefile
gauss.o: gauss.h Makefile
