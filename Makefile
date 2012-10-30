
CC = mpixlc_r -qarch=qp -qtune=qp -qsmp=omp:auto:opt
#CC = mpixlc_r -qarch=qp -qtune=qp -qsmp=omp:auto:opt -q64 -g
#CC = mpixlc_r -qarch=qp -qtune=qp -pg -qnostaticlink -g

#CC = /bgsys/drivers/ppcfloor/comm/gcc/bin/mpicxx -fopenmp -DNDEBUG -Wall -O6

LIBS  = -lm

#LIBS += -L /bgsys/linux/RHEL6.2_V1R1M2-12/opt/ibmmath/essl/5.1/lib64/libesslsmpbg.a -lesslsmpbg.a
#-L /bgsys/linux/RHEL6.2_V1R1M2-12/opt/ibmmath/lib64/libesslsmpbg.a


#LIBS += -L/bgsys/ibmhpc/ppedev.hpct/lib64/ -lxlsmp_pomp 
#LIBS += -L/bgsys/ibmhpc/ppedev.hpct/lib64/ -lpomprof_probe 

#LIBS += -pg

#LIBS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lhpc_r 
#LIBS += -L/bgsys/drivers/ppcfloor/bgpm/lib/ -lbgpm 

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
