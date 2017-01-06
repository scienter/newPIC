EXEC = show
CC = mpicc
OBJS = resolChange.o subBoundary.o hdfModule.o particleSave.o clean.o restoreDenPart.o findAverage.o mathUtil.o
INCL = mesh.h
LIBS = -lm -lhdf5 -lgsl -lgslcblas
$(EXEC):$(OBJS)
#	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(EXEC)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
