EXEC = show
CC = h5pcc
OBJS = main.o parameterSetting.o findparam.o boundary.o saveFile.o fieldSolve.o loadLaser.o fieldShareX.o fieldShareY.o saveDumpHDF.o restoreDumpHDF.o clean.o movingDomain.o loadPlasma.o loadMovingPlasma.o removeEdge.o rearrangeParticles.o interpolation.o particlePush.o updateCurrent.o particleShareX.o particleShareY.o saveDensityHDF.o savePMap.o pml.o saveParticleHDF.o saveFieldHDF.o 
#fieldShareZ.oparticleShareZ.o  clean.o loadPlasma_crystal.o densityShare.o track.o saveRamanHDF.o saveDumpHDF.o fieldTransfer.o

INCL = constants.h laser.h mesh.h particle.h plasma.h
LIBS = -lm -lhdf5
$(EXEC):$(OBJS)
#	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(EXEC)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
