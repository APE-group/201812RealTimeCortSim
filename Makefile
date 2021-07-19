# FILE: Makefile
# DESCRIPTION: 
# Makefile for the DPSNN and DiSyGen project - MPI version
# (Distributed Polychronous Spiking Neural Network)
# AUTHOR: Pier Stanislao Paolucci - Roma - Italy (2011-...)
# AUTHOR: Elena Pastorelli - Roma - Italy (2014-....)
# AUTHOR: ...
# AUTHOR: plus other members of INFN Lab, Roma, Italy

CC	 =	mpic++
GCC	 =	g++

#makeOverallTopology:=
# 0 old style explicit stencils and connectivity
# 1 homogeneous connectivity,
# 2 homogeneous conn, but local inhibition
# 3 ...
CXXLIFCAFLAGS = -g -O3 -Wall -fno-omit-frame-pointer -fno-PIC \
	-D miniLocN=32000 \
	-D makeMaxGlobH=256 \
	-D makeStencil_Max=96 \
	-D makeActiveLTP \
	-D makeActiveLTD \
	-D miniM=7500 \
	-D miniSimSpikes=6000 \
	-D makeMaxCMFract=1 \
	-D LIFCAneuron \
	-D makeInsertBarrierBeforeSendRec \
	-mcmodel=large
#	-D makeOverallTopology=1 \
#	-D makeMaxHomogeneousSide=8 \
#	-D stencilDim=8 \
#	-D PerseoDyn 
#	-D DelayUniform
#	-D fixedNumSyn
#	-D periodicBoundaryConditions
#	-D makeActiveLTP
#	-D makeActiveLTD 
#	-D makeInitExcitWeight=0.5
#	-D makeInitInhibWeight=1.5
#	-D makeInsertBarrierAfterSendRec

#CXXFLAGS = -O3 -save-temps -fverbose-asm
#CXXFLAGS = -O0 -g

HFILES = localNetProcess.h dal.h DPSNN_environmentSelection.h		\
	DPSNN_parameters.h DPSNN_dataStructDims.h DPSNN_synapse.h	\
	DPSNN_neuron.h DPSNN_spike.h DPSNN_messagePassing.h		\
	DPSNN_localNet.h DPSNN_random.h DPSNN_stopWatch.h DPSNN_stat.h	\
	DPSNN_chrono.h DPSNN_memMeasure.h DPSNN_getParameters.h \
	DPSNN_messageMeasure.h erflib.h \
	DPSNN_LIFCAconnectome.h randdev.h

SNN_HFILES = snn_singleNet.h\
	snn_dataStructDims.h\
	DPSNN_parameters.h\
	DPSNN_stat.h

CFILES	= dalEmulationUnderMPI.c\
	localNetProcess.c\
	DPSNN_synapse.c\
	DPSNN_neuron_init.c\
	DPSNN_neuron_sim.c\
	DPSNN_spike.c\
	DPSNN_messagePassing.c\
	DPSNN_localNet_init.c\
	DPSNN_synListSort.c\
	DPSNN_localNet_sim.c\
	DPSNN_mainUnderMPI.c\
	DPSNN_random.c\
	DPSNN_stopWatch.c\
	DPSNN_stat.c\
	DPSNN_chrono.c\
	DPSNN_memMeasure.c\
	DPSNN_getParameters.c\
	DPSNN_corticalModuleGen.cpp\
	DPSNN_connectome.cpp\
	DPSNN_localNet_chrono.cpp\
	DPSNN_messageMeasure.c \
	erflib.c \
	DPSNN_LIFCAcorticalModuleGen.cpp \
	DPSNN_LIFCAconnectome.cpp \
	DPSNN_writeIniFiles.c \
	randdev.c

all:    mpi-LIFCA-DPSNN

clean:  
	rm ./LIFCA-DPSNN.out *.o

LIFCAclean:  
	rm ./LIFCA-DPSNN.out \
	*LIFCA*.o

mpi-LIFCA-DPSNN: LIFCA-DPSNN_mainUnderMPI.o \
	LIFCA-DPSNN_synapse.o \
	LIFCA-DPSNN_neuron_init.o \
	LIFCA-DPSNN_neuron_sim.o \
	LIFCA-DPSNN_spike.o \
	LIFCA-DPSNN_messagePassing.o \
	LIFCA-DPSNN_localNet_init.o \
	LIFCA-DPSNN_localNet_sim.o \
	LIFCA-DPSNN_synListSort.o \
	LIFCA-localNetProcess.o \
	LIFCA-DPSNN_getParameters.o \
	LIFCA-dalEmulationUnderMPI.o \
	LIFCA-DPSNN_random.o \
	LIFCA-DPSNN_stopWatch.o \
	LIFCA-DPSNN_stat.o \
	LIFCA-DPSNN_chrono.o \
	LIFCA-DPSNN_memMeasure.o \
	LIFCA-DPSNN_LIFCAcorticalModuleGen.o \
	LIFCA-DPSNN_LIFCAconnectome.o \
	LIFCA-DPSNN_writeIniFiles.o \
	LIFCA-DPSNN_localNet_chrono.o \
	LIFCA-DPSNN_messageMeasure.o \
	LIFCA-erflib.o \
	LIFCA-randdev.o
	$(CC) $(CXXLIFCAFLAGS) LIFCA-DPSNN_mainUnderMPI.o \
	LIFCA-DPSNN_synapse.o \
	LIFCA-DPSNN_neuron_init.o \
	LIFCA-DPSNN_neuron_sim.o \
	LIFCA-DPSNN_spike.o \
	LIFCA-DPSNN_messagePassing.o \
	LIFCA-DPSNN_localNet_init.o \
	LIFCA-DPSNN_localNet_sim.o \
	LIFCA-DPSNN_synListSort.o \
	LIFCA-localNetProcess.o \
	LIFCA-DPSNN_getParameters.o \
	LIFCA-dalEmulationUnderMPI.o \
	LIFCA-DPSNN_random.o \
	LIFCA-DPSNN_stopWatch.o \
	LIFCA-DPSNN_stat.o \
	LIFCA-DPSNN_chrono.o \
	LIFCA-DPSNN_memMeasure.o \
	LIFCA-DPSNN_LIFCAcorticalModuleGen.o \
	LIFCA-DPSNN_LIFCAconnectome.o \
	LIFCA-DPSNN_writeIniFiles.o \
	LIFCA-DPSNN_localNet_chrono.o \
	LIFCA-DPSNN_messageMeasure.o \
	LIFCA-erflib.o \
	LIFCA-randdev.o \
	-o ./LIFCA-DPSNN.out

LIFCA-dalEmulationUnderMPI.o: $(HFILES) dalEmulationUnderMPI.c 
	$(CC) $(CXXLIFCAFLAGS) dalEmulationUnderMPI.c -c -o LIFCA-dalEmulationUnderMPI.o

LIFCA-DPSNN_synapse.o: $(HFILES) DPSNN_synapse.c 
	$(CC) $(CXXLIFCAFLAGS) DPSNN_synapse.c -c -o LIFCA-DPSNN_synapse.o

LIFCA-DPSNN_neuron_init.o: $(HFILES) DPSNN_neuron_init.c 
	$(CC) $(CXXLIFCAFLAGS) DPSNN_neuron_init.c -c -o LIFCA-DPSNN_neuron_init.o

LIFCA-DPSNN_neuron_sim.o: $(HFILES) DPSNN_neuron_sim.c 
	$(CC) $(CXXLIFCAFLAGS) DPSNN_neuron_sim.c -c -o LIFCA-DPSNN_neuron_sim.o

LIFCA-DPSNN_spike.o: $(HFILES) DPSNN_spike.c 
	$(CC) $(CXXLIFCAFLAGS) DPSNN_spike.c -c -o LIFCA-DPSNN_spike.o

LIFCA-DPSNN_messagePassing.o: $(HFILES) DPSNN_messagePassing.c
	$(CC) $(CXXLIFCAFLAGS) DPSNN_messagePassing.c -c -o LIFCA-DPSNN_messagePassing.o

LIFCA-DPSNN_localNet_init.o: $(HFILES) DPSNN_localNet_init.c
	$(CC) $(CXXLIFCAFLAGS)  DPSNN_localNet_init.c \
	-c -o LIFCA-DPSNN_localNet_init.o

LIFCA-DPSNN_synListSort.o: $(HFILES) DPSNN_synListSort.c
	$(CC) $(CXXLIFCAFLAGS) DPSNN_synListSort.c -c -o LIFCA-DPSNN_synListSort.o

LIFCA-DPSNN_localNet_sim.o: $(HFILES) DPSNN_localNet_sim.c
	$(CC) $(CXXLIFCAFLAGS) DPSNN_localNet_sim.c -c -o LIFCA-DPSNN_localNet_sim.o

LIFCA-localNetProcess.o: $(HFILES) localNetProcess.c
	$(CC) $(CXXLIFCAFLAGS) localNetProcess.c -c -o LIFCA-localNetProcess.o

LIFCA-DPSNN_getParameters.o: $(HFILES) DPSNN_getParameters.c
	$(CC) $(CXXLIFCAFLAGS) DPSNN_getParameters.c -c -o LIFCA-DPSNN_getParameters.o

LIFCA-DPSNN_mainUnderMPI.o: $(HFILES) DPSNN_mainUnderMPI.c
	$(CC) $(CXXLIFCAFLAGS) DPSNN_mainUnderMPI.c -c -o LIFCA-DPSNN_mainUnderMPI.o

LIFCA-DPSNN_random.o: $(HFILES) DPSNN_random.c
	$(CC) $(CXXLIFCAFLAGS) DPSNN_random.c \
	-c -o LIFCA-DPSNN_random.o

LIFCA-DPSNN_stopWatch.o: $(HFILES) DPSNN_stopWatch.c
	$(CC) $(CXXLIFCAFLAGS) DPSNN_stopWatch.c -c -o LIFCA-DPSNN_stopWatch.o

LIFCA-DPSNN_stat.o: $(HFILES) DPSNN_stat.c
	$(CC) $(CXXLIFCAFLAGS) DPSNN_stat.c -c -o LIFCA-DPSNN_stat.o

LIFCA-DPSNN_chrono.o: $(HFILES) DPSNN_chrono.c
	$(CC) $(CXXLIFCAFLAGS) DPSNN_chrono.c -c -o LIFCA-DPSNN_chrono.o

LIFCA-DPSNN_memMeasure.o: $(HFILES) DPSNN_memMeasure.c
	$(CC) $(CXXLIFCAFLAGS) DPSNN_memMeasure.c -c -o LIFCA-DPSNN_memMeasure.o

LIFCA-DPSNN_LIFCAcorticalModuleGen.o: $(HFILES) DPSNN_LIFCAcorticalModuleGen.cpp
	$(CC) $(CXXLIFCAFLAGS) DPSNN_LIFCAcorticalModuleGen.cpp \
	-c -o LIFCA-DPSNN_LIFCAcorticalModuleGen.o

LIFCA-DPSNN_LIFCAconnectome.o: $(HFILES) DPSNN_LIFCAconnectome.cpp
	$(CC) $(CXXLIFCAFLAGS) DPSNN_LIFCAconnectome.cpp \
	-c -o LIFCA-DPSNN_LIFCAconnectome.o

LIFCA-DPSNN_writeIniFiles.o: $(HFILES) DPSNN_writeIniFiles.c
	$(CC) $(CXXLIFCAFLAGS) DPSNN_writeIniFiles.c \
	-c -o LIFCA-DPSNN_writeIniFiles.o

LIFCA-DPSNN_localNet_chrono.o: $(HFILES) DPSNN_localNet_chrono.cpp
	$(CC) $(CXXLIFCAFLAGS) DPSNN_localNet_chrono.cpp \
	-c -o LIFCA-DPSNN_localNet_chrono.o

LIFCA-DPSNN_messageMeasure.o: $(HFILES) DPSNN_messageMeasure.c
	$(CC) $(CXXLIFCAFLAGS) DPSNN_messageMeasure.c -c -o LIFCA-DPSNN_messageMeasure.o

LIFCA-erflib.o: $(HFILES) erflib.c
	$(CC) $(CXXLIFCAFLAGS) erflib.c -c -o LIFCA-erflib.o

LIFCA-randdev.o: $(HFILES) randdev.c
	$(CC) $(CXXLIFCAFLAGS) randdev.c -c -o LIFCA-randdev.o

snn_singleNet.o: $(SNN_HFILES) snn_singleNet.cpp
	g++ $(CXXFLAGS) snn_singleNet.cpp -c -o snn_singleNet.o

snn_main.o : snn_main.cpp $(SNN_HFILES)
	g++ $(CXXFLAGS) snn_main.cpp -c -o snn_main.o

snn : snn_main.o snn_singleNet.o DPSNN_stat.o
	g++ $(CXXFLAGS) snn_main.o snn_singleNet.o DPSNN_stat.o -o ./s/snn.out


