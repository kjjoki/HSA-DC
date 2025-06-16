# Makefile for the HSA-DC method

FF = gfortran
FFLAGS = -fbounds-check -Wall -pg
OPEN =  

all: execution

execution: HSA-DC_method.o 
	$(FF) -o execution $(FFLAGS) $(OPEN) HSA-DC_method.o

HSA-DC_method.o: HSA-DC_method.for
	$(FF) -c $(FFLAGS) $(OPEN) HSA-DC_method.for 

clean:	
	rm execution HSA-DC_method.o 
	echo Clean done
