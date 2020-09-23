EXE    = fort2maise
FC     = gfortran -fcheck=all -fbacktrace 

M_PATH = ../../

LIB    = $(M_PATH)/lib/libmaise.a -fopenmp

./bin/$(EXE): $(EXE).f90 maise.o maisedef.o
	@$(FC) -o $@ $^ $(LIB);rm ./*.o;

maise.o:  maise.f90 maisedef.o
	@$(FC) -c maise.f90;

maisedef.o: maisedef.f90 
	@$(FC) -c maisedef.f90;

clean:
	@rm -rf ./bin/$(EXE){,.o} ./*.{o,mod} 2>&1;
	@rm -rf ./maisedef.{o.mod} 2>&1;
