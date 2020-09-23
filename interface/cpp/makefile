EXE    = cpp2maise
CC     = gcc 

M_PATH = ../../

LIB    = $(M_PATH)/lib/libmaise.a -lm -fopenmp;
INC    = -I$(M_PATH)/inc -I$(M_PATH)/lib/include

./bin/$(EXE): $(EXE).cpp
	@$(CC) -o $@ $^ $(INC) $(LIB) 

clean:
	@rm -rf ./bin/$(EXE)
