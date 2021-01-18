#=============================================================================#
EXE        = maise
MLIB       = libmaise.a
CC         = gcc
VER        = maise.2.8.09
CFLAGS     = -O3 -fno-strict-overflow
GSL_H      = $(shell ./lib/gsl-config --cflags 2> /dev/null)
GSL_LIB    = $(shell ./lib/gsl-config --libs   2> /dev/null)
#=============================================================================#
#                 !!! DO NOT CHANGE BELOW THIS LINE!!!                        # 
#=============================================================================#
SDIR       = ./src
ODIR       = ./obj
LDIR       = ./lib
DDIR      := $(ODIR)/dep
#=============================================================================#
_SRC      := $(shell cd $(SDIR); ls *.c; cd ../)
_OBJ       = $(subst .c,.o,$(_SRC))
OBJ        = $(patsubst %,$(ODIR)/%,$(_OBJ))
LIB        = $(GSL_LIB)
#=============================================================================#
LIB       += ./lib/libsymspg.a -lm -lstdc++
INCLUDE    = -I./inc -I./lib/include $(GSL_H) $(GSL_LIB) 
LDFLAGS    = -lm -fopenmp
CPPFLAGS   = -DVERSION='"$(VER)"' 
CFLAGS    += $(LDFLAGS) $(INCLUDE) $(CPPFLAGS)
DFLAGS     = -MT $@ -MMD -MP -MF 
#=============================================================================#
$(MLIB): $(OBJ) $(EXE)
	@ld -o maise.o -r $(shell ls ./obj/*.o | grep -v main.o ) $(GSL_H) ./lib/libgsl.a ./lib/libgslcblas.a ./lib/libsymspg.a;
	@ar rcs $@ maise.o; ranlib $@;
	@if [ -e $@ ]; then mv ./$@ ./lib/$@;fi; rm -rf ./maise.o; 

$(EXE): $(OBJ) 
	@printf "Linking $(_OBJ)...\n";
	@$(CC) $(CFLAGS) -o $@ $^ $(LIB);
	@if [ -e $@ ]; then cp ./$@ ./bin/$@;fi; printf "\033[0;32m$@ is ready!\n\033[0m";

$(ODIR)/%.o: $(SDIR)/%.c $(LDIR)/gsl-config
	@mkdir -p $(DDIR); printf "Compiling $(@F)...\n";
	@$(CC) $(DFLAGS) $(DDIR)/$*.tmpd $(CFLAGS) -c -o $@ $<;
	@mv -f $(DDIR)/$*.tmpd $(DDIR)/$*.d && touch $@; 

$(LDIR)/gsl-config: $(LDIR)/libsymspg.a
	@if [ ! -e $@ ]; then mkdir -p $(LDIR); printf "Getting $(@F)...\n"; bash dep-gsl;fi;

$(LDIR)/libsymspg.a: 
	@mkdir -p $(LDIR); printf "Getting $(@F)...\n";
	@bash dep-spg;

$(DDIR)/%.d = ;

.PRECIOUS: $(DDIR)/%.d 
include $(wildcard $(patsubst %,$(DDIR)/%.d,$(basename $(_SRC)))) 
#=============================================================================#
.PHONY: clean clean-lib clean-all

clean-all: clean clean-lib

clean:
	@printf "Cleaning object files, dependencies, and maise executables...";
	@rm -rf test/parse/maise test/relax/maise test/train/maise;
	@rm -rf $(ODIR)/* ./maise ./bin/maise ./lib/libmaise.a ; printf "Done!\n";
	@rm -rf obj;

clean-lib: 
	@printf "Cleaning lib..."
	@rm -f  $(LDIR)/*.a; 
	@rm -f  $(LDIR)/include/*.h; 
	@rm -f  $(LDIR)/gsl-config;
	@rm -f  $(LDIR)/include/gsl/*.h; 
	@rm -rf dep-ext;
	@rm -rf lib;
	@printf "Done!\n"
