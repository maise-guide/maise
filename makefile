#=============================================================================#
EXE        = maise
CC         = gcc
VER        = maise.2.2.03
CFLAGS     = -O3 -w
GSL_H      = $(shell ./lib/gsl-config --cflags)
GSL_LIB    = $(shell ./lib/gsl-config --libs)
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
SERIAL    ?= 0
LIB       += ./lib/libsymspg.a -lm -lstdc++
INCLUDE    = -I./inc -I./lib/include $(GSL_H) $(GSL_LIB) 
LDFLAGS    = -lm -fopenmp
CPPFLAGS   = -DSERIAL=$(SERIAL) -DVERSION='"$(VER)"' 
CFLAGS    += $(LDFLAGS) $(INCLUDE) $(CPPFLAGS)
DFLAGS     = -MT $@ -MMD -MP -MF 
#=============================================================================#
$(EXE): $(OBJ) 
	@printf "Linking $(_OBJ)...\n";
	@$(CC) $(CFLAGS) -o $@ $^ $(LIB);
	@if [ -e $@ ]; then cp ./$@ ./bin/$@;fi; printf "$@ is ready!\n";

$(ODIR)/%.o: $(SDIR)/%.c $(LDIR)/gsl-config
	@mkdir -p $(DDIR); printf "Compiling $(@F)...\n";
	@$(CC) $(DFLAGS) $(DDIR)/$*.tmpd $(CFLAGS) -c -o $@ $<;
	@mv -f $(DDIR)/$*.tmpd $(DDIR)/$*.d && touch $@; 

$(LDIR)/gsl-config: $(LDIR)/libsymspg.a
	@mkdir -p $(LDIR); printf "Getting $(@F)...\n";
	@bash gsl-dep;

$(LDIR)/libsymspg.a: 
	@mkdir -p $(LDIR); printf "Getting $(@F)...\n";
	@bash spg-dep;

$(DDIR)/%.d = ;

.PRECIOUS: $(DDIR)/%.d 
include $(wildcard $(patsubst %,$(DDIR)/%.d,$(basename $(_SRC)))) 
#=============================================================================#
.PHONY: clean clean-lib clean-all

clean-all: clean clean-lib

clean:
	@printf "Cleaning object files, dependencies, and maise executables...";
	@rm -rf examples/parse/maise examples/relax/maise examples/train/maise;
	@rm -rf $(ODIR)/* ./maise ./bin/maise; printf "Done!\n";

clean-lib: 
	@printf "Cleaning lib..."
	@rm -f  $(LDIR)/*.a; 
	@rm -f  $(LDIR)/include/*.h; 
	@rm -f  $(LDIR)/gsl-config;
	@rm -f  $(LDIR)/include/gsl/*.h; 
	@rm -rf ext-dep;
	@printf "Done!\n"
