#=============================================================================#
EXE        = maise
CC         = gcc
VER        = 2.0
IPATH     := ~/bin/isotropy
CFLAGS     = -O3 -Wall 
#=============================================================================#
#                 !!! DO NOT CHANGE BELOW THIS LINE!!!                        # 
#=============================================================================#
SDIR       = ./src
ODIR       = ./obj
LDIR       = ./lib
ADIR       = $(LDIR)/spgr
DDIR      := $(ODIR)/dep
#=============================================================================#
_SRC      := $(shell cd $(SDIR); ls *.c; cd ../)
_OBJ       = $(subst .c,.o,$(_SRC))
_LIB       = libgsl.a libgslcblas.a libspgr.a
OBJ        = $(patsubst %,$(ODIR)/%,$(_OBJ))
LIB        = $(patsubst %,$(LDIR)/%,$(_LIB)) 
#=============================================================================#
SERIAL    ?= 0
LIB       += -lstdc++
INCLUDE    = -I./inc -I./lib
LDFLAGS    = -lm -fopenmp
CPPFLAGS   = -DISOPATH='"$(IPATH)"' -DSERIAL=$(SERIAL) -DVERSION='"$(VER)"'
CFLAGS    += $(LDFLAGS) $(INCLUDE) $(CPPFLAGS)
DFLAGS     = -MT $@ -MMD -MP -MF 
#=============================================================================#
$(EXE): $(LDIR)/libspgr.a $(OBJ) 
	@printf "Linking $(_OBJ) and $(_LIB)...\n";
	@$(CC) $(CFLAGS) -o $@ $^ $(LIB);
	@if [ -e $@ ]; then cp ./$@ ./bin/$@;fi; printf "$@ is ready!\n";

$(LDIR)/libspgr.a: $(ADIR)/spgr.c
	@mkdir -p $(ADIR); printf "Compiling $(@F)...\n";
	@$(CC) $(DFLAGS) $(ADIR)/spgr.tmpd  $(CFLAGS) -c -o $(ADIR)/spgr.o $<;
	@ar rcs $@ $(ADIR)/spgr.o; printf "$(@F) compiled successfully!\n";

$(ODIR)/%.o: $(SDIR)/%.c
	@mkdir -p $(DDIR); printf "Compiling $(@F)...\n";
	@$(CC) $(DFLAGS) $(DDIR)/$*.tmpd $(CFLAGS) -c -o $@ $<;
	@mv -f $(DDIR)/$*.tmpd $(DDIR)/$*.d && touch $@; 

$(DDIR)/%.d = ;

.PRECIOUS: $(DDIR)/%.d 
include $(wildcard $(patsubst %,$(DDIR)/%.d,$(basename $(_SRC)))) 
#=============================================================================#
.PHONY: clean clean-spgr clean-all 

clean-all: clean clean-spgr

clean-spgr: 
	@printf "Cleaning lib/spgr/spgr.o and lib/libspgr.a..."
	@rm -f  $(ADIR)/spgr.o $(LDIR)/libspgr.a; printf "Done!\n"
clean:
	@printf "Cleaning object files, dependencies, and maise executables...";
	@rm -rf $(ODIR)/* ./maise ./bin/maise; printf "Done!\n";
