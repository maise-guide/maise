#=============================================================================#
.RECIPEPREFIX := >
EXE        = maise
MLIB       = libmaise.a
CC         = gcc
VER        = maise.3.0.02
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
>@ld -o maise.o -r $(shell ls ./obj/*.o | grep -v main.o ) $(GSL_H) ./lib/libgsl.a ./lib/libgslcblas.a ./lib/libsymspg.a;
>@ar rcs $@ maise.o; ranlib $@;
>@if [ -e $@ ]; then mv ./$@ ./lib/$@;fi; rm -rf ./maise.o;

$(EXE): $(OBJ)
>@printf "Linking $(_OBJ)...\n";
>@$(CC) $(CFLAGS) -o $@ $^ $(LIB);
>@if [ -e $@ ]; then mkdir -p ./bin; cp ./$@ ./bin/$@;fi; printf "\033[0;32m$@ is ready!\n\033[0m";

$(ODIR)/%.o: $(SDIR)/%.c $(LDIR)/gsl-config
>@mkdir -p $(DDIR); printf "Compiling $(@F)...\n";
>@$(CC) $(DFLAGS) $(DDIR)/$*.tmpd $(CFLAGS) -c -o $@ $<;
>@mv -f $(DDIR)/$*.tmpd $(DDIR)/$*.d && touch $@;

$(LDIR)/gsl-config: $(LDIR)/libsymspg.a
>@if [ ! -e $@ ]; then mkdir -p $(LDIR); printf "Getting $(@F)...\n"; bash dep-gsl || { rc=$$?; printf "\033[0;31mERROR: dep-gsl failed\033[0m\n"; if [ -f dep-ext/gsl.log ]; then tail -80 dep-ext/gsl.log; fi; exit $$rc; }; fi;
>@if [ ! -x "$(LDIR)/gsl-config" ]; then \
>    printf "\033[0;31mERROR: dep-gsl did not create executable $(LDIR)/gsl-config\033[0m\n"; \
>    if [ -f dep-ext/gsl.log ]; then \
>        printf "\nLast 80 lines of dep-ext/gsl.log:\n"; \
>        tail -80 dep-ext/gsl.log; \
>    fi; \
>    exit 1; \
>fi;

$(LDIR)/libsymspg.a:
>@mkdir -p $(LDIR) $(LDIR)/include; printf "Getting $(@F)...\n";
>@bash dep-spg || { rc=$$?; printf "\033[0;31mERROR: dep-spg failed\033[0m\n"; if [ -f dep-ext/spglib.log ]; then printf "\nLast 80 lines of dep-ext/spglib.log:\n"; tail -80 dep-ext/spglib.log; fi; exit $$rc; }
>@if [ ! -s "$(LDIR)/libsymspg.a" ]; then \
>    printf "\033[0;31mERROR: dep-spg did not create $(LDIR)/libsymspg.a\033[0m\n"; \
>    printf "This is a SPGLIB build/install failure, not a final MAISE link failure.\n"; \
>    if [ -f dep-ext/spglib.log ]; then \
>        printf "\nLast 80 lines of dep-ext/spglib.log:\n"; \
>        tail -80 dep-ext/spglib.log; \
>    fi; \
>    printf "\nSPGLIB files found under dep-ext:\n"; \
>    find dep-ext \( -name 'libsymspg.a' -o -name 'libsymspg.so*' -o -name 'spglib.h' \) -print 2>/dev/null || true; \
>    exit 1; \
>fi;
>@if [ ! -s "$(LDIR)/include/spglib.h" ]; then \
>    printf "\033[0;31mERROR: dep-spg did not create $(LDIR)/include/spglib.h\033[0m\n"; \
>    printf "SPGLIB header files found under dep-ext:\n"; \
>    find dep-ext -name 'spglib.h' -print 2>/dev/null || true; \
>    exit 1; \
>fi;

$(DDIR)/%.d = ;

.PRECIOUS: $(DDIR)/%.d
include $(wildcard $(patsubst %,$(DDIR)/%.d,$(basename $(_SRC))))
#=============================================================================#
.PHONY: clean clean-lib clean-all

clean-all: clean clean-lib

clean:
>@printf "Cleaning object files, dependencies, and maise executables...";
>@rm -rf test/parse/maise test/relax/maise test/train/maise;
>@rm -rf $(ODIR)/* ./maise ./bin/maise ./lib/libmaise.a ; printf "Done!\n";
>@rm -rf obj;

clean-lib:
>@printf "Cleaning lib..."
>@rm -f  $(LDIR)/*.a;
>@rm -f  $(LDIR)/include/*.h;
>@rm -f  $(LDIR)/gsl-config;
>@rm -f  $(LDIR)/include/gsl/*.h;
>@rm -rf dep-ext;
>@rm -rf lib;
>@printf "Done!\n"
