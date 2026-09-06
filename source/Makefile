COMPIL = gfortran
# Mode le plus sur (voir toute erreur et warning) A FAIRE AU DEBUT POUR ETRE SUR
# CFLAGS = -g -O0 -fcheck=all -fbounds-check -Wall -fbacktrace
CFLAGS = -g -O0 -fcheck=all -fbounds-check -fbacktrace -fstack-check
# Mode Performance (bcp plus rapide)
# CFLAGS = -O3
# Meilleur compromis 
#CFLAGS = -O2
LFLAGS = -o
# LFLAGS = -O2 -L$(HOME)
COPTS = $(CFLAGS)
LOPTS = $(LFLAGS)
OBJS = Module_parametres.o Module_sortie.o Module_lecture.o Module_riemann.o Module_preprocess.o Module_solve.o main.o
PROG = solve

.SUFFIXES: .f90 .o
SUFF = .f90.o

%.o: %.f90
	$(COMPIL) -c $(CFLAGS) $< -I.

$(PROG): $(OBJS)
	@echo "Compilation terminee"
	@echo "Creation de l'executable $(PROG)"
	$(COMPIL) $(LOPTS) $(PROG) $(OBJS)

clean:
	@echo "Nettoyage"
	@rm -f core *.o *.mod $(PROG)




