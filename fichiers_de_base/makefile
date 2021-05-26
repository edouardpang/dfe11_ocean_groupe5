# Choix des clefs cpp

#CPPFLAGS      = -DrealWorld -DbetaPlane -DdynZonalFlow
CPPFLAGS      = -DdynGill

# ensemble des options cpp disponibles
# Geometrie et plan f/beta : -DrealWorld -Dsill -DbetaPlane 
# Forcage par le vent : -DwindStress
# Etat initial : -DdynWarmPool -DdynKelvinWave -DdynRossbyWave -DdynZonalFlow -DdynGill -DdynFreeConfig -Dnonlinear -DdynOBC
# Restart : -DdynRestart
# Diagnostics integraux (energie, enstrophie, ...) : -DdiaTrends
# Schémas numériques : -Dnonlinear -Denstrophy


FC = gfortran
OBJS = main.o init.o step.o diag.o onc_open.o onc_apfltmsk.o onc_defstdgrd.o onc_putstdvar.o onc_newstep.o
SRCS = main.F90 init.F90 step.F90 diag.F90 onc_open.F onc_apfltmsk.F onc_defstdgrd.F onc_putstdvar.F onc_newstep.F

HOMELIB     = $(shell echo ~mortier/lib/)
INCLUDES    = -I$(HOMELIB)
#INCLUDES    = $(HOMEDIR)
#FFLAGS      =  -fimplicit-none -DopaCDF_debug $(INCLUDES) -DNDEBUG
FFLAGS      =  -fimplicit-none $(INCLUDES) -DNDEBUG
FFLAGS8     = $(CPPFLAGS) $(FFLAGS) -fdefault-real-8
LFLAGS    = -L $(HOMELIB) -l cdfpack4_3 -L /usr/lib64/ -l netcdf -l netcdff

main:	$(OBJS)
	$(FC)  $(FFLAGS) -o main.exe  $(OBJS)  $(LFLAGS) 
main.o:	main.F90 parameter.h common.h makefile
	$(FC) $(FFLAGS8) -c main.F90
init.o:	init.F90 parameter.h common.h inicoriolis.h  inidepth.h  inidyn.h  inigrid.h  inirelax.h makefile
	$(FC) $(FFLAGS8) -c init.F90
step.o:	step.F90 parameter.h common.h  makefile
	$(FC) $(FFLAGS8) -c step.F90
diag.o:	diag.F90 parameter.h common.h  makefile
	$(FC) $(FFLAGS8) -c diag.F90
onc_putstdvar.o: onc_putstdvar.F parameter.h common.h  
	$(FC) $(FFLAGS) -c -ffixed-line-length-132 onc_putstdvar.F
onc_apfltmsk.o:	onc_apfltmsk.F parameter.h common.h 
	$(FC) $(FFLAGS) -c -ffixed-line-length-132 onc_apfltmsk.F
onc_defstdgrd.o: onc_defstdgrd.F parameter.h common.h
	$(FC) $(FFLAGS) -c -ffixed-line-length-132 onc_defstdgrd.F
onc_newstep.o:	onc_newstep.F parameter.h common.h
	$(FC) $(FFLAGS) -c -ffixed-line-length-132 onc_newstep.F
onc_open.o:	onc_open.F parameter.h common.h
	$(FC) $(FFLAGS) -c -ffixed-line-length-132 onc_open.F

zip:	
	zip mf210_F90.v0.3.zip       makefile signals.conf epic.key namelist.txt *.sh *.h *.m *.F90 *.F

tar:	
	tar -cvf mf210_F90.v0.3.tar  makefile signals.conf epic.key namelist.txt *.sh *.h *.m *.F90 *.F

echo: 
	@echo $(HOMELIB)

clean: 
	rm -f *.o main.exe channel.ascii mf210_F90.*
