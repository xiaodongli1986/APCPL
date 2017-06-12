#gfortran -o LSS_test LSS_main_LSS_ximu.f90 LSS_main_program.f90

#gfortran -o LSS_MCMC LSS_main_LSS_ximu.f90 LSS_main_ChisqLike.f90 $delm


#lib = lib/libde.a
EXE1 = AP_test
EXE2 = AP_MCMC
EXE3 = AP_grid

OBJS = AP_tools.o

F90C = gfortran  
F90FLAGS = 

default: $(EXE1) $(EXE2) $(EXE3)

AP_test.o : $(OBJS) 
AP_MCMC.o : $(OBJS) 
AP_test: $(OBJS) AP_test.o
AP_MCMC: $(OBJS) AP_MCMC.o


$(EXE1): $(OBJS) AP_test.o
	$(F90C)  -o $(EXE1) $(OBJS) AP_test.o 
$(EXE2): $(OBJS) AP_MCMC.o
	$(F90C)  -o $(EXE2) $(OBJS) AP_MCMC.o -I/home/xiaodongli/software/delib/mods -lde
$(EXE3): $(OBJS) AP_grid.o
	$(F90C)  -o $(EXE3) $(OBJS) AP_grid.o -I/home/xiaodongli/software/delib/mods -lde

%.o: %.f90
	$(F90C) -c $*.f90 -lde -I/home/xiaodongli/software/delib/mods 
clean :
	rm -rf *.o *.mod *.a *.so  $(EXE1) $(EXE2) $(EXE3)
	
