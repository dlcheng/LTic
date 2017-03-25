#---------------------------------------Transfer functions
#OPT += -DBBKS_T_K
OPT += -DEH_T_K
#OPT += -DWD_ST                          # warm dark mattter, sterile neutrino actrually
#----------------------------------------For test
OPT += -DCUT_OFF
#--------------------------------------- Select target computer
#SYSTYPE="dlcheng"
#SYSTYPE="ITSC"
SYSTYPE="Mac"
#--------------------------------------- Adjust settings for target computer

ifeq ($(SYSTYPE),"dlcheng")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall 
GSL_INCL =  -I/home/dalong/Install/gsl/include
GSL_LIBS =  -L/home/dalong/Install/gsl/lib
endif

ifeq ($(SYSTYPE),"Mac")
CC       =   gcc-4.9   
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/usr/local/Cellar/gsl/1.16/lib 
GSL_INCL =  -I/usr/local/Cellar/gsl/1.16/include
endif

ifeq ($(SYSTYPE),"ITSC")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/usr/local/gsl-1.14/lib  -Wl,"-R /usr/local/gsl-1.14/lib" 
GSL_INCL =  -I/usr/local/gsl-1.14/include
endif



OPTIONS =  $(OPTIMIZE) $(OPT) 

EXEC   = LTic

OBJS   = allvars.o eh_tk.o initialization.o local_shape_transfer.o main.o memory_control.o\
         normalization.o output_convert.o warning.o zedolvech.o

INCL   = allvars.h proto.h define.h Makefile


CFLAGS = $(OPTIONS) $(GSL_INCL) -fopenmp


LIBS   = $(GSL_LIBS) -lgsl -lgslcblas -lm -fopenmp

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC) *.gch
