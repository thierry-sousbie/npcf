#Choose PC or MAC:
ARCH=PC
#Properties can be set to TRUE or FALSE

DEBUG=FALSE
USE_THREADS=TRUE
#TREE_SHUFFLE improves speed, only set to FALSE if you are memory limited
TREE_SHUFFLE=TRUE
#To use floats instead of doubles (floats are good enough)
SIMPLEPRECISION=TRUE
#To use long int indexes (false unless there are really a lot of particles)
#this is not implemented anyway :)
LONGINT=FALSE
#use only on 32 bit systems
FORCE64=FALSE

# Don't use this, threads alone are better
#USE_OPENMP=FALSE 

CC=gcc
CPP=g++

ifeq ($(DEBUG),TRUE)
ifeq ($(FORCE64),TRUE)
	CFLAGS= -Wall -g -DDEBUG=1 -m64
	LDFLAGS= -Wall -g -m64 
else
	CFLAGS=-Wall -g -DDEBUG=1 
	LDFLAGS=-Wall -g 
endif
else
ifeq ($(FORCE64),TRUE)
	CFLAGS= -O3 -g -m64 #-pg
	LDFLAGS= -O3 -g -m64 #-pg 
else
	CFLAGS= -O3 -g #-pg
	LDFLAGS= -O3 -g #-pg
endif
endif



#########################################################
############ Change the libraries paths here ############
#########################################################

LIBS = -lm 

#generic path
LIBS_PATH = -L/usr/X11R6/lib -L/usr/X11R6/lib64 -L/sw/lib -L/usr/lib -L/usr/lib64 -L/usr/local/lib  -L/usr/local/lib64 -L/home/sousbieth/lib/lib -L/home/sousbie/apps/cfitsio/lib/
INCLUDES = -I/usr/X11R6/include -I/sw/include -I/home/sousbieth/lib/include -I/home/sousbie/apps/cfitsio/include/

ifeq ($(USE_THREADS),TRUE)
LIBS += -lpthread
INCLUDES += 
CFLAGS += -DUSE_THREADS -DPARALLEL
endif

ifeq ($(USE_OPENMP),TRUE)
CFLAGS += -DUSE_OPENMP -fopenmp -DPARALLEL -DUSE_THREADS
LDFLAGS += -fopenmp
endif

#64 bit path (only for 32 bit systems)
ifeq ($(FORCE64),TRUE)
LIBS_PATH += -L/Users/sousbie/apps/lib64/lib 
INCLUDES += -I/Users/sousbie/apps/lib64/include
endif

########################################################
######### Nothing should be changed below this #########
########################################################

ifeq ($(ARCH),MAC)
endif

ifeq ($(ARCH),PC)
endif

ifeq ($(SIMPLEPRECISION),TRUE)
CFLAGS += -DSIMPLEPRECISION
endif

ifeq ($(LONGINT),TRUE)
CFLAGS += -DUSELONGINT
endif

ifeq ($(TREE_SHUFFLE),TRUE)
CFLAGS += -DTREE_SHUFFLE
endif

LIBS += $(LIBS_PATH)
INCLUDES += 

CFLAGS+= $(INCLUDES)
LDFLAGS+= $(LIBS)

NPCF = ./bin/npcf.exe

NPCF_OBJS = npcf.o endian.o gadget_io.o NDfield.o mystring.o bbtree.o \
	correl.o correl_parms.o tools.o 2points.o 3points.o 

ifeq ($(USE_THREADS),TRUE)
NPCF_OBJS += correl_threads.o
endif

all: $(NPCF)
clean: 
	@rm -vf  *~ core* "#"*"#" ./bin/*.exe *.o

$(NPCF): $(NPCF_OBJS)

%.o: %.c %.h
	$(CC) -o $@ -c $< $(CFLAGS)

%.o: %.cpp %.h
	$(CPP) -o $@ -c $< $(CFLAGS)

%.exe:
	$(CC) -o $@ $^ $(LDFLAGS)


