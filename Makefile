CC     = cc -O2
CFLAGS = -ansi -pedantic -Wall
LIBS   = -L$(HOME)/lib -lbiop -lgen -lgdbm
INC    = -I$(HOME)/include
OFILES = nr.o
DEFS   = 

nr : $(OFILES)
	$(CC) -o $@ $(OFILES) $(LIBS)

.c.o :
	$(CC) $(DEFS) $(CFLAGS) -c $(INC) $<

clean :
	\rm -f $(OFILES)
