
FC = gfortran
CC = gcc
CFLAGS = -O2 -g
FFLAGS = -O2 -g
FPPFLAGS = -I`nf-config --includedir`
CPPFLAGS = -I`nc-config --includedir`
FLIBS = `nf-config --flibs`
CLIBS = `nc-config --libs`

all:: make_emis make_newcur odcmpr2nc

odcmpr2nc: odcmpr2nc.f90
	$(FC) $(FFLAGS) $(FPPFLAGS) -o $@ $< $(FLIBS)

make_newcur: make_newcur.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(CLIBS)

make_emis: make_emis.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(CLIBS)

clean:
	rm -f odcmpr2nc make_newcur make_emis
