# Time-stamp: <2014/01/13 13:04 baron>

F90C=gfortran
FFLAGS =-O3 -m64
LDFLAGS=

all: advdif blast cholesky cluster conjug drunkard findata finelem2 fluidyn gaussor heatcrni heatelem heatex heatri hotplate leapdf mac magelec markov3 mcplate metropol oscillat pic polywalk radbar reactor satellit shootemp trojans twodwave wave

advdif: advdif.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

blast: blast.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

cholesky: cholesky.o iccg.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

cluster: cluster.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

conjug: conjug.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

drunkard: drunkard.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

findata: findata.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

finelem2: finelem2.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

fluidyn: fluidyn.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

gaussor: gaussor.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

heatcrni: heatcrni.o tridiag.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

heatelem: heatelem.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

heatex: heatex.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

heatri: heatri.o tridiag.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

hotplate: hotplate.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

leapdf: leapdf.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

mac: mac.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

magelec: magelec.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

markov3: markov3.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

mcplate: mcplate.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

metropol: metropol.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

oscillat: oscillat.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

pic: pic.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

polywalk: polywalk.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

radbar: radbar.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

reactor: reactor.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

satellit: satellit.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

shootemp: shootemp.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

trojans: trojans.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

twodwave: twodwave.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

wave: wave.o
	$(F90C) $(FFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm *.o
	for i in *.f;\
		do \
		rm $${i%%.f}; \
		done

# compile
%.o: %.f90
	$(F90C)	 -c $(FFLAGS) $^ -o $*.o
