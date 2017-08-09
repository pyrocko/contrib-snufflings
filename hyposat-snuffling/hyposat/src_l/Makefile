##############################################################
#                                                            #
# Makefile for the program   hyposat                         #
#                                                            #
# October 2002, Johannes Schweitzer, NORSAR                  #
#                                                            #
#  johannes@norsar.no                                        #
#                                                            #
#  tested for Linux f77 compiler                             #
#                                                            #
##############################################################


#F77 = /usr/bin/g77 -fno-automatic -W -Wunused -O -fbounds-check
F77 = f77 -fno-automatic 
PREFIX = /usr/local


CC  = /usr/bin/gcc

C-FLAG = 
#C-FLAG =  -g

L-FLAG = -lm

#
# Users using their own tau-spline travel-time tables have to link
# with their own version of the tau-spline software!
#
TTLIB = libtau_h.o 
DEST = ../bin_l
prog = hyposat 

obj1 = hyposat.o  
obj2 = hyposat_lib.o hyposat_geo.o hyposat_name.o hyposat_clib.o hyposat_loc.o
obj3 = hyposat_crust.o hyposat_crust_mod.o hyposat_mag.o
obj4 = hyposat_cross.o hyposat_time.o hyposat_gmi.o hyposat_file.o $(TTLIB)


$(prog) : $(obj1) $(obj2) $(obj3) $(obj4)
	$(F77) $(L-FLAG) -o $(prog) $(obj1) $(obj2) $(obj3) $(obj4) 

	mv hyposat $(DEST)/hyposat

hyposat_cross.o : hyposat_cross.f
	$(F77) $(F-FLAG) -c hyposat_cross.f

hyposat_time.o : hyposat_time.f
	$(F77) $(F-FLAG) -c hyposat_time.f

hyposat_geo.o : hyposat_geo.f
	$(F77) $(F-FLAG) -c hyposat_geo.f

hyposat_gmi.o : hyposat_gmi.f
	$(F77) $(F-FLAG) -c hyposat_gmi.f

hyposat_name.o : hyposat_name.f
	$(F77) $(F-FLAG) -c hyposat_name.f

hyposat_loc.o : hyposat_loc.f
	$(F77) $(F-FLAG) -c hyposat_loc.f

hyposat_crust.o : hyposat_crust.f
	$(F77) $(F-FLAG) -c hyposat_crust.f

hyposat_crust_mod.o : hyposat_crust_mod.f
	$(F77) $(F-FLAG) -c hyposat_crust_mod.f

hyposat_lib.o : hyposat_lib.f
	$(F77) $(F-FLAG) -c hyposat_lib.f

hyposat_file.o : hyposat_file.f
	$(F77) $(F-FLAG) -c hyposat_file.f

hyposat_mag.o : hyposat_mag.f
	$(F77) $(F-FLAG) -c hyposat_mag.f

hyposat.o : hyposat.f
	$(F77) $(F-FLAG) -c hyposat.f

libtau_h.o : libtau_h.f
	$(F77) $(F-FLAG) -c libtau_h.f

hyposat_clib.o : hyposat_clib.c
	$(CC) $(C-FLAG) -c hyposat_clib.c

.PYONY: clean

clean:
	rm -f *.o ../bin_l/hyposat

PREFIX = /usr/local

.PHONY: install
install: hyposat 
	mkdir -p $(DESTDIR)$(PREFIX)/bin
		    cp $(DEST)/$(prog) $(DESTDIR)$(PREFIX)/bin/$(prog)

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/$(prog)
