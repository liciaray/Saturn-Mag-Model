#makefile for the SaturnMagModel

default: saturn_test

ifeq ($(F77),gfortran)
$(info gfortran)
F77 = /opt/homebrew/Cellar/gcc/13.2.0/bin/gfortran
FLIBS = -L/opt/homebrew/lib
F77FLAGS = --fast-math
#FFLAGS=-O1 --unroll-loops --fast-math -ffixed-line-length-0 -fno-automatic -frecord-marker=4 -fconvert=swap 
endif

OBJECTS = find_footprint.o saturn_int.o saturn_ext.o \
	cart_field_to_sphere.o dipoleshield.o igrf_geo_08.o \
	read_pdyn.o tailsheet.o tailsheet_sym.o tailsheet_asym.o \
	bowldeform.o cartharmonic.o cartharmonic_alt.o \
	igrf_geo_cart.o rotate_about_y.o bessjj.o utc_to_tdt.o \
	cart_field_to_cyl.o find_left_index.o ksun.o bowldist.o \
        line_trace.o create_filename.o read_events.o

.SUFFIXES: .for .f

.for.o:
	$(F77) $(F77FLAGS) -c -o $@ $<
.f.o:
	$(F77) $(F77FLAGS) -c -o $@ $< 

saturn_test: $(OBJECTS)
	$(F77) -o $@ $(F77FLAGS) $(OBJECTS)

clean:
	rm -f *.o saturn_test

#gfortran example_saturn.for saturn_int.f saturn_ext.f cart_field_to_sphere.f dipoleshield.f igrf_geo_08.f read_pdyn.f tailsheet.f tailsheet_sym.f tailsheet_asym.f bowldeform.f cartharmonic.f cartharmonic_alt.f igrf_geo_cart.f rotate_about_y.f bessjj.f utc_to_tdt.f cart_field_to_cyl.f find_left_index.f ksun.f bowldist.f -o example_saturn
