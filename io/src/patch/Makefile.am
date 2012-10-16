## Process w/ automake.  Or, autoreconf; make ##
AUTOMAKE_OPTIONS = foreign
SUBDIRS = . colors doc examples
INCLUDES = -I. -I$(srcdir) -I$(prefix)/include -I$(prefix)/toolkit/include

includedir = $(prefix)/include 
colordir = $(libdir)/colors

lib_LTLIBRARIES = librsl.la

librsl_la_LDFLAGS = -version-info 1:43
librsl_la_SOURCES = \
$(rapic_c) $(radtec_c)\
dorade.c dorade_print.c dorade_to_radar.c\
lassen.c lassen_to_radar.c \
 radar.c volume.c image_gen.c cappi.c fraction.c read_write.c farea.c \
 range.c radar_to_uf.c uf_to_radar.c wsr88d_to_radar.c \
 carpi.c cube.c sort_rays.c toga_to_radar.c gts.c histogram.c \
 ray_indexes.c anyformat_to_radar.c get_win.c endian.c mcgill_to_radar.c \
 mcgill.c interp.c toga.c wsr88d.c wsr88d_get_site.c wsr88d_m31.c \
 gzip.c prune.c reverse.c fix_headers.c \
 nsig_to_radar.c nsig.c nsig2_to_radar.c \
 africa_to_radar.c africa.c \
 radar_to_hdf_2.c hdf_to_radar.c toolkit_memory_mgt.c \
 radar_to_hdf_1.c rainbow.c rainbow_to_radar.c FourDD.c  \
 findRay.c firstGuess.c firstGuess_noread.c \
 prepVolume.c unfoldVolume.c window.c $(headers)

librsl_la_DEPENDENCIES = $(build_headers)

build_headers = rsl.h wsr88d.h toolkit_1BC-51_appl.h

headers = africa.h dorade.h lassen.h \
          mcgill.h nsig.h radtec.h rainbow.h \
          rapic_routines.h toga.h \
		  $(build_headers)

rapic_c =  rapic_to_radar.c rapic.y rapic-lex.l rapic_routines.c
radtec_c = radtec_to_radar.c radtec.c


rsl.h: Makefile
	@for h in $(build_headers); do \
	echo -n "Checking substitutions in header file $$h ... "; \
	cp $$h $$h.in; \
	sed -e 's/RSL_VERSION_STR.*/RSL_VERSION_STR \"$(VERSION)\"/' \
         -e 's|#define COLORDIR.*|#define COLORDIR \"$(colordir)\"|' \
         -e 's|#define WSR88D_SITE_INFO_FILE.*|#define WSR88D_SITE_INFO_FILE \"$(libdir)/wsr88d_locations.dat\"|' \
    < $$h.in > $$h.new; \
	if cmp -s $$h $$h.new; then \
       rm $$h.new; \
       echo "$$h remains untouched."; \
    else \
       mv $$h.new $$h; \
       echo "substitutions made in $$h."; \
    fi; \
	rm -f $$h.in; \
    done

PREFIX = rapic
LFLAGS = -P$(PREFIX)
YFLAGS = -d -p $(PREFIX)
LEX_OUTPUT_ROOT = lex.$(PREFIX)


install-exec-hook:
	$(INSTALL) -d  $(includedir)
	$(INSTALL) -m 644 rsl.h $(includedir)
	$(INSTALL) -m 644 toolkit_1BC-51_appl.h $(includedir)
	$(INSTALL) -m 644 wsr88d_locations.dat $(libdir)

EXTRA_DIST = CHANGES Copyright GPL LGPL wsr88d_locations.dat rapic.h

DISTCLEANFILES = rapic.c rapic-lex.c
