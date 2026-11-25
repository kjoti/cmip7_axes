#
# cmip7_axes -- Generate GTAXLOC files for CMIP7 experiments
#
# Required libraries:
#   * libgtool3
#
SHELL	= /bin/sh
RM	= rm -f

#PREFIX	= /usr/local
PREFIX	= $(HOME)/.local
GTAXDIR = $(PREFIX)/share/gtool3/axes

CC	= gcc
CFLAGS	= -Wall -pedantic -O2 -I $(PREFIX)/include
LDFLAGS = -L $(PREFIX)/lib -Wl,-rpath=$(PREFIX)/lib

PROGRAM = cmip7_axes.exe

OBJS	= \
	hybrid.o \
	main.o \
	ocean.o \
	plevs.o \
	utils.o

LIBS	= -lgtool3 -lm

$(PROGRAM): $(OBJS)
	$(CC) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)

all: $(PROGRAM)

install: $(PROGRAM)
	@./$(PROGRAM)
	@echo "Copy axis files into $(GTAXDIR)"
	@for x in GTAXLOC.*; do \
	    cp $${x} $(GTAXDIR); \
            echo "Copy $${x} ..."; \
	done

clean:
	@$(RM) $(OBJS) $(PROGRAM)

distclean: clean
	@$(RM) -f GTAXLOC.*
