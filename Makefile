PLATFORM = $(shell uname)

CFLAGS= -g -Wall -std=c99
LDFLAGS=

ifeq ($(PLATFORM),Darwin)
## Mac OS X
CFLAGS += -m64 -isystem/usr/local/include
LDFLAGS+= -m64 -lc -framework AGL -framework OpenGL -framework GLUT \
	-framework Foundation
else
## Linux
INCLUDEPATH  = -I/usr/include/GL/
#LIBPATH = -L/usr/lib #-L/usr/X11R6/lib
LDFLAGS+=  -lGL -lglut -lm #-lrt -lGLU -lX11 -lm  -lXmu -lXext -lXi
endif

CC = gcc
GRAPHICS = $(LIBPATH) $(LDFLAGS)

default: compile

compile: pqueue.o grid.o simplification.o display.o list.o
	$(CC) $? $(GRAPHICS) -o simplify

simp.o: simplification.c
	$(CC) $(CFLAGS) -c $< -o $@

pQueue.o: pqueue.c
	$(CC) $(CFLAGS) -c $< -o $@

grid.o:	grid.c
	$(CC) $(CFLAGS) -c $< -o $@

display.o: display.c
	$(CC) $(CFLAGS) -c $< -o $@

list.o: llist.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) *.o *~ simplify
