#
CFLAGS=	-O
INSTDIR=/usr/local/lib

OBJS=	arc.o box.o circle.o close.o dot.o erase.o label.o \
	line.o linemod.o move.o open.o point.o space.o subr.o

lib4014.a:	${OBJS}
	ar rc lib4014.a ${OBJS}
	-ranlib lib4014.a

install:	lib4014.a
	cp lib4014.a $(INSTDIR)

clean:
	rm -f $(OBJS) lib4014.a
