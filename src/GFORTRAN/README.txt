Sep 2024 FEHM V3.6

These files have platform relevant PC lines commented out for for UNIX
and are used by src/Makefile to compile xfehm for unix machines.

These files in GFORTRAN are used instead of the files in src.
inrestart.f  insptr.f  mainrip.f

These files are included here to handle gfortran builds where binary is not supported.
Makefile will use these files instead of the files in src

Makefile:

inrestart.o : ${ALTDIR}inrestart.f
	${FC} ${FFLAGS} $< -c -o ${SRCDIR}$@

insptr.o : ${ALTDIR}insptr.f
	${FC} ${FFLAGS} $< -c -o ${SRCDIR}$@


Note: change version and date in dated.template
Used in Makefile to set version and date
.PHONY: dated
dated:
	rm -f ${SRCDIR}dated.f
	rm -f ${SRCDIR}dated.o
	sed s/OS\ DATE/${OSTAG}\ ${DATETAG}/ ${SRCDIR}dated.template > ${SRCDIR}dated.f

