EXE = polymer.exe
FILES = variables.f90 util_random.f90 util.f90 sumup.f90 moves.f90 main.f90
LOC = .
FFLAG = O3
#  ---------------------
intel: $(FILES)
	ifort $(FILES) -o $(LOC)/$(EXE)
#	ifort -c variables.f90
#	ifort -c util.f90
#	ifort -c sumup.f90
#	ifort -c moves.f90
#	ifort -c main.f90 
#	ifort *.o -o $(LOC)/$(EXE)

#  ---------------------
gnu:	$(FILES)
##	/opt/local/bin/gfortran-mp-4.7  $(FILES) -o $(LOC)/$(EXE)
	/usr/bin/gfortran -g -$(FFLAG) $(FILES) -o $(LOC)/$(EXE)
#	gfortran -c -g -$(FFLAG) sumup.f90
#	gfortran -c -g -$(FFLAG) moves.f90
#	gfortran -c -g -$(FFLAG) main.f90
#	gfortran *.o -o -g -$(FFLAG) $(LOC)/$(EXE)

#  ---------------------
clean:
	rm -f $(LOC)/$(EXE) *.mod *.o $(LOC)/*.out movie.xyz config.dat results.dat dum.dat
