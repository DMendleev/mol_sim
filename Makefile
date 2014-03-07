EXE = polymer.exe
FILES = variables.f90 util_random.f90 util.f90 sumup.f90 moves.f90 main.f90
LOC = .

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
	/usr/bin/gfortran -ggdb $(FILES) -o $(LOC)/$(EXE)
#	gfortran -c -ggdb sumup.f90
#	gfortran -c -ggdb moves.f90
#	gfortran -c -ggdb main.f90
#	gfortran *.o -o -ggdb $(LOC)/$(EXE)

#  ---------------------
clean:
	rm -f $(LOC)/$(EXE) *.mod *.o $(LOC)/*.out config.dat results.dat
