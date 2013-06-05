all: dos_kpt.f90
	ifort dos_kpt.f90 -o exe.x -O3 
#-stand f95
dft: dos_kpt_fixed.f90
	ifort dos_kpt_fixed.f90 -o fixed.x -O3

clean:
	rm *.x
