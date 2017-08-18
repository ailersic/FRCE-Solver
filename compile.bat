del /Q .\mod\*
del /Q .\obj\*
del /Q .\bin\*

gfortran .\src\pentadiag.f95 -c -J .\mod -o .\obj\pentadiag.o
gfortran .\src\xyEquations.f95 -c -J .\mod -o .\obj\xyEquations.o
gfortran .\src\zEquation.f95 -c -J .\mod -o .\obj\zEquation.o
gfortran .\src\arcLength.f95 -c -J .\mod -o .\obj\arcLength.o
gfortran .\src\ropeCoil.f95 -c -J .\mod -o .\obj\ropeCoil.o

gfortran .\obj\*.o -o .\bin\ropeCoil.exe