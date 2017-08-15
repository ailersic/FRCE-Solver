gfortran .\src\xyEquations.f95 -c -J .\mod -o .\obj\xyEquations.o
gfortran .\src\zEquation.f95 -c -J .\mod -o .\obj\zEquation.o
gfortran .\src\zEquationSOR.f95 -c -J .\mod -o .\obj\zEquationSOR.o
gfortran .\src\arcLength.f95 -c -J .\mod -o .\obj\arcLength.o
gfortran .\src\ropeCoil.f95 -c -J .\mod -o .\obj\ropeCoil.o
gfortran .\obj\*.o -o .\bin\ropeCoil.exe