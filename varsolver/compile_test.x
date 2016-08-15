module load intel impi netcdf

mpif90 -f90=ifort -openmp -I/apps/netcdf/4.3.0-intel/include -I/contrib/gptl/gptl-v5.4.4_impi_noomp/include -I/contrib/gptl/gptl-v5.4.4_impi_noomp/include -L/contrib/gptl/gptl-v5.4.4_impi_noomp/lib -lgptl -openmp -O2 -m64 -heap-arrays -vec-report0 -assume buffered_io -c -g -finstrument-functions varsolver_12aug.f90 -openmp

mpif90 -f90=ifort varsolver_12aug.o -I/apps/netcdf/4.3.0-intel/include -I/contrib/gptl/gptl-v5.4.4_impi_noomp/include -L/apps/netcdf/4.3.0-intel/lib -lnetcdff -lnetcdf -mkl=sequential -L/contrib/gptl/gptl-v5.4.4_impi_noomp/lib -lgptl -openmp -o a.out

mv a.out test.x

./test.x > test.out

./hex2name.pl ./test.x ./timing.0 > ./timing.test.txt
./parsetiming.rb -t 1 -d 6 timing.test.txt > timing.test.parse.txt
