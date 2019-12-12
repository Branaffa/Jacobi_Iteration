# Jacobi_Iteration

For running and compiling the serial code run:
gcc -fopenmp jacobi_serial.c -o jacobi_serial //
./jacobi_serial


For compiling and running the openmp code, run the following:
gcc -fopenmp jacobi_parralel.c -o jacobi_parallel //
./jacobi_parallel

For running the openACC code, run the following:
gcc -fopenacc jacobi_acc.c -o jacobi_acc //
./jacobi_acc

For when the prompt asks for a file to run for the matrix use any of the following:
m32.txt
m64.txt
m98.txt
m128.txt

