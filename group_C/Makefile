
all: mtrx_test sho_print sho_test print_gaulag sho_test_gauss hydrogen print_gauher sho_test_gauher hydrogen-bench

mtrx_test: mtrx_test.c
	gcc -o mtrx_test mtrx_test.c -llapacke

sho_print: sho_print.c sho.c
	gcc -o sho_print sho_print.c sho.c -lm -lgsl -lgslcblas

sho_test: sho_test.c sho.c
	gcc -o sho_test sho_test.c sho.c -lm -lgsl -lgslcblas

print_gaulag: print_gaulag.c gaulag.c
	gcc -o print_gaulag print_gaulag.c gaulag.c -lm

sho_test_gauss: sho_test_gauss.c sho.c gaulag.c
	gcc -o sho_test_gauss sho_test_gauss.c sho.c gaulag.c -lm -lgsl -lgslcblas

hydrogen: hydrogen.c sho.c eigen.c gaulag.c
	gcc -o hydrogen hydrogen.c sho.c eigen.c gaulag.c -lm -llapacke -lgsl -lgslcblas

print_gauher: print_gauher.c gauher.c
	gcc -o print_gauher print_gauher.c gauher.c -lm

sho_test_gauher: sho_test_gauher.c sho.c gauher.c
	gcc -o sho_test_gauher sho_test_gauher.c sho.c gauher.c -lm -lgsl -lgslcblas

hydrogen-bench: hydrogen-bench.c sho.c eigen.c gaulag.c
	gcc -o hydrogen-bench hydrogen-bench.c sho.c eigen.c gaulag.c -lm -llapacke -lgsl -lgslcblas
