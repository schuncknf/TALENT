
all: mtrx_test sho_print sho_test print_gaulag sho_test_gauss hydrogen print_gauher sho_test_gauher

mtrx_test: mtrx_test.c
	gcc -o mtrx_test mtrx_test.c -llapacke

sho_print: sho_print.c sho.c
	gcc -o sho_print sho_print.c sho.c -lm

sho_test: sho_test.c sho.c
	gcc -o sho_test sho_test.c sho.c -lm

print_gaulag: print_gaulag.c gaulag.c
	gcc -o print_gaulag print_gaulag.c gaulag.c -lm

sho_test_gauss: sho_test_gauss.c sho.c gaulag.c
	gcc -o sho_test_gauss sho_test_gauss.c sho.c gaulag.c -lm

hydrogen: hydrogen.c sho.c eigen.c gaulag.c
	gcc -o hydrogen hydrogen.c sho.c eigen.c gaulag.c -lm -llapacke

print_gauher: print_gauher.c gauher.c
	gcc -o print_gauher print_gauher.c gauher.c -lm

sho_test_gauher: sho_test_gauher.c sho.c gauher.c
	gcc -o sho_test_gauher sho_test_gauher.c sho.c gauher.c -lm
