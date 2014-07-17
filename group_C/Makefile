
all: mtrx_test sho_print sho_test gaulag sho_test_gauss hydrogen

mtrx_test: mtrx_test.c
	gcc -o mtrx_test mtrx_test.c -llapacke

sho_print: sho_print.c sho.c
	gcc -o sho_print sho_print.c sho.c -lm

sho_test: sho_test.c sho.c
	gcc -o sho_test sho_test.c sho.c -lm

gaulag: gaulag.c
	gcc -o gaulag gaulag.c -lm

sho_test_gauss: sho_test_gauss.c sho.c
	gcc -o sho_test_gauss sho_test_gauss.c sho.c -lm

hydrogen: hydrogen.c sho.c
	gcc -o hydrogen hydrogen.c sho.c -lm -llapacke
