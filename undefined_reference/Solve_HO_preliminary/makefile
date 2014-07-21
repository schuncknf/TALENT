# SHELL=/bin/sh

a:   
	c++ MY_HO.cpp -o HO.x -larmadillo -lgsl  

b:
	c++ MY_HO.cpp -o HO.x -larmadillo -lgsl -DTeasy
	
help:
	@echo 
	@echo "-----------------------------------------------------"
	@echo "-  type: <<make a>> to solve the H1 with the brute diagonalization of the Kinetic energy"
	@echo "-  type: <<make b>> to solve the H1 with th smart, pre-integrated Kinetic energy"
	@echo "-  to change the variational coefficients and the size of the basis go inside 'MY_HO.cpp' and change the defines... I know you can do that."
	@echo "-----------------------------------------------------"
	@echo ""
