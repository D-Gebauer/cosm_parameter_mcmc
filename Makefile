mcmc: mcmcMethods.c ./data/zUnion.txt ./data/mbUnion.txt ./data/emBUnion.txt
	gcc -Wall main.c mcmcMethods.c -o mcmc -lgsl -lgslcblas -lm -ffast-math

run: mcmcMethods.c ./data/zUnion.txt ./data/mbUnion.txt ./data/emBUnion.txt
	gcc -Wall main.c mcmcMethods.c -o mcmc -lgsl -lgslcblas -lm -ffast-math
	./mcmc

assembler: mcmcMethods.c ./data/zUnion.txt ./data/mbUnion.txt ./data/emBUnion.txt
	gcc -Wall main.c mcmcMethods.c -S -lgsl -lgslcblas -lm -ffast-math

clean:
	rm ./mcmc
