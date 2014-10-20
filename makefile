all: gen_source.c
	mpicc -g gen_source.c -o gensource
