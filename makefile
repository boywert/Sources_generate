all: gen_source gen_source_nphot

gen_source: gen_source.c
	gcc -g gen_source.c -o gen_source -lm
gen_source_nphot: gen_source_nphot.c
	gcc -g gen_source_nphot.c -o gen_source_nphot -lm
clean:
	rm -f gen_source gen_source_nphot
