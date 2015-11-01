all: gen_source gen_source_nphot gen_source_halo

gen_source: gen_source.c
	gcc -g gen_source.c -o gen_source -lm
gen_source_nphot: gen_source_nphot.c
	gcc -g gen_source_nphot.c -o gen_source_nphot -lm
gen_source_halo: gen_source_halo.c
	gcc -g gen_source_halo.c -o gen_source_halo -lm
clean:
	rm -f gen_source gen_source_nphot gen_source_halo
