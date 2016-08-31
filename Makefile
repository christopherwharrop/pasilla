all:
	make -C utilities
	make -C models
	make -C varsolver
	make -C varsolveroo

clean:
	make clean -C utilities
	make clean -C models
	make clean -C varsolver
	make clean -C varsolveroo
