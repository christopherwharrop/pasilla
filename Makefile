all:
	make -C utilities
	make -C models
	make -C varsolver
	make -C varsolveroo

models: utilities
	make -C models

varsolver: utilities
	make -C varsolver

varsolveroo: utilities
	make -C varsolveroo

utilties:
	make -C utilities

clean:
	make clean -C utilities
	make clean -C models
	make clean -C varsolver
	make clean -C varsolveroo
