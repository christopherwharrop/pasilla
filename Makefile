all:
	make -C utilities
	make -C models
	make -C varsolver

clean:
	make clean -C utilities
	make clean -C models
	make clean -C varsolver
