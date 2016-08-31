all:
	make -C src
	cp src/*.a lib
	cp src/*.mod include

clean:
	make clean -C src
	rm -f include/*.mod
	rm -f lib/*.a
