all:
	make -C src
	cp src/*.exe exe
	cp src/*.mod include

clean:
	make clean -C src
	rm -f exe/*.exe
	rm -f include/*.mod
