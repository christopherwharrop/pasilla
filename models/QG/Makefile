all:
	make -C src/tapenade/ADFirstAidKit
	make -C src
	cp src/*.exe exe
	cp src/*.mod include
	cp src/*.a lib

clean:
	make clean -C src/tapenade/ADFirstAidKit
	make clean -C src
	rm -f exe/*.exe
	rm -f include/*.mod
	rm -f lib/*.a
