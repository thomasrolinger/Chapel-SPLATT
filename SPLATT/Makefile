INC_DIRS=${CBLAS_DIR}
LIB_DIRS=${BLAS_LIBS}
LIBS=-lblas
MISC_FLAGS=--permit-unhandled-module-errors

CHPL=chpl
EXEC=./bin/splatt

SOURCES=./src/*.chpl

${EXEC}: ${SOURCES}
	${CHPL} -o ${EXEC} ${MISC_FLAGS} -I${INC_DIRS} -L${LIB_DIRS} ${LIBS} ${SOURCES}

clean:
	rm -rf ./bin/splatt
