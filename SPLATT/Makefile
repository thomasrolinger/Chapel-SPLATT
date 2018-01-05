INC_DIRS=${CBLAS_DIR}
LIB_DIRS=${BLAS_LIBS}
LIBS=-lopenblas -lpthread -lgfortran
MISC_FLAGS=--permit-unhandled-module-errors -O --specialize

# This will turn off all runtime checks. It makes the code a lot faster
# but defeats the purpose of Chapel, in a sense. But it is useful to use
# to provide a "goal" for the code to reach
FAST_FLAGS=--no-checks



CHPL=chpl
EXEC=./bin/splatt



SOURCES=./src/*.chpl

${EXEC}: ${SOURCES}
	${CHPL} -o ${EXEC} ${MISC_FLAGS} -I${INC_DIRS} -L${LIB_DIRS} ${LIBS} ${SOURCES}

clean:
	rm -rf ./bin/splatt
