INC_DIRS=${CBLAS_DIR}
LIB_DIRS=${BLAS_LIBS}
LIBS=-lopenblas -lpthread -lgfortran
CHPL_FLAGS=--permit-unhandled-module-errors -O --specialize --fast --no-llvm --local --mem=cstdlib

PROFILE=0
GEN_C_CODE=0

ifeq ($(PROFILE), 1)
	CHPL_FLAGS += --ccflags="-pg" --ldflags="-pg"
endif

ifeq ($(GEN_C_CODE), 1)
	CHPL_FLAGS += --savec=${PWD}/generated_C/src
endif

CHPL=chpl
EXEC=./bin/splatt
SOURCES=./src/*.chpl

all:
	$(MAKE) ${EXEC}
	if [ $(GEN_C_CODE) -eq 1 ]; then \
		mkdir -p generated_C/bin; \
		mv generated_C/src/Makefile generated_C; \
	fi

${EXEC}: ${SOURCES}
	mkdir -p generated_C
	${CHPL} -o ${EXEC} ${CHPL_FLAGS} -I${INC_DIRS} -L${LIB_DIRS} ${LIBS} ${SOURCES}

clean:
	rm -rf ./bin/splatt ./generated_C/bin/splatt
