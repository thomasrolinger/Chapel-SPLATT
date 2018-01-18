INC_DIRS=${CBLAS_DIR}
LIB_DIRS=${BLAS_LIBS}
LIBS=-lopenblas -lgfortran -lpthread
CHPL_FLAGS=--permit-unhandled-module-errors -O --specialize --fast --no-llvm --local --mem=${CHPL_MEM} --tasks=${CHPL_TASKS} --ccflags="-std=c99 -funroll-loops -fgnu89-inline -fstrict-aliasing -fPIC"

PROFILE=0
GEN_C_CODE=0
DEBUG=0
ifeq ($(PROFILE), 1)
	CHPL_FLAGS += --ccflags="-pg -g" --ldflags="-pg"
endif

ifeq ($(GEN_C_CODE), 1)
	CHPL_FLAGS += --savec=${PWD}/generated_C/src
endif

ifeq ($(DEBUG), 1)
    CHPL_FLAGS += --ccflags="-g" --ldflags="-g"
endif

CHPL=chpl
EXEC=./bin/splatt
SOURCES=./src/*.chpl

all:
	if [ $(GEN_C_CODE) -eq 1 ]; then \
        mkdir -p generated_C; \
    fi
	$(MAKE) ${EXEC}
	if [ $(GEN_C_CODE) -eq 1 ]; then \
		mkdir -p generated_C/bin; \
		mv generated_C/src/Makefile generated_C; \
	fi

${EXEC}: ${SOURCES}
	${CHPL} -o ${EXEC} ${CHPL_FLAGS} -I${INC_DIRS} -L${LIB_DIRS} ${LIBS} ${SOURCES}

clean:
	rm -rf ./bin/splatt ./generated_C/bin/splatt
