# SPLATT + Chapel

Port of SPLATT to Chapel

Currently how to compile:

chpl -I${CBLAS_DIR} -L${BLAS_LIBS} -lblas --permit-unhandled-module-errors mats.chpl args.chpl base.chpl CSF.chpl splatt_IO.chpl splatt_sort.chpl sptensor.chpl stats.chpl util.chpl kruskal.chpl CPD.chpl main.chpl
