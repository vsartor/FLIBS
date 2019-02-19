## Makefile
## Copyright (C) 2018 Victhor Sartorio
## * This Source Code Form is part of the FLIBS project.
## * This Source Code Form is subject to the terms of the Mozilla Public
##   License, v. 2.0. If a copy of the MPL was not distributed with this
##   file, You can obtain one at http://mozilla.org/MPL/2.0/.
## * This Source Code Form may use code originally licensed under possibly
##   different open software licenses. For copying, read through subroutine
##   and function description comments to check for those instances.

F2PY = f2py
LIBS = -lopenblas -lpthread -lgfortran
FC   = gfortran
FF   = -march=native -O3 -Wall -Wextra -flto

all: flibs-lib

flibs-lib:
	$(FC) -c $(FF) src/math.f90
	$(FC) -c $(FF) src/random.f90
	$(F2PY) $(LIBS) -c --f90exec=$(FC) --f90flags="$(FF)" \
	    src/random.f90 src/math.f90 \
	    -m flibs
	rm -rf *.mod *.o
