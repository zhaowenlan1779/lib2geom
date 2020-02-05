#!/bin/bash

mkdir -p build-as-subproject
(cd build-as-subproject; cmake ..  -D2GEOM_AS_SUBPROJECT=ON)
(cd build-as-subproject; make main -j 2)
(cd build-as-subproject; ./main)

mkdir -p build-with-find-package
(cd build-with-find-package; cmake .. -D2GEOM_AS_SUBPROJECT=OFF -D2Geom_DIR="$PWD/../../../opt/lib/cmake/2Geom")
(cd build-with-find-package; make main -j 2)
(cd build-with-find-package; ./main)
