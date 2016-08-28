#!/bin/bash
#COBALT -t 20
#COBALT -n 256
#COBALT -A Photovoltaics
#COBALT --disable_preboot

echo $COBALT_PARTNAME

/soft/interpreters/python-2.7.9/powerpc64-linux-gnu/bin/python2.7 ../src/GAtor_master.py ./Cetus.conf

