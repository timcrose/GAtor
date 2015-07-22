#!/bin/bash
#COBALT -t 20
#COBALT -n 256
#COBALT -A GAtor
#COBALT --disable_preboot

echo $COBALT_PARTNAME
python /path/to/gator/src/core/master.py -i -f sample.conf  

