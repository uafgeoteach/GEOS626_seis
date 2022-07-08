#!/bin/bash

set -ex

docker run -it --rm -p 8888:8888 -v $(pwd):/home/jovyan/GEOS626_seis geos626_seis:latest

# if [ ! -d "home" ] 
# then
#     mkdir -p home/GEOS626_seis
# fi

# docker run -it --rm -p 8888:8888 -v $(pwd)/home:/home/jovyan geos626_seis:latest