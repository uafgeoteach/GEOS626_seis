#!/bin/bash

set -ex

docker run -it --rm -p 8888:8888 -v $(pwd)/..:/home/jovyan/GEOS626_seis geos626_seis:latest