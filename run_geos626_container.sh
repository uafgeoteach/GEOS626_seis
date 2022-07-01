#!/bin/bash

set -ex

docker run -it --rm -p 8888:8888 -v $(pwd)/home:/home/jovyan geos626_seis:latest