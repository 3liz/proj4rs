#!/bin/bash

(
  set -e; 
  cd .docker; 
  docker build -t ol-proj4rs-test .;
)

docker run \
    --entrypoint=bash \
    --name ol-proj4rs-test \
    --network host \
    --rm \
    -it \
    -w /src \
    -v $(pwd):/src -u $UID:$UID \
    ol-proj4rs-test .docker/ol-proj4rs-test.sh

