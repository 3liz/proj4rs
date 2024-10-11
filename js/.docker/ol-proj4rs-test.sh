#!/bin/bash

export npm_config_cache=$(pwd)/.npm

cd ol-proj4rs-demo-app
npm --loglevel=verbose update
echo "Starting ol-proj4rs-demo-app"
npm --loglevel=verbose start -c vite.docker.js


