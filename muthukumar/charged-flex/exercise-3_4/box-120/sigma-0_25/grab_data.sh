# !/bin/bash

rsync -R samhoover@muthugroup:$(pwd)/./RUN_*/{radius*,atom*,dump*} .
