#!/bin/bash
mkdir -p deps-prefix
cmake -B./deps-prefix -H../deps
cd deps-prefix
make
cd -

