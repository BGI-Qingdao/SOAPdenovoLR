#!/usr/bin/env bash

make clean all
make debug=1 clean all
make 127mer=1 clean all
make 127mer=1 debug=1 clean all
