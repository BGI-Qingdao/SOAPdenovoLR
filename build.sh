#!/usr/bin/env bash

tag=`date +"%Y-%m-%d_%H:%M:%S"`
rename SOAPdenovo "SOAPdenovo_"$tag"_" *
make clean all
make debug=1 clean all
make 127mer=1 clean all
make 127mer=1 debug=1 clean all
