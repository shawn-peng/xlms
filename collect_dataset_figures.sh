#!/bin/bash

mkdir -p $2_figures
for f in figures*/$1.png;
do
	cp $f $2_figures/$(echo $f | sed -e "s/figures_\(.*\)\/$1\.png/\1.png/") ;
done
