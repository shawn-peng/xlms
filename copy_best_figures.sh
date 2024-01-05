#!/bin/bash

cd $1 ;
i=1 ;
for d in $(find -mindepth 1 -maxdepth 1 -type d);
do
   	cp $d/best.png $i.png ;
   	((i++)) ;
done

