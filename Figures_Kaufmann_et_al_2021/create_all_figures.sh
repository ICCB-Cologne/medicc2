#!/bin/bash

for f in $(ls | grep -e "Fig.*py"); do 
    python $f & 
done
wait
