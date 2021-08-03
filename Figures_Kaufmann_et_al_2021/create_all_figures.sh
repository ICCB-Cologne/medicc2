#!/bin/bash

for f in $(ls | grep -e "Fig.*py"); do 
    (python $f >> /dev/null; echo "Finished file $f") &
done
wait
