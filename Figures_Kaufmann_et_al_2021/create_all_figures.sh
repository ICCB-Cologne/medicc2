#!/bin/bash

for f in $(ls | egrep "^(Fig|Supp).*py"); do 
    (python $f >> /dev/null; echo "Finished file $f") &
done
wait
