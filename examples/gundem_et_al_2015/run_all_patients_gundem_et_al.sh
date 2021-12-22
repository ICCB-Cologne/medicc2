#!/bin/bash
set -euo pipefail

output_dir="../output_gundem_et_al_2015"

for file in $(ls | grep PTX)
do
    patient=$(echo $file | egrep "PTX0.." -o)
    echo "$patient"
    python ../../medicc2 $file $output_dir --prefix $patient
done
