#!/bin/bash

python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --wgd --prefix wgd
python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --prefix no_wgd
python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --wgd --max-num-wgds 1 --prefix wgd_1
python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --wgd --max-num-wgds 2 --prefix wgd_2
