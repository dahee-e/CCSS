# Code for CCSS Algorithm

This repository contains a reference implementation of the algorithms for the paper:
“CCSS: Towards Conductance-based Community Search with Size Constraints”.

## Experimental Environment
This code can be run on a Linux environment. Our code is compiled on a Linux server with Intel(R) Xeon(R) Gold 6348 @ 2.60GHz CPU and 256GB RAM running CentOS 7.7.
## Dataset Description
We guarantee that the maximum connected component of the original data is taken. The IDs of the vertex are arranged consecutively starting from 0, with the format:
```shell
 first_v  second_v
 ```
All public datasets can be downloaded from http://snap.stanford.edu/data/index.html

## Running Format
python3 -u CCSS.py [1. graph name] [2. iterations] [3. query vid: q][4. size LB: l] [5. size UB: h] 

## Running Instance
python3 -u CCSS.py graph.txt 2 3 3 4
