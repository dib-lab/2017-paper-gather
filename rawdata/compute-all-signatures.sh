#! /bin/bash
for i in ?.fa ??.fa
do
    sourmash compute --scaled 10000 -k 21,31,41,51 $i --name-from-first
done
