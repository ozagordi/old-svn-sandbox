#!/bin/bash

nr=100

## A simple example (ignoring error correction)

## in.popl is the initial population (in amino acids)
## 1. sample 1000 pyrosequencing reads of length 30 with no error

../pyroseq.pl  -e 0 -k 30 -n $nr -A sample
mv sample.$nr.aln.cor.fas sample.$nr.fas

## output: sample.1000.fas


## 2. translate to read format
## output: sample.$nr.read
../fas2read.pl -f sample.$nr.fas

## 3. eliminate redundant reads
## output: sample.$nr.rest
../contain  -f sample.$nr

## 4. run maximum matching, output up to 200 haplotypes
## output: sample.$nr.geno
../mm.py sample.$nr.rest 200
#python2.5 -m cProfile ../mm.py sample.$nr.rest 200

## 5. run EM
## output: sample.$nr.popl
../freqEst -f sample.$nr

