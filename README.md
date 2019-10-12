# Afann
Afann (Alignment-Free methods Adjusted by Neural Network) is an alignment-free software that supports fast calculation of different dissimilarity measures including d2star, d2shepp, CVtree, Manhattan, Euclidean and d2. It also adjusts the bias of d2star and d2shepp calculated from NGS samples.

## Prerequisites:
### Program:
* [Python3](https://www.python.org/downloads/release/python-363/) or [Anaconda3](https://www.anaconda.com/download/)
### Packages:
* Required Python3 packages: numpy, sklearn-learn.
* We recommend use [Anaconda3](https://www.anaconda.com/download/) to install all required packages:
```
conda install numpy scipy scikit-learn
```

## Installing:
Clone this repository into your local directory:
```
git clone https://github.com/GeniusTang/Afann.git
```

Unix setup:
```
CC=g++ python setup.py install --install-platlib=./src/
```
Mac setup:
```
MACOSX_DEPLOYMENT_TARGET=10.9 CC=g++ python setup.py install --install-platlib=./src/
```

## Runing the tests:
### Example1: 
Calculate pairwise d2star,d2shepp,CVtree,Ma,Eu,d2 distances among all samples listed in test_file.txt, using kmer length 5, Markovian order 0. 
* -r: Consider reverse complement of kmers. 
* -t: Use 8 threads.
* -d: Save kmer counts in test_count/
* -o: Save outputs in test_result/
```
python afann.py -r -a d2star,d2shepp,CVtree,Ma,Eu,d2 -k 5 -m 0 -f test_file.txt -t 8 -d test_count/ -o test_result/
```
### Example2: 
Calculate pairwise d2star,d2shepp distances among all samples listed in test_file.txt, using kmer length 5, Markovian order 0 with bias adjustment.
* -r: Consider reverse complement of kmers.
* -t: Use 8 threads.
* -d: Save kmer counts in test_count/
* -o: Save outputs in test_result/
* --adjust: Bias adjustment for NGS samples.
```
python afann.py -r -a d2star,d2shepp -k 5 -m 0 -f test_file.txt -t 8 -d test_count/ -o test_result/ --adjust
```
### Example3:
Calculate pairwise d2star,d2shepp,CVtree,Ma,Eu,d2 distances among all samples listed in test_file_1.txt and all samples listed in test_file_2.txt, using kmer length 5, Markovian order 0.
```
python afann.py --slow -r -a d2star,d2shepp,CVtree,Ma,Eu,d2 -k 5 -m 0 -f1 test_file_1.txt -f2 test_file_2.txt -t 8 -d test_count/ -o test_result/test
```
* -o: Save outputs in test_result/ with prefix test
* --slow: Calculate with less memory usage, but slower
### Example4:
Calculate pairwise d2star,d2shepp distances among all sequences listed in test_samples/crm.fa, using kmer length 5, Markovian order 1.
```
python afann.py -r -a d2star,d2shepp -k 5 -m 1 -s test_samples/crm.fa -t 8 -d test_count/ -o test_result/crm
```
* -o: Save outputs in test_result/ with prefix crm
### Example5:
Calculate the Markovian orders of all sequences listed in test_file.txt.
```
python afann.py -r --BIC -k 5 -f test_file.txt -t 8 -d test_count/ -o test_result/test
```
## Usage:
```
usage: afann.py [-h] [-a METHOD] -k K [-m M] [-f FILENAME]
                        [-s SEQUENCE_FILE] [-f1 FILENAME1] [-f2 FILENAME2]
                        [-s1 SEQUENCE_FILE_1] [-s2 SEQUENCE_FILE_2] [-d DIR]
                        [-o OUTPUT] [-t THREADS] [-r] [--adjust] [--BIC]
                        [--slow]
```

Optional arguments:
```
  -h, --help           show this help message and exit
  -a METHOD            A list of alignment-free method, separated by comma:
                       d2star,d2shepp,CVtree,Ma,Eu,d2
  -k K                 Kmer length
  -m M                 Markovian Order, required for d2star, d2shepp and
                       CVtree
  -f FILENAME          A file that lists the paths of all samples, cannot be
                       used together with -f1, -f2, -s, -s1, -s2
  -s SEQUENCE_FILE     A fasta file that lists the sequences of all samples,
                       cannot be used together with -f, -f1, -f2, -s1, -s2
  -f1 FILENAME1        A file that lists the paths of the first group of
                       samples, must be used together with -f2, cannot be used
                       together with -f, -s, -s1, -s2
  -f2 FILENAME2        A file that lists the paths of the second group of
                       samples, must be used together with -f1, cannot be used
                       together with -f, -s, -s1, -s2
  -s1 SEQUENCE_FILE_1  A fasta file that lists the sequences of the first
                       group of samples, must be used together with -s2,
                       cannot be used together with -f, -f1, -f2, -s
  -s2 SEQUENCE_FILE_2  A fasta file that lists the sequences of the second
                       group of samples, must be used together with -s1,
                       cannot be used together with -f, -f1, -f2, -s
  -d DIR               A directory that saves kmer count
  -o OUTPUT            Prefix of output (defualt: Current directory)
  -t THREADS           Number of threads
  -r                   Count the reverse complement of kmers (default: False)
  --adjust             Adjust d2star and/or d2shepp distances for NGS samples,
                       -r will be set automatically
  --BIC                Use BIC to estimate the Markovian orders of sequences
  --slow               Use slow mode for calculation with less memory usage
                       (default: False)
```

## Copyright and License Information:
Copyright (C) 2019 University of Southern California

Authors: Kujin Tang, Jie Ren, Fengzhu Sun

This program is freely available at <https://github.com/GeniusTang/Afann> under the terms of USC-RL v1.0.

Commercial users should contact Dr. Sun at <fsun@usc.edu>, copyright at the University of Southern California.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
