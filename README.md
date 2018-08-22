# AlignmentFree
Fast alignment-free software

## Prerequisites:
### Program:
* [Python3](https://www.python.org/downloads/release/python-363/) or [Anaconda3](https://conda.io/docs/user-guide/install/download.html)
### Packages:
* Required Python3 packages: numpy, scipy, sklearn-learn.
* We recommend use [Anaconda3](https://conda.io/docs/user-guide/install/download.html) to install all required packages:
```
conda install numpy scipy scikit-learn
```

## Installing:
Clone this repository into your local directory:
```
git clone https://github.com/GeniusTang/AlignmentFree.git
```

Unix setup:
```
CC=g++ python setup.py install --install-platlib=./src/
```
Mac setup:
```
MACOSX_DEPLOYMENT_TARGET=10.9 CC=g++ python setup.py install --install-platlib=./src/
```

## Runing the tests 
Example1:
```
python alignmentfree.py -r -a d2star,d2shepp,CVtree,Ma,Eu,d2 -k 5 -m 0 -f test_file.txt -t 12 -d test_count/ -o test_result/
```
Example2:
```
python alignmentfree.py -r -a d2star,d2shepp,CVtree,Ma,Eu,d2 -k 5 -m 0 -f1 test_file_1.txt -f2 test_file_2.txt -t 12 -d test_count/ -o test_result/
```
