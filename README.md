# AlignmentFree
Fast alignment-free software


Installation:

git clone https://github.com/GeniusTang/AlignmentFree.git

For unix:

CC=g++ python setup.py install --install-platlib=./src/

For Mac:

MACOSX_DEPLOYMENT_TARGET=10.9 CC=g++ python setup.py install --install-platlib=./src/

Example1:

python alignmentfree.py -r -a d2star,d2shepp,CVtree,Ma,Eu,d2 -k 5 -m 0 -f test_file.txt -t 12 -d test_count/ -o test_result/

Example2:

python alignmentfree.py -r -a d2star,d2shepp,CVtree,Ma,Eu,d2 -k 5 -m 0 -f1 test_file_1.txt -f2 test_file_2.txt -t 12 -d test_count/ -o test_result/
