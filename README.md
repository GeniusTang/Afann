# AlignmentFree
Fast alignment-free software


Installation:
For unix:
CC=g++ python setup.py install --install-platlib=./src/

For Mac:
MACOSX_DEPLOYMENT_TARGET=10.9 CC=g++ python setup.py install --install-platlib=./src/

Example:
python alignmentfree.py -r -a d2star,d2shepp,CVtree,Ma,Eu,d2 -k 5 -m 0 -f test_file.txt -t 12 -d test_count/ -o test_result/
