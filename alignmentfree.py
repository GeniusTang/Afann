from src._count import kmer_count
from src._count import kmer_count_m_k
from scipy.spatial.distance import cosine
from sklearn.metrics.pairwise import cosine_similarity
from scipy import stats
import numpy as np
import time
import os
import method
import argparse

Suffix = ['fna', 'fa', 'fasta']
Alphabeta = ['A', 'C', 'G', 'T']
Alpha_dict = dict(zip(Alphabeta, range(4)))

def num2nuc(num, k):
    num_bin = format(num, '0%sb'%(k*2))
    return ''.join([Alphabeta[int(num_bin[i:i+2], 2)] for i in range(0, len(num_bin), 2)])

def shift(num, K, M, x):
    mask = 2**(2*M)-1
    num = num >> ((K-M)*2 - 2*x)
    num &= mask
    return num

def get_sequence_from_file(filename):
    sequence_list = []
    with open(filename) as f:
        for line in f.readlines():
            line = line.strip()
            for suffix in Suffix:
                if line.endswith(suffix):
                    sequence_list.append(line)
    return sequence_list 
   
def get_matrix(a_method):
    if a_method not in ['d2star', 'd2shepp', 'cvtree', 'ma', 'eu', 'd2']:
        print('Invalid method %s'%a_method)
        print('Only d2star,d2shepp,CVtree,Ma,Eu,d2 are supported') 
        raise NameError
    else:
        if a_method == 'd2star':
            return method.d2star_matrix_pairwise
        elif a_method == 'd2shepp':
            return method.d2shepp_matrix_pairwise
        elif a_method == 'cvtree':
            return method.CVTree_matrix_pairwise
        elif a_method == 'ma':
            return method.Ma_matrix_pairwise
        elif a_method == 'eu':
            return method.Eu_matrix_pairwise
        else:
            return method.d2_matrix_pairwise
 
def get_matrix_group(a_method):
    if a_method not in ['d2star', 'd2shepp', 'cvtree', 'ma', 'eu', 'd2']:
        print('Invalid method %s'%a_method)
        print('Only d2star,d2shepp,CVtree,Ma,Eu,d2 are supported')
        raise NameError
    else:
        if a_method == 'd2star':
            return method.d2star_matrix_groupwise
        elif a_method == 'd2shepp':
            return method.d2shepp_matrix_groupwise
        elif a_method == 'cvtree':
            return method.CVTree_matrix_groupwise
        elif a_method == 'ma':
            return method.Ma_matrix_groupwise
        elif a_method == 'eu':
            return method.Eu_matrix_groupwise
        else:
            return method.d2_matrix_groupwise

def write_phylip(output, a_method, sequence_list, matrix):
    if output.endswith('/'):
        filename = output + a_method + '.' + 'phylip'
    else:
        filename = '.'.join([output, a_method, 'phylip'])
    num = len(sequence_list)
    with open(filename, 'wt') as f:
        f.write('%d\n'%num)
        for i in range(num):
            f.write('.'.join(os.path.basename(sequence_list[i]).split('.')[:-1]))
            for value in matrix[i]:
                f.write('\t%.4f'%value)
            f.write('\n')

def write_plain(output, a_method, sequence_list, matrix):
    if output.endswith('/'):
        filename = output + a_method + '.' + 'plain'
    else:
        filename = '.'.join([output, a_method, 'plain'])
    num = len(sequence_list)
    with open(filename, 'wt') as f:
        for i in range(num):
            seq_1 = '.'.join(os.path.basename(sequence_list[i]).split('.')[:-1])
            for j in np.argsort(matrix[i]):
                seq_2 = '.'.join(os.path.basename(sequence_list[j]).split('.')[:-1])
                if i != j:
                    f.write('%s\t%s\t%.4f\n'%(seq_1, seq_2, matrix[i][j]))

def write_phylip_group(output, a_method, sequence_list_1, sequence_list_2, matrix):
    if output.endswith('/'):
        filename = output + a_method + '.' + 'phylip'
    else:
        filename = '.'.join([output, a_method, 'phylip'])
    num1 = len(sequence_list_1)
    num2 = len(sequence_list_2)
    with open(filename, 'wt') as f:
        f.write('%d,%d\n'%(num1, num2))
        for j in range(num2):
            f.write('\t'+'.'.join(os.path.basename(sequence_list_2[j]).split('.')[:-1])) 
        f.write('\n')
        for i in range(num1):
            f.write('.'.join(os.path.basename(sequence_list_1[i]).split('.')[:-1]))
            for value in matrix[i]:
                f.write('\t%.4f'%value)
            f.write('\n')

def write_plain_group(output, a_method, sequence_list_1, sequence_list_2, matrix):
    if output.endswith('/'):
        filename = output + a_method + '.' + 'plain'
    else:
        filename = '.'.join([output, a_method, 'plain'])
    num1 = len(sequence_list_1)
    with open(filename, 'wt') as f:
        for i in range(num1):
            seq_1 = '.'.join(os.path.basename(sequence_list_1[i]).split('.')[:-1])
            for j in np.argsort(matrix[i]):
                seq_2 = '.'.join(os.path.basename(sequence_list_2[j]).split('.')[:-1])
                f.write('%s\t%s\t%.4f\n'%(seq_1, seq_2, matrix[i][j]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Example: python alignmentfree.py -r -a d2star,d2shepp,CVtree -k 12 -m 10 -f filename -d dir -o output') 
    parser.add_argument('-a', dest='method', required = True, help='A list of alignment-free method, separated by comma: d2star,d2shepp,CVtree,Ma,Eu,d2')
    parser.add_argument('-k', dest='K', required = True, type = int, help='Kmer length')
    parser.add_argument('-m', dest='M', type = int, default=0, help='Markovian Order, required for d2star, d2shepp and CVtree')
    parser.add_argument('-f', dest='filename', help='A file that lists the paths of all samples, cannot be used together with -f1, -f2')
    parser.add_argument('-f1', dest='filename1', help='A file that lists the paths of the first group of samples, must be used together with -f2, cannot be used together with -f')
    parser.add_argument('-f2', dest='filename2', help='A file that lists the paths of the second group of samples, must be used together with -f1, cannot be used together with -f')
    parser.add_argument('-d', dest='Dir', default='None', help='A directory that saves kmer count')
    parser.add_argument('-o', dest='output', help='Prefix of output (defualt: Current directory)', default='./')
    parser.add_argument('-t', dest='threads', type = int, default=1, help='Number of threads')
    parser.add_argument('-r', dest='reverse', action='store_const',
                    const=True, default=False, help='Count the reverse complement of kmers (default: False)')
    args = parser.parse_args()
    K = args.K 
    M = args.M + 1
    filename = args.filename
    filename1 = args.filename1
    filename2 = args.filename2
    Reverse = args.reverse
    P_dir = args.Dir
    if P_dir == 'None':
        print('Warning: Using -d option to save kmer counts in a directory can save you a lot of counting time.')
    Num_Threads = args.threads
    output = args.output
    if output == './':
        print('Warning: Using -o option to change output directory and prefix. Otherwise, output will be generated in the current directory.')
    methods = [x.strip().lower() for x in args.method.split(',')]
    if filename and not filename1 and not filename2:
        sequence_list = get_sequence_from_file(filename) 
        for a_method in methods:
            print('Calculating %s.'%a_method)
            matrix = get_matrix(a_method)(sequence_list, M, K, Num_Threads, Reverse, P_dir)
            write_plain(output, a_method, sequence_list, matrix)
            write_phylip(output, a_method, sequence_list, matrix)
    elif not filename and filename1 and filename2:
        sequence_list_1 = get_sequence_from_file(filename1)
        sequence_list_2 = get_sequence_from_file(filename2)
        for a_method in methods:
            print('Calculating %s.'%a_method)
            matrix = get_matrix_group(a_method)(sequence_list_1, sequence_list_2, M, K, Num_Threads, Reverse, P_dir)
            write_phylip_group(output, a_method, sequence_list_1, sequence_list_2, matrix)
            write_plain_group(output, a_method, sequence_list_1, sequence_list_2, matrix)
    else:
        print('Cannot use -f, -f1, -f2 at the same time!')
