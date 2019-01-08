from src._count import kmer_count
from src._count import kmer_count_seq
from src._count import kmer_count_m_k
from src._count import kmer_count_m_k_seq
from scipy.spatial.distance import cosine
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics.pairwise import manhattan_distances
from sklearn.metrics.pairwise import euclidean_distances
from scipy import stats
from functools import partial
import numpy as np
import os
from numpy import linalg as LA

Suffix = ['fna', 'fa', 'fasta']
Alphabeta = ['A', 'C', 'G', 'T']
Alpha_dict = dict(zip(Alphabeta, range(4)))

def rev_comp(num, K):
    nuc_rc = 0
    for i in range(K):
        shift = 2 * (K-i-1)
        nuc_rc += (3 - (num>>shift)&3) * (4**i)
    return nuc_rc

def rev_count(count, K):
    new_count = np.zeros_like(count)
    for i in range(4**K):
        rev = rev_comp(i, K)
        new_count[i] = count[i] + count[rev]
    return new_count

def check_count(seqfile, K_count):
    if K_count[0] == -1:
        raise Exception('Sequence file %s is not in the correct fasta format!'%seqfile)
    if np.sum(K_count) == 0:
        raise Exception('Sequence file %s is empty!'%seqfile)

def count_pickle(seqfile, K, Reverse, P_dir):
    seq_count_p = os.path.join(P_dir, os.path.basename(seqfile) + '.%s_%d_cnt.npy'%('R' if Reverse else 'NR', K))
    return seq_count_p

def get_sequences(seqfile):
    seq_old_name_list = []
    seq_new_name_list = []
    sequence_list = []
    sequence = ''
    first = True
    with open(seqfile) as f:
        for line in f.readlines():
            if line.startswith('>'):
                seq_old_name = line.strip().split()[0][1:]
                seq_new_name = seq_old_name.replace('/', '_slash_')
                seq_old_name_list.append(seq_old_name)
                seq_new_name_list.append(seq_new_name)
                if first:
                    first = False
                else:
                    sequence_list.append(sequence)
                    sequence = ''
            else:
                sequence += line.strip()
    sequence_list.append(sequence)
    return seq_old_name_list, seq_new_name_list, sequence_list

def get_K(seqfile, K, Num_Threads, Reverse, P_dir, sequence = '', from_seq=False):
    seq_count_K_p = count_pickle(seqfile, K, Reverse, P_dir)
    if os.path.exists(seq_count_K_p):
        K_count = np.load(seq_count_K_p)
    else:
        #print('Counting kmers of %s.'%seqfile)
        if not Reverse or K>= 6:
            if from_seq:
                K_count = np.copy(kmer_count_seq(sequence, K, Num_Threads, Reverse))
            else:
                K_count = np.copy(kmer_count(seqfile, K, Num_Threads, Reverse))
            check_count(seqfile, K_count)
        else:
            if from_seq:
                K_count = np.copy(kmer_count_seq(sequence, K, Num_Threads, False))
            else:
                K_count = np.copy(kmer_count(seqfile, K, Num_Threads, False))
            check_count(seqfile, K_count)
            K_count = rev_count(K_count, K)   
        if P_dir != 'None':
            np.save(seq_count_K_p, K_count)
    return K_count

def get_M_K(seqfile, M, K, Num_Threads, Reverse, P_dir, sequence = '', from_seq=False):
    if M >= K:
        raise ValueError('Markovian order cannot be greater than K-2!') 
    seq_count_M_p = count_pickle(seqfile, M, Reverse, P_dir)
    seq_count_K_p = count_pickle(seqfile, K, Reverse, P_dir)
    if os.path.exists(seq_count_M_p) and os.path.exists(seq_count_K_p):
        M_count = np.load(seq_count_M_p)
        K_count = np.load(seq_count_K_p)
    else:
        print('Counting kmers of %s.'%seqfile)
        if not Reverse or M>=6:
            if from_seq:
                count = np.copy(kmer_count_m_k_seq(sequence, M, K, Num_Threads, Reverse))
            else:
                count = np.copy(kmer_count_m_k(seqfile, M, K, Num_Threads, Reverse))
            check_count(seqfile, count) 
            M_count = count[:4**M]
            K_count = count[4**M:]
        else:
            if from_seq:
                M_count = np.copy(kmer_count_seq(sequence, M, Num_Threads, False))
            else:
                M_count = np.copy(kmer_count(seqfile, M, Num_Threads, False))
            check_count(seqfile, M_count)
            M_count = rev_count(M_count, M)
            if K>= 6:
                if from_seq:
                    K_count = np.copy(kmer_count_seq(sequence, K, Num_Threads, Reverse))
                else:
                    K_count = np.copy(kmer_count(seqfile, K, Num_Threads, Reverse))
            else:
                if from_seq:
                    K_count = np.copy(kmer_count_seq(sequence, K, Num_Threads, False)) 
                else:
                    K_count = np.copy(kmer_count(seqfile, K, Num_Threads, False))
                K_count = rev_count(K_count, K)
        if P_dir != 'None':
            np.save(seq_count_M_p, M_count)
            np.save(seq_count_K_p, K_count)
    return M_count, K_count

def get_transition(count_array):
    shape = len(count_array)
    transition_array = count_array.reshape(shape//4, 4)
    with np.errstate(divide='ignore', invalid='ignore'):
        transition_array = (transition_array / np.sum(transition_array, 1)[:, np.newaxis])
        transition_array[np.isnan(transition_array)] = 0
    return transition_array

def get_expect(seqfile, M, K, Num_Threads, Reverse, P_dir, sequence = '', from_seq=False):
    M_count, K_count = get_M_K(seqfile, M, K, Num_Threads, Reverse, P_dir, sequence, from_seq)
    seqfile_e_p = os.path.join(P_dir, os.path.basename(seqfile) + '.%s_M%d_K%d_e.npy'%('R' if Reverse else 'NR', M-1, K))
    if os.path.exists(seqfile_e_p):
        expect = np.load(seqfile_e_p)
    else:
        trans = get_transition(M_count)
        expect = M_count
        for _ in range(K-M):
            trans = np.tile(trans, (4, 1))
            expect = (expect[:,np.newaxis] * trans).flatten()
        if P_dir != 'None':
            np.save(seqfile_e_p, expect)
    return K_count, expect

def get_expect_reverse(seqfile, M, K, Num_Threads, P_dir):
    a_M_count, a_K_count = get_M_K(seqfile, M, K, Num_Threads, False, P_dir)
    b_M_count, b_K_count = get_M_K(seqfile, M, K, Num_Threads, True, P_dir)
    M_count = b_M_count - a_M_count
    del a_M_count
    del b_M_count
    K_count = b_K_count - a_K_count
    del a_K_count
    del b_K_count
    trans = get_transition(M_count)
    expect = M_count
    for _ in range(K-M):
        trans = np.tile(trans, (4, 1))
        expect = (expect[:,np.newaxis] * trans).flatten()
    return K_count, expect

def BIC(seqfile, K, Num_Threads, Reverse, P_dir, sequence = '', from_seq=False):
    M = K - 2
    M_count = get_K(seqfile, M+1, Num_Threads, Reverse, P_dir, sequence, from_seq)
    S = []
    for i in range(M+1):
        M_count = M_count.reshape(4**M, 4)
        with np.errstate(divide='ignore', invalid='ignore'):
            log = M_count * np.log(M_count/np.sum(M_count, axis=1)[:,np.newaxis])
            log[np.isnan(log)] = 0
            log = np.sum(log)
            bic = -2 * log + 3 * 4**M * np.log(np.sum(M_count))
        S = [bic] + S
        M_count = np.sum(M_count, axis=1)
        M -= 1
    return S.index(min(S))

def all_BIC(seqname_list, K, Num_Threads, Reverse, P_dir, sequence_list = [], from_seq=False):
    order = []
    if from_seq:
        for seqname, sequence in zip(seqname_list, sequence_list):
            order.append(BIC(seqname, K, Num_Threads, Reverse, P_dir, sequence, from_seq))
    else:
        for seqfile in seqname_list:
            order.append(BIC(seqfile, K, Num_Threads, Reverse, P_dir))
    return order

def get_d2star_f(seqfile, M, K, Num_Threads, Reverse, P_dir, sequence = '', from_seq=False):
    seqfile_f_p = os.path.join(P_dir, os.path.basename(seqfile) + '.%s_M%d_K%d_d2star_f.npy'%('R' if Reverse else 'NR', M-1, K))
    if os.path.exists(seqfile_f_p):
        d2star_f = np.load(seqfile_f_p)
    else:
        K_count, expect = get_expect(seqfile, M, K, Num_Threads, Reverse, P_dir, sequence, from_seq)
        with np.errstate(divide='ignore', invalid='ignore'):
            d2star_f = (K_count-expect)/np.sqrt(expect)
            d2star_f[np.isnan(d2star_f)]=0
        if P_dir != 'None':
            np.save(seqfile_f_p, d2star_f)
    return d2star_f

def get_d2star_all_f(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list = [], from_seq=False):
    N = len(seqname_list)
    f_matrix = np.ones((N, 4**K))
    for i in range(N):
        if from_seq:
            sequence = sequence_list[i]
        else:
            sequence = ''
        seqfile = seqname_list[i]
        f_matrix[i] = get_d2star_f(seqfile, M, K, Num_Threads, Reverse, P_dir, sequence, from_seq)
    return f_matrix

def get_CVTree_f(seqfile, M, K, Num_Threads, Reverse, P_dir, sequence = '', from_seq=False):
    M = K - 1
    seqfile_f_p = os.path.join(P_dir, os.path.basename(seqfile) + '.%s_M%d_K%d_CVTree_f.npy'%('R' if Reverse else 'NR', M-1, K))
    if os.path.exists(seqfile_f_p):
        CVTree_f = np.load(seqfile_f_p)
    else:
        K_count, expect = get_expect(seqfile, M, K, Num_Threads, Reverse, P_dir, sequence, from_seq)
        with np.errstate(divide='ignore', invalid='ignore'):
            CVTree_f = (K_count-expect)/expect
            CVTree_f[np.isnan(CVTree_f)]=0
        if P_dir != 'None':
            np.save(seqfile_f_p, CVTree_f)
    return CVTree_f   

def get_CVTree_all_f(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list = [], from_seq=False):
    N = len(seqname_list)
    f_matrix = np.ones((N, 4**K))
    for i in range(N):
        if from_seq:
            sequence = sequence_list[i]
        else:
            sequence = ''
        seqfile = seqname_list[i]
        f_matrix[i] = get_CVTree_f(seqfile, M, K, Num_Threads, Reverse, P_dir, sequence, from_seq)
    return f_matrix

def get_all_f(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list = [], from_seq=False):
    N = len(seqname_list)
    f_matrix = np.ones((N, 4**K))
    for i in range(N):
        if from_seq:
            sequence = sequence_list[i]
        else:
            sequence = ''
        seqfile = seqname_list[i]
        a_K = get_K(seqfile, K, Num_Threads, Reverse, P_dir, sequence, from_seq)
        f_matrix[i] = a_K/np.sum(a_K)
    return f_matrix

def get_all_diff(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list = [], from_seq=False):
    N = len(seqname_list)
    diff_matrix = np.ones((N, 4**K))
    for i in range(N):
        if from_seq:
            sequence = sequence_list[i]
        else:
            sequence = ''
        seqfile = seqname_list[i]
        a_K_count, a_expect = get_expect(seqfile, M, K, Num_Threads, Reverse, P_dir, sequence, from_seq)
        a_diff = a_K_count - a_expect
        diff_matrix[i] = a_diff
    return diff_matrix

'''
def get_all_K(sequence_list, M, K, Num_Threads, Reverse, P_dir):
    K_matrix = np.ones((len(sequence_list), 4**K))
    for i, seqfile in enumerate(sequence_list):
        M_count, K_count = get_M_K(seqfile, M, K, Num_Threads, Reverse, P_dir)
        K_matrix[i] = K_count/np.sum(K_count)
    return K_matrix
'''

def Ma(seqfile_1, seqfile_2, M, K, Num_Threads, Reverse, P_dir, sequence_1 = '', sequence_2 = '', from_seq=False):
    a_K = get_K(seqfile_1, K, Num_Threads, Reverse, P_dir, sequence_1, from_seq)
    b_K = get_K(seqfile_2, K, Num_Threads, Reverse, P_dir, sequence_2, from_seq) 
    diff = a_K/np.sum(a_K) - b_K/np.sum(b_K)
    return LA.norm(diff, 1) 

def Eu(seqfile_1, seqfile_2, M, K, Num_Threads, Reverse, P_dir, sequence_1 = '', sequence_2 = '', from_seq=False):
    a_K = get_K(seqfile_1, K, Num_Threads, Reverse, P_dir, sequence_1, from_seq)
    b_K = get_K(seqfile_2, K, Num_Threads, Reverse, P_dir, sequence_2, from_seq)
    diff = a_K/np.sum(a_K) - b_K/np.sum(b_K)
    return LA.norm(diff)

def d2(seqfile_1, seqfile_2, M, K, Num_Threads, Reverse, P_dir, sequence_1 = '', sequence_2 = '', from_seq=False):
    a_K = get_K(seqfile_1, K, Num_Threads, Reverse, P_dir, sequence_1, from_seq)
    b_K = get_K(seqfile_2, K, Num_Threads, Reverse, P_dir, sequence_2, from_seq)
    return 0.5 * cosine(a_K/np.sum(a_K), b_K/np.sum(b_K))

def d2star(seqfile_1, seqfile_2, M, K, Num_Threads, Reverse, P_dir, sequence_1 = '', sequence_2 = '', from_seq=False):
    if from_seq:
        a_f = get_d2star_f(seqfile_1, M, K, Num_Threads, Reverse, P_dir, sequence_1, from_seq)
        b_f = get_d2star_f(seqfile_2, M, K, Num_Threads, Reverse, P_dir, sequence_2, from_seq)
    else:
        a_f = get_d2star_f(seqfile_1, M, K, Num_Threads, Reverse, P_dir)
        b_f = get_d2star_f(seqfile_2, M, K, Num_Threads, Reverse, P_dir)
    return 0.5 * cosine(a_f, b_f)

def CVTree(seqfile_1, seqfile_2, M, K, Num_Threads, Reverse, P_dir, sequence_1 = '', sequence_2 = '', from_seq=False):
    if from_seq:
        a_f = get_CVTree_f(seqfile_1, M, K, Num_Threads, Reverse, P_dir, sequence_1, from_seq)
        b_f = get_CVTree_f(seqfile_2, M, K, Num_Threads, Reverse, P_dir, sequence_2, from_seq)
    else:
        a_f = get_CVTree_f(seqfile_1, M, K, Num_Threads, Reverse, P_dir)
        b_f = get_CVTree_f(seqfile_2, M, K, Num_Threads, Reverse, P_dir)
    return 0.5 * cosine(a_f, b_f)

def d2shepp(seqfile_1, seqfile_2, M, K, Num_Threads, Reverse, P_dir, sequence_1 = '', sequence_2 = '', from_seq=False):
    a_K_count, a_expect = get_expect(seqfile_1, M, K, Num_Threads, Reverse, P_dir, sequence_1, from_seq)
    a_diff = a_K_count - a_expect
    del a_K_count
    del a_expect
    b_K_count, b_expect = get_expect(seqfile_2, M, K, Num_Threads, Reverse, P_dir, sequence_2, from_seq)
    b_diff = b_K_count - b_expect
    del b_K_count
    del b_expect
    denom = np.power(a_diff**2 + b_diff**2, 0.25)
    with np.errstate(divide='ignore', invalid='ignore'):
        a_f = a_diff/denom
        a_f[np.isnan(a_f)]=0
        b_f = b_diff/denom
        b_f[np.isnan(b_f)]=0 
    return 0.5 * cosine(a_f, b_f)

def d2shepp_error(seqfile, M, K, Num_Threads, P_dir):
    a_K_count, a_expect = get_expect(seqfile, M, K, Num_Threads, False, P_dir)
    a_diff = a_K_count - a_expect
    del a_K_count
    del a_expect
    b_K_count, b_expect = get_expect_reverse(seqfile, M, K, Num_Threads, P_dir)
    b_diff = b_K_count - b_expect
    del b_K_count
    del b_expect
    denom = np.power(a_diff**2 + b_diff**2, 0.25)
    with np.errstate(divide='ignore', invalid='ignore'):
        a_f = a_diff/denom
        a_f[np.isnan(a_f)]=0
        b_f = b_diff/denom
        b_f[np.isnan(b_f)]=0
    return 0.5 * cosine(a_f, b_f)

def cosine_matrix(f1_matrix, f2_matrix=None):
    if f2_matrix:
        matrix = 0.5 * (1 - cosine_similarity(f1_matrix, f2_matrix))
    else:
        matrix = 0.5 * (1 - cosine_similarity(f1_matrix))
    np.fill_diagonal(matrix, 0)
    return matrix

def Ma_matrix(f1_matrix, f2_matrix=None):
    if f2_matrix:
        matrix = manhattan_distances(f1_matrix, f2_matrix)
    else:
        matrix = manhattan_distances(f1_matrix)
    np.fill_diagonal(matrix, 0)
    return matrix

def Eu_matrix(f1_matrix, f2_matrix=None):
    if f2_matrix:
        matrix = euclidean_distances(f1_matrix, f2_matrix)
    else:
        matrix = euclidean_distances(f1_matrix)
    np.fill_diagonal(matrix, 0)
    return matrix
 
def dist_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list = [], from_seq=False, method = None):
    #print('Slow mode')
    N = len(seqname_list)
    matrix = np.zeros((N, N))
    sequence_1 = ''
    sequence_2 = ''
    for i in range(N):
        if from_seq:
            sequence_1 = sequence_list[i]
        seqfile_1 = seqname_list[i]
        for j in range(i+1, N):
            if from_seq:
                sequence_2 = sequence_list[j]
            seqfile_2 = seqname_list[j]
            matrix[i][j] = method(seqfile_1, seqfile_2, M, K, Num_Threads, Reverse, P_dir, sequence_1, sequence_2, from_seq)
            matrix[j][i] = matrix[i][j]
    return matrix

def d2star_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list = [], from_seq=False, slow=False):
    if not slow:
        f_matrix = get_d2star_all_f(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list, from_seq)
        return cosine_matrix(f_matrix)
    else:
        return dist_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list, from_seq, method = d2star)

def CVTree_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list = [], from_seq=False, slow=False):
    if not slow:
        f_matrix = get_CVTree_all_f(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list, from_seq)
        return cosine_matrix(f_matrix)
    else:
        return dist_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list, from_seq, method = CVTree)

def d2_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list = [], from_seq=False, slow=False):
    if not slow:
        f_matrix = get_all_f(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list, from_seq)
        return cosine_matrix(f_matrix)
    else:
        return dist_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list, from_seq, method = d2)

def Ma_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list = [], from_seq=False, slow=False):
    if not slow:
        f_matrix = get_all_f(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list, from_seq)
        return Ma_matrix(f_matrix)
    else:
        return dist_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list, from_seq, method = Ma)

def Eu_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list = [], from_seq=False, slow=False):
    if not slow:
        f_matrix = get_all_f(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list, from_seq)
        return Eu_matrix(f_matrix)
    else:
        return dist_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list, from_seq, method = Eu)

def d2shepp_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list = [], from_seq=False, slow=False):
    if not slow:
        diff_matrix = get_all_diff(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list, from_seq) 
        N = len(seqname_list)
        matrix = np.zeros((N, N))  
        for i in range(N):
            a_diff = diff_matrix[i]
            for j in range(i+1, N):
                b_diff = diff_matrix[j]
                denom = np.power(a_diff**2 + b_diff**2, 0.25)
                with np.errstate(divide='ignore', invalid='ignore'):
                    a_f = a_diff/denom
                    a_f[np.isnan(a_f)]=0
                    b_f = b_diff/denom
                    b_f[np.isnan(b_f)]=0 
                matrix[i][j] = 0.5 * cosine(a_f, b_f)
                matrix[j][i] = matrix[i][j]
        return matrix 
    else:
        return dist_matrix_pairwise(seqname_list, M, K, Num_Threads, Reverse, P_dir, sequence_list, from_seq, method = d2shepp)
 
def dist_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1 = [], sequence_list_2 = [], from_seq=False, method=None):
    #print('Slow mode')
    N1 = len(seqname_list_1)
    N2 = len(seqname_list_2)
    matrix = np.zeros((N1, N2))
    sequence_1 = ''
    sequence_2 = ''
    for i in range(N1):
        if from_seq:
            sequence_1 = sequence_list_1[i]
        seqfile_1 = seqname_list_1[i]
        for j in range(N2):
            if from_seq:
                sequence_2 = sequence_list_2[j]
            seqfile_2 = seqname_list_2[j]
            matrix[i][j] = method(seqfile_1, seqfile_2, M, K, Num_Threads, Reverse, P_dir, sequence_1, sequence_2, from_seq)
    return matrix

def d2star_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1 = [], sequence_list_2 = [], from_seq=False, slow=False):
    if not slow:
        f1_matrix = get_d2star_all_f(seqname_list_1, M, K, Num_Threads, Reverse, P_dir, sequence_list_1, from_seq)
        f2_matrix = get_d2star_all_f(seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_2, from_seq)
        return cosine_matrix(f1_matrix, f2_matrix)
    else:
        return dist_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1, sequence_list_2, from_seq, method = d2star)

def CVTree_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1 = [], sequence_list_2 = [], from_seq=False, slow=False):
    if not slow:
        f1_matrix = get_CVTree_all_f(seqname_list_1, M, K, Num_Threads, Reverse, P_dir, sequence_list_1, from_seq)
        f2_matrix = get_CVTree_all_f(seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_2, from_seq)
        return cosine_matrix(f1_matrix, f2_matrix)
    else:
        return dist_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1, sequence_list_2, from_seq, method = CVTree)

def d2_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1 = [], sequence_list_2 = [], from_seq=False, slow=False):
    if not slow:
        f1_matrix = get_all_f(seqname_list_1, M, K, Num_Threads, Reverse, P_dir, sequence_list_1, from_seq)
        f2_matrix = get_all_f(seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_2, from_seq)
        return cosine_matrix(f1_matrix, f2_matrix)
    else:
        return dist_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1, sequence_list_2, from_seq, method = d2)

def Ma_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1 = [], sequence_list_2 = [], from_seq=False, slow=False):
    if not slow:
        f1_matrix = get_all_f(seqname_list_1, M, K, Num_Threads, Reverse, P_dir, sequence_list_1, from_seq)
        f2_matrix = get_all_f(seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_2, from_seq)
        return Ma_matrix(f1_matrix, f2_matrix)
    else:
        return dist_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1, sequence_list_2, from_seq, method = Ma)

def Eu_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1 = [], sequence_list_2 = [], from_seq=False, slow=False):
    if not slow:
        f1_matrix = get_all_f(seqname_list_1, M, K, Num_Threads, Reverse, P_dir, sequence_list_1, from_seq)
        f2_matrix = get_all_f(seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_2, from_seq)
        return Eu_matrix(f1_matrix, f2_matrix)
    else:
        return dist_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1, sequence_list_2, from_seq, method = Eu)

def d2shepp_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1 = [], sequence_list_2 = [], from_seq=False, slow=False):
    if not slow:
        a_diff_matrix = get_all_diff(seqname_list_1, M, K, Num_Threads, Reverse, P_dir, sequence_list_1, from_seq)
        b_diff_matrix = get_all_diff(seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_2, from_seq) 
        N1 = len(seqname_list_1)
        N2 = len(seqname_list_2)
        matrix = np.zeros((N1, N2))
        for i in range(N1):
            a_diff = a_diff_matrix[i]
            for j in range(N2):
                b_diff = b_diff_matrix[j]
                denom = np.power(a_diff**2 + b_diff**2, 0.25)
                with np.errstate(divide='ignore', invalid='ignore'):
                    a_f = a_diff/denom
                    a_f[np.isnan(a_f)]=0
                    b_f = b_diff/denom
                    b_f[np.isnan(b_f)]=0 
                matrix[i][j] = 0.5 * cosine(a_f, b_f)
        return matrix
    else:
        return dist_matrix_groupwise(seqname_list_1, seqname_list_2, M, K, Num_Threads, Reverse, P_dir, sequence_list_1, sequence_list_2, from_seq, method = d2shepp)

def error_array(sequence_list, M, K, Num_Threads, P_dir, method):
    N = len(sequence_list)
    array = np.zeros(N)
    for i in range(N):
        seqfile = sequence_list[i]
        array[i] = method(seqfile, M, K, Num_Thread, P_dir)
    return array

d2shepp_error_array = partial(error_array, method = d2shepp_error)

