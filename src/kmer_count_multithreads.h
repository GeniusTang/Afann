#include <iostream>
#include <unordered_map>
#include "ctpl_stl.h"
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <thread>
#include <unistd.h>
#include <atomic>
#include <time.h>


std::atomic<int> X;
bool VALID = true;
std::unordered_map<char, int> nuc2num = {
        {'N', -1}, {'A', 0}, {'C', 1 }, {'G', 2}, {'T', 3},
	{'n', -1}, {'a', 0}, {'c', 1}, {'g', 2}, {'t', 3}, {'$', -2}
};

int revcomp(int num, int K){
    int nuc_rc = 0;
    int shift;
    for (int i=0; i<K; i++){
        shift = 2 * (K-i-1);
        nuc_rc += (3 - ((num>>shift)&3)) * pow(4, i);
    }
    return nuc_rc;
}


void count_one_read(int id, int K, std::string one_read, std::vector<std::atomic<int>> &count_array, bool Reverse) {
    int length = one_read.length();
    int mask = pow(2, (2*(K-1)))-1;
    int num = 0;
    int nuc_num = 0;
    std::unordered_map<char, int>::iterator search;
    int j = 0;
    char nuc;
    int rev = 0;
    for (int i=0;i<length;i++){ 
        if (!::VALID) {
            return;
        }
        nuc = one_read[i];
        search = nuc2num.find(nuc);
        if (search == nuc2num.end()){
            ::VALID = false;
            return;
        }
        else {
            nuc_num = search -> second;
        }
        if (nuc_num == -1){
            num = 0;
	    rev = 0;
            j = 0;
	}
        else{
            if (j < (K-1)){
                num = num * 4 + nuc_num;
		if (Reverse) rev = (rev >> 2) + (3-nuc_num) * pow(4, K-1);
                j += 1;
	    }
            else{
                num = (num&mask)<<2;
                num += nuc_num;
		count_array[num] ++;
		if (Reverse) {
		    /*
		    rev = revcomp(num, K);
		    count_array[rev] ++;
		    */
	            rev = (rev >> 2) + (3-nuc_num) * pow(4, K-1);
		    count_array[rev] ++;
		}
	   }
        }
    }
}

void count_one_read_M_K(int id, int M, int K, std::string one_read, std::vector<std::atomic<int>> &count_array, bool Reverse) {
    int length = one_read.length();
    int mask_K = pow(2, (2*(K-1)))-1;
    int mask_M = pow(2, 2*M)-1;
    int num_K = 0;
    int num_M = 0;
    int nuc_num = 0;
    std::unordered_map<char, int>::iterator search;
    int j = 0;
    int i = 0;
    int M_start = 0;
    int rev_K = 0;
    char nuc;
    nuc = one_read[0];
    search = nuc2num.find(nuc);
    if (search == nuc2num.end()){
        ::VALID = false;
        return;
    }
    if (search -> second == -2){
        i = 1;
	M_start = K-M;
    }
    for (;i<length;i++){
        if (!::VALID) {
            return;
        }
        nuc = one_read[i];
        search = nuc2num.find(nuc);
        if (search == nuc2num.end()){
            ::VALID = false;
            return;
        }
        else {
            nuc_num = search -> second;
        }
        if (nuc_num == -1){
            num_K = 0;
	    rev_K = 0;
            j = 0;
	    M_start = 0;
        }
        else{
	    if (j < (M-1)){
		num_K = num_K * 4 + nuc_num;
		if (Reverse) rev_K = (rev_K >> 2) + (3-nuc_num) * pow(4, K-1);
		j += 1;
	    }
	    else {
                if (j < (K-1)){
                    num_K = num_K * 4 + nuc_num;
		    if (Reverse) rev_K = (rev_K >> 2) + (3-nuc_num) * pow(4, K-1);
                    j += 1;
                }
                else{
                    num_K = (num_K&mask_K)<<2;
                    num_K += nuc_num;
                    count_array[num_K+mask_M+1] ++;
                    if (Reverse) {
			rev_K = (rev_K >> 2) + (3-nuc_num) * pow(4, K-1);
			count_array[rev_K+mask_M+1]++;
                    }
               }
	       if (M_start == 0) {
	           num_M = num_K&mask_M;
	           count_array[num_M] ++;
	           if (Reverse) count_array[(rev_K)>>(2*(K-M))]++;
	       }
	       else M_start--;
	    }
        }
    }
}

std::vector<std::atomic<int>> count(std::string filename, int K, int Num_Threads, bool Reverse) {
    const int SIZE = pow(4, K);
    std::vector<std::atomic<int>> count_array(SIZE);
    std::ifstream fs(filename);
    char temp_char[5000];
    std::string temp;
    std::string one_read;
    //int Num_Threads =  std::thread::hardware_concurrency();
    ctpl::thread_pool p(Num_Threads);
    while (::VALID) {
	fs.getline(temp_char, 5000, '\n');
	std::string temp(temp_char);
	if(temp[0] == '>'){
	    one_read = one_read + "N";
	}
	else{
	    one_read = one_read + temp;
	}
	if (one_read.length() >= 5000){
	    p.push([one_read, K, Reverse, &count_array](int id){count_one_read(id, K, one_read, count_array, Reverse);});
	    one_read = one_read.substr(one_read.length()-K+1);
        }
	if (fs.eof()) break;
	fs.clear();
    }
    p.push([one_read, K, Reverse, &count_array](int id){count_one_read(id, K, one_read, count_array, Reverse);});
    fs.close();
    if (!::VALID) {
        count_array[0] = -1;
    }
    return count_array;
}

std::vector<std::atomic<int>> count_M_K(std::string filename, int M, int K, int Num_Threads, bool Reverse) {
    const int SIZE = pow(4, M) + pow(4, K);
    const unsigned int READ_LENGTH = 5000;
    std::vector<std::atomic<int>> count_array(SIZE);
    std::ifstream fs(filename);
    char temp_char[READ_LENGTH];
    std::string temp;
    std::string one_read;
    //int Num_Threads =  std::thread::hardware_concurrency();
    ctpl::thread_pool p(Num_Threads);
    while (::VALID) {
        fs.getline(temp_char, READ_LENGTH, '\n');
        std::string temp(temp_char);
        if(temp[0] == '>'){
            one_read = one_read + "N";
        }
        else{
            one_read = one_read + temp;
        }
        if (one_read.length() >= READ_LENGTH){
            p.push([one_read, M, K, Reverse, &count_array](int id){count_one_read_M_K(id, M, K, one_read, count_array, Reverse);});
            one_read = "$" + one_read.substr(one_read.length()-K+1);
        }
        if (fs.eof()) break;
        fs.clear();
    }
    p.push([one_read, M, K, Reverse, &count_array](int id){count_one_read_M_K(id, M, K, one_read, count_array, Reverse);});
    fs.close();
    if (!::VALID) count_array[0] = -1;
    return count_array;
}
