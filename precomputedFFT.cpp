#include "precomputedFFT.h"
#include <math.h>

#include <iostream>
using namespace std;

void precompute_sincos() {
    for(int i = 1 ; i <= SUBCARRIER_NUM ; i++) {
        for(int j = 0 ; j < SUBCARRIER_NUM ; j++) {
            float angle = 2 * M_PI / i * j;
            precomputed_sin[i][j] = sin(angle);
            precomputed_cos[i][j] = cos(angle);
        }
    }
}

void precompute_twiddle_factor() {
    for(int i = 1 ; i <= SUBCARRIER_NUM ; i++){
        for(int j = 0 ; j < SUBCARRIER_NUM ; j++) {
            twiddle_factor[i][j] = Complex(precomputed_cos[i][j], -precomputed_sin[i][j]);
            twiddle_factor_inverse[i][j] = Complex(precomputed_cos[i][j], precomputed_sin[i][j]);
        }
    }
}

void precompute_step_number() {
    int radix_num = sizeof(radices) / sizeof(radices[0]);
    for(int i = 1 ; i <= SUBCARRIER_NUM ; i++) {
        int ii = i;
        for(int j = 0 ; j < radix_num ; j++){
            while(ii % radices[j] == 0) {
                ii /= radices[j];
                step_radix[i][step_number[i]] = radices[j];
                step_number[i]++;
            }
        }
    }
}

void precompute_permutation() {
    for(int i = 1 ; i <= SUBCARRIER_NUM ; i++) {
        for(int j = 0 ; j < SUBCARRIER_NUM ; j++) {
            int ii = i, jj = j;
            int radix_inverse = 0;
            for(int k = 0 ; k < step_number[i] ; k++) {
                radix_inverse += (jj % step_radix[i][k]) * (ii / step_radix[i][k]);
                jj /= step_radix[i][k];
                ii /= step_radix[i][k];
            }
            precomputed_permutation[i][j] = radix_inverse;
        }
    }
}

void precompute_all() {
    precompute_sincos();
    precompute_twiddle_factor();
    precompute_step_number();
    precompute_permutation();
}

void precomputed_fft(int size, Complex *list) {
    Complex tmp_list[size];
    for(int i = 0 ; i < size ; i++) tmp_list[i] = list[i];
    for(int i = 0 ; i < size ; i++) list[precomputed_permutation[size][i]] = tmp_list[i];

    int step_size = 1, previous_size = 1;
    for(int i = step_number[size] - 1 ; i >= 0 ; i--){
        int radix = step_radix[size][i];
        previous_size = step_size;
        step_size *= radix;

        for(int j = 0 ; j < size ; j++) {
            tmp_list[j] = list[j];
            list[j].real = 0;
            list[j].image = 0;
        }

        for(int j = 0 ; j < size ; j++) {
            int base = j / step_size * step_size;
            int step_offset = j - base;
            int prev_step_offset = j % previous_size;
            for(int k = 0 ; k < radix ; k++) {
                int tmp_list_index = base + k * previous_size + prev_step_offset;
                list[j] += tmp_list[tmp_list_index] * twiddle_factor[step_size][step_offset * k];
                //debug(j, tmp_list[tmp_list_index], twiddle_factor[step_size][k]);
            }
        }
    }
}

void precomputed_inverse_fft(int size, Complex *list) {
    Complex tmp_list[size];
    for(int i = 0 ; i < size ; i++) tmp_list[i] = list[i];
    for(int i = 0 ; i < size ; i++) list[precomputed_permutation[size][i]] = tmp_list[i];

    int step_size = 1, previous_size = 1;
    for(int i = step_number[size] - 1 ; i >= 0 ; i--){
        int radix = step_radix[size][i];
        previous_size = step_size;
        step_size *= radix;
        
        for(int j = 0 ; j < size ; j++) {
            tmp_list[j] = list[j];
            list[j].real = 0;
            list[j].image = 0;
        }

        for(int j = 0 ; j < size ; j++) {
            int base = j / step_size * step_size;
            int step_offset = j - base;
            int prev_step_offset = j % previous_size;
            for(int k = 0 ; k < radix ; k++) {
                int tmp_list_index = base + k * previous_size + prev_step_offset;
                list[j] += tmp_list[tmp_list_index] * twiddle_factor_inverse[step_size][step_offset * k];
                //debug(j, tmp_list[tmp_list_index], twiddle_factor[step_size][k]);
            }
        }
    }

    for(int i = 0 ; i < size ; i++) list[i] /= size;
}


int sizes[] = {8};


int main() {
    Complex list[SUBCARRIER_NUM];

    precompute_all();

    for(int size : sizes) {
        //for(int i = 0 ; i < size ; i++) cout << twiddle_factor[size][i] << ' ';
        //cout << '\n';

        for(int i = 0 ; i < size ; i++) list[i].real = 2;
        
        precomputed_fft(size, list);
        precomputed_inverse_fft(size, list);
    }

    return 0;
}