#include "precomputedFFT.h"
#include "timer.h"
#include "/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/include/mkl_dfti.h"
#include "/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/include/mkl.h"

#include <math.h>

#include <iostream>
#include <assert.h>
using namespace std;

void precompute_sincos() {
    for(int i = 1 ; i <= SUBCARRIER_NUM ; i++) {
        for(int j = 0 ; j < SUBCARRIER_NUM ; j++) {
            float angle = 2 * M_PI * j / i;
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
            int base = j - j % step_size;
            int step_offset = j - base;
            int prev_step_offset = j % previous_size;
            for(int k = 0 ; k < radix ; k++) {
                int tmp_list_index = base + k * previous_size + prev_step_offset;
                list[j] += tmp_list[tmp_list_index] * twiddle_factor[step_size][step_offset * k % step_size];
                //debug(step_offset * k, SUBCARRIER_NUM);
                //assert(step_offset * k < SUBCARRIER_NUM);
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
            }
        }
    }

    for(int i = 0 ; i < size ; i++) list[i] /= size;
}
void validate_check(int size, Complex *arr) {
    Complex a[size];
    float real[size], image[size];

    double eps = 1e-2;
    for(int i = 0 ; i < size ; i++) {
        a[i] = arr[i];
        real[i] = arr[i].real;
        image[i] = arr[i].image;
    }

    precomputed_fft(size, a);

    DFTI_DESCRIPTOR_HANDLE my_desc2_handle = NULL;
    MKL_LONG status;

    status = DftiCreateDescriptor(&my_desc2_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, size);
    status = DftiSetValue(my_desc2_handle, DFTI_COMPLEX_STORAGE, DFTI_REAL_REAL);
    status = DftiCommitDescriptor(my_desc2_handle);

    status = DftiComputeForward(my_desc2_handle, real, image);

    for(int i = 0 ; i < size ; i++) {
        float diff_real = (a[i].real - real[i]);
        float diff_image = (a[i].image - image[i]);

        assert(abs(diff_real) < eps);
        assert(abs(diff_image) < eps);
    }
}

int sizes[] = {512};

#define REPEAT 200

int main() {
    Complex list[SUBCARRIER_NUM];
    float r2c_data_real[SUBCARRIER_NUM], r2c_data_image[SUBCARRIER_NUM];

    precompute_all();

    DFTI_DESCRIPTOR_HANDLE my_desc2_handle = NULL;
    MKL_LONG status;

    for(int size : sizes) {
        for(int i = 0 ; i < size ; i++) list[i].real = r2c_data_real[i] = rand() % 100;

        validate_check(size, list);

        DftiCreateDescriptor(&my_desc2_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, size);
        DftiSetValue(my_desc2_handle, DFTI_COMPLEX_STORAGE, DFTI_REAL_REAL);
        DftiCommitDescriptor(my_desc2_handle);

        Timer m1, m2;

        precomputed_fft(size, list);
        DftiComputeForward(my_desc2_handle, r2c_data_real, r2c_data_image);

        m1.start();
        for(int i = 0 ; i < REPEAT ; i++) precomputed_fft(size, list);
        m1.stop();

        m2.start();
        for(int i = 0 ; i < REPEAT ; i++) DftiComputeForward(my_desc2_handle, r2c_data_real, r2c_data_image);
        m2.stop();

        status = DftiFreeDescriptor(&my_desc2_handle);
        //precomputed_inverse_fft(size, list);

        cout << m1.elapsednanoseconds() / REPEAT << ' ' << m2.elapsednanoseconds() / REPEAT << '\n';
        cout << m1.elapsednanoseconds() / m2.elapsednanoseconds() << '\n';
    }

    return 0;
}