#include <iostream>

using namespace std;

class Complex {
    public :
        float real, image;

    public :
    Complex(float real = 0, float image = 0) : real(real), image(image) {}
    Complex operator + (const Complex &o) {
        return Complex(this->real + o.real, this->image + o.image);
    }
    Complex operator - (const Complex &o) {
        return Complex(this->real - o.real, this->image - o.image);
    }
    Complex operator * (const Complex &o) {
        return Complex(this->real * o.real - this->image * o.image, this->real * o.image + this->image * o.real);
    }
    Complex & operator += (const Complex &o) {
        this->real += o.real;
        this->image += o.image;
        return *this;
    }
    Complex & operator /= (const int &o) {
        this->real /= o;
        this->image /= o;
        return *this;
    }
};

#define debug(...) cout << " [-] ", _dbg(#__VA_ARGS__, __VA_ARGS__)
template<class TH> void _dbg(const char *sdbg, TH h){ cout << sdbg << '=' << h << endl; }
template<class TH, class... TA> void _dbg(const char *sdbg, TH h, TA... a) {
    while(*sdbg != ',') cout << *sdbg++;
    cout << '=' << (h) << ','; 
    _dbg(sdbg+1, a...);
}

ostream & operator << (ostream &os, Complex x){
    cout << '{'
    << x.real << ", "
    << x.image << '}';
    return os;
}

#define SUBCARRIER_NUM 600
#define MAX_STEP 10

class Complex twiddle_factor[SUBCARRIER_NUM + 1][SUBCARRIER_NUM];
class Complex twiddle_factor_inverse[SUBCARRIER_NUM + 1][SUBCARRIER_NUM];
void precompute_twiddle_factor();

float precomputed_sin[SUBCARRIER_NUM + 1][SUBCARRIER_NUM + 1];
float precomputed_cos[SUBCARRIER_NUM + 1][SUBCARRIER_NUM];
void precompute_sincos();

int precomputed_permutation[SUBCARRIER_NUM + 1][SUBCARRIER_NUM];
void precompute_permutation();

int step_number[SUBCARRIER_NUM + 1];
int step_radix[SUBCARRIER_NUM + 1][MAX_STEP];
void precompute_step_number();

int radices[] = {59, 57, 53, 47, 43, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3, 2};

void precompute_all();
void precomputed_fft(int size, Complex *list);