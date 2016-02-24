#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include "gmp.h"
#include "gmpxx.h"
using namespace std;
void Rand_x_modp (mpz_t x, mpz_t p, mpz_t n, unsigned long c, int m2exp); //rand() 0<x<p
int ZeroStep(mpz_t d, mpz_t p);
bool Miller_Rabin(mpz_t p, int k);
bool Step_2_point_2(mpz_t p, mpz_t x, mpz_t d, int s);
void Horner(mpz_t result, mpz_t x, mpz_t exp, mpz_t mod);
void Rand_x (mpz_t x, mpz_t n, unsigned long c, int m2exp, int min); //rand() from interval [2^min;2^(min+1)]
void decrypt(mpz_t M, mpz_t C, mpz_t d, mpz_t n);
void encrypt(mpz_t C, mpz_t M, mpz_t e, mpz_t n);
void signat(mpz_t S, mpz_t M, mpz_t d, mpz_t n);
void RSA(mpz_t e, mpz_t d,  mpz_t n, mpz_t phi, mpz_t e1, mpz_t d1, mpz_t n1, mpz_t phi1);




