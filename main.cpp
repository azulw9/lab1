#include "ascrypt2.h"

int main(int argc, char *argv[])
{
    mpz_t e, n, d, phi, e1, n1, d1, phi1, k,k1,S,S1,temp;
    mpz_init_set_str (temp, "157863", 0); //for Random func
    mpz_init_set_str(e1, "10001", 16);   mpz_init_set_str(e, "10001", 16);
    mpz_init(n); mpz_init(d);  mpz_init(phi);
    mpz_init(n1); mpz_init(d1);  mpz_init(phi1);
    mpz_init(k);  mpz_init(k1);
    mpz_init(S);  mpz_init(S1);
    RSA(e,d,n,phi,e1,d1,n1,phi1);
    FILE* key;
    key=fopen("KEY_EXCHANGE.txt", "w");
    fputs("A:\n public key\n e= ", key);
    mpz_out_str(key, 16,e);
    fputs("\n n= ", key);
    mpz_out_str(key,16,n);
    fputs("\n private key\n d= ", key);
    mpz_out_str(key,16,d);
    fputs("\n B:\n public key\n e1= ", key);
    mpz_out_str(key, 16,e1);
    fputs("\n n1= ", key);
    mpz_out_str(key,16,n1);
    fputs("\n private key\n d1= ", key);
    mpz_out_str(key,16,d1);
    Rand_x_modp(k,n,temp,2345,78);
    fputs("\n\n  k= ", key);
    mpz_out_str(key,16,k);
    Horner(k1,k,e1,n1);
    fputs("\n k1= ", key);
    mpz_out_str(key,16,k1);
    Horner(S,k,d,n);
    fputs("\n S= ", key);
    mpz_out_str(key,16,S);
    Horner(S1,S,e1,n1);
    fputs("\n S1= ", key);
    mpz_out_str(key,16,S1);
    Horner(k,k1,d1,n1);
    fputs("\n k'= ", key);
    mpz_out_str(key,16,k);
    Horner(S,S1,d1,n1);
    fputs("\n S' = ", key);
    mpz_out_str(key,16,S);
    Horner(k,S,e,n);
    fputs("\n Check signature of A:\n S^emodn= ", key);
    mpz_out_str(key,16,k);
    system("PAUSE");
    return EXIT_SUCCESS;
}
