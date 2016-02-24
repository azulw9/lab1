#include "ascrypt2.h"

void Rand_x_modp (mpz_t x, mpz_t p, mpz_t n, unsigned long c, int m2exp) //c, m2exp ìåíÿòü ïåðåä íîâûì âûçîâîì
{
     mpz_t temp; mpz_init(temp);
     gmp_randstate_t state; 
     // gmp_randinit_lc_2exp(state, n, c, m2exp);
     mpz_urandomb(n,state,257);
	 // mpz_urandomb(n,state,256);
     
}

int ZeroStep(mpz_t d, mpz_t p) // p=p-1 ïåðåäàâàòü ñðàçó p-1
{
     mpz_t temp, p_1;
     mpz_init_set(p_1,p);     
     int s=0;
     mpz_init(temp);     
     while (mpz_divisible_ui_p(p_1,2)!=0)
     {
           mpz_cdiv_q_ui(temp,p_1,2);
           mpz_set(p_1,temp);
           s++;
     }
     mpz_set(d,p_1);
     return s;
 }
 
bool Step_2_point_2(mpz_t p, mpz_t x, mpz_t d, int s)
{
     mpz_t temp, exp, temp1, temp2;
     mpz_init(temp); mpz_init(exp); mpz_init(temp1); mpz_init(temp2);
     mpz_sub_ui(temp2, p, 1);
     int r=1;
     Horner(temp,x,d,p);
     if (mpz_cmp_ui(temp,1)==0) {return 1;}
     else if (mpz_cmp(temp,temp2)==0) {return 1;}
     else
     {
      while(r!=s)
      {
           if (r==1)
              {
                      mpz_mul_2exp(exp,d,r);
                      Horner(temp,x,exp,p);
                      mpz_set(temp1,temp); 
              }
          else
          {
           Horner(temp,temp1,exp,p);
           mpz_set(temp1,temp);
          }
          if (mpz_cmp_ui(temp,1)==0) {return 0;}
          else if (mpz_cmp(temp,temp2)==0) {return 1;}
          else r++;             
     } }
     return 0;
 }

bool Miller_Rabin(mpz_t p, int k)
{
     int i=0, m2exp=0, s=0, check=0;
     unsigned long c = 0;
     mpz_t temp, p_1, d, x, n, temp1; mpz_init(temp1); mpz_init(temp); mpz_init(p_1); mpz_init(d); mpz_init(x); mpz_init_set_str (n, "34", 0);
     mpz_sub_ui(p_1,p,1);
     s=ZeroStep(d,p_1);
     while (i<k)
     {
           m2exp=rand()%100+2;
           c = rand()%10000+1;
           Rand_x_modp (x, p, n, c, m2exp);
           mpz_gcd(temp1,x,p);
           if (mpz_cmp_ui(temp1,1)==0) 
           { 
               bool t = Step_2_point_2(p,x,d,s);
               if (t==false) return false;
               else i++;
           }
            else {return false;}
     }
     if (i==k) return true;
     else return false;
}

void Horner(mpz_t result, mpz_t x, mpz_t exp, mpz_t mod)
{
     mpz_t temp, temp1;
     mpz_init_set_str(temp, "0", 0);
     mpz_init_set_str(temp1, "0", 0);
     int size=(mpz_sizeinbase (exp, 2)-1); 
     mpz_set_str(result, "1", 0);
     for (int i=size; i>=0; i--)
     {
         mpz_powm_ui(temp,result,2,mod);
         if (mpz_tstbit(exp, i)==1) 
         {
                  mpz_mul(temp1, temp, x); 
                  mpz_mod(result, temp1, mod);
         }  
         else mpz_set(result,temp);
     }
 }
 
void Rand_x (mpz_t x, mpz_t n, unsigned long c, int m2exp, int min)
{
      int k=0;
      mpz_t temp; mpz_init(temp);
      mpz_t two; mpz_init_set_str(two,"2",0);
      mpz_t r_min; mpz_init(r_min); 
      mpz_pow_ui(r_min, two, min);
      gmp_randstate_t state; 
      gmp_randinit_lc_2exp(state, n, c, m2exp);
      mpz_urandomb(temp,state,min);
      mpz_add(x, temp, r_min);
      if (mpz_divisible_ui_p(x,2)==1) {mpz_add_ui(temp, x, 1); mpz_init_set(x,temp);}
      k=5+rand()%10;
      bool h=Miller_Rabin(x,k);
      while (h!=1) {mpz_add_ui(temp,x,2); mpz_set(x,temp); h=Miller_Rabin(x,k);}
}

void RSA(mpz_t e, mpz_t d,  mpz_t n, mpz_t phi, mpz_t e1, mpz_t d1, mpz_t n1, mpz_t phi1)
{
     FILE* rsa;
     string s;
     rsa = fopen("RSA.txt", "w");
     mpz_t p, q, p1, q1, temp1, temp2; mpz_init(p); mpz_init(q); mpz_init(p1); mpz_init(q1); mpz_init(temp1); mpz_init(temp2);  
     unsigned long c=0;
     int mp2exp=0;
     mpz_set_str(e,"65537",0);
     mpz_set_str(e1,"65537",0);
     mpz_set_str(n,"345678908765432345678654",0);
     cout<<"Input randomize constants\n C and mp2exp\n";
     cin>>c;
     cin>>mp2exp;
     Rand_x(p,n,c,mp2exp,256);
     Rand_x(q,n,c+123,mp2exp+105,256);
     Rand_x(p1,n,c-5,mp2exp+5,256);
     Rand_x(q1,n,c-8,mp2exp+89,256);
     mpz_sub_ui(temp1,p,1);
     mpz_sub_ui(temp2,q,1);
     mpz_mul(phi,temp1, temp2);
     mpz_gcdext(n,temp1, temp2, e, phi);
     mpz_mod(temp2,temp1,phi);
     mpz_set(d,temp2);
     mpz_mul(n,p,q);
     mpz_sub_ui(temp1,p1,1);
     mpz_sub_ui(temp2,q1,1);
     mpz_mul(phi1,temp1, temp2);
     mpz_gcdext(n1,temp1, temp2, e1, phi1);
     mpz_mod(temp2,temp1,phi1);
     mpz_set(d1,temp2);
     mpz_mul(n1,p1,q1);
     mpz_t C,M;
     mpz_init(C); mpz_init(M);
     if (mpz_cmp(n,n1)>0)
     {
      fputs("A:\n public key\n e= ", rsa);
      mpz_out_str(rsa, 16,e);
      fputs("\n n= ", rsa);
      mpz_out_str(rsa,16,n);
      fputs("\n private key\n d= ", rsa);
      mpz_out_str(rsa,16,d);
      fputs("\n phi(n)= ", rsa);
      mpz_out_str(rsa,16,phi);
      fputs("\np= ", rsa);
      mpz_out_str(rsa,16,p);
      fputs("\nq= ", rsa);
      mpz_out_str(rsa,16,q);
      fputs("\n decrypt/encrypt for A", rsa);
      Rand_x_modp(M,n,temp1,c,mp2exp);
      fputs("\n M= ", rsa);
      mpz_out_str(rsa,16,M);
      encrypt(C,M,e,n);
      fputs("\n C= ", rsa);
      mpz_out_str(rsa,16,C);
      decrypt(M,C,d,n);
      fputs("\n M= ", rsa);
      mpz_out_str(rsa,16,M);  
      fputs("\n B:\n public key\n e1= ", rsa);
      mpz_out_str(rsa, 16,e1);
      fputs("\n n1= ", rsa);
      mpz_out_str(rsa,16,n1);
     fputs("\n private key\n d1= ", rsa);
      mpz_out_str(rsa,16,d1);
      fputs("\n phi(n1)= ", rsa);
      mpz_out_str(rsa,16,phi1);  
      fputs("\np1= ", rsa);
      mpz_out_str(rsa,16,p1);
      fputs("\nq1= ", rsa);
      mpz_out_str(rsa,16,q1); 
      fputs("\nSignature. B is a sender", rsa);
       Rand_x_modp(M,n,temp1,c,mp2exp);
      fputs("\n M= ", rsa);
      mpz_out_str(rsa,16,M);
      signat(C,M,d1,n1);
      fputs("\n S= ", rsa);
      mpz_out_str(rsa,16,C);
      encrypt(M,C,e1,n1);
      fputs("\n M= ", rsa);  
      mpz_out_str(rsa,16,M); 
      }
     else
     {
     fputs("A:\n public key\n e= ", rsa);
      mpz_out_str(rsa, 16,e1);
      fputs("\n n= ", rsa);
      mpz_out_str(rsa,16,n1);
     fputs("\n private key\n d= ", rsa);
      mpz_out_str(rsa,16,d1);
      fputs("\n phi(n)= ", rsa);
      mpz_out_str(rsa,16,phi1); 
      fputs("\np= ", rsa);
      mpz_out_str(rsa,16,p1);
      fputs("\nq= ", rsa);
      mpz_out_str(rsa,16,q1);
      fputs("\n decrypt/encrypt for A", rsa);
      Rand_x_modp(M,n,temp1,c,mp2exp);
      fputs("\n M= ", rsa);
      mpz_out_str(rsa,16,M);
      encrypt(C,M,e,n);
      fputs("\n C= ", rsa);
      mpz_out_str(rsa,16,C);
      decrypt(M,C,d,n);
      fputs("\n M= ", rsa);
      mpz_out_str(rsa,16,M); 
     fputs("\nB:\n public key\n e1= ", rsa);
      mpz_out_str(rsa, 16,e);
      fputs("\n n1= ", rsa);
      mpz_out_str(rsa,16,n);
     fputs("\n private key\n d1= ", rsa);
      mpz_out_str(rsa,16,d);
      fputs("\n phi(n1)= ", rsa);
      mpz_out_str(rsa,16,phi);
      fputs("\np1= ", rsa);
      mpz_out_str(rsa,16,p);
      fputs("\nq1= ", rsa);
      mpz_out_str(rsa,16,q);
      fputs("\nSignature. B is a sender", rsa);
      Rand_x_modp(M,n,temp1,c,mp2exp);
      fputs("\n M= ", rsa);
      mpz_out_str(rsa,16,M);
      signat(C,M,d,n);
      fputs("\n S= ", rsa);
      mpz_out_str(rsa,16,C);
      encrypt(M,C,e,n);
      fputs("\n M= ", rsa);
       mpz_out_str(rsa,16,M);
     }
}
 
//void encrypt(mpz_t C, mpz_t M, mpz_t e, mpz_t n)
 //{    Horner(C,M,e,n);  }

void encrypt1(mpz_t C, mpz_t M, mpz_t e, mpz_t n)
{    Horner(C,M,e,n);  }
 
void decrypt(mpz_t M, mpz_t C, mpz_t d, mpz_t n)
{     Horner(M,C,d,n); }

void signat(mpz_t S, mpz_t M, mpz_t d, mpz_t n)
{     Horner(S,M,d,n); }
