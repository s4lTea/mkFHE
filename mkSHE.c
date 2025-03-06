#include "flint/fmpz_mod_poly.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_vec.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/* #include "fmpz.h" */

// LD_LIBRARY_PATH=/usr/local/lib
// export LD_LIBRARY_PATH

double PI = 3.1415926;
int bit_num;
int sigma;       // 离散高斯分布标准差
int B;           //[-B, B] 近似离散高斯分布
int d;           // 同态加密f = x^d + 1
int lamda;       // 系统安全参数
int Max = 10000; // 以概率ro 输出x的方法：随机生成 【0， 10000】
                 // 之间的数z，如果z小于（小于等于）ro*10000则输出x，否则不输出
int delta;       // 系统安全参数
int l;
fmpz_t q;
fmpz_t neg_Delta;
fmpz_t t;
fmpz_t Delta_2;
fmpz_t Delta;
fmpz_t T;
fmpz_t T2;
fmpz_t T3;
/* fmpz_t* T; */

fmpz_poly_t f;
fmpz_mod_ctx_t qmod;

// 可达到128bit安全
void System_Param() {
  sigma = 3.17;
  B = 10 * sigma;
  /* d = 8192; */
  d = 2048;
  bit_num = 128;
  /* char *str1 = "100000000000000000000000000000000"; */
  /* char *str1 = "7f001"; */
  /* char *str1 = "8000000000000000000000000000001d"; */
  char *str1 = "2a5d00b";
  /* char *str2 = "8000"; */
  char *str2 = "2";
  fmpz_set_str(
      q, str1,
      16); // q = 34028236692093846346337460743176821145665536 = pow(2, 128)
  fmpz_mod_ctx_init(qmod, q);
  fmpz_set_str(t, str2, 16); // t = pow(2, 15)=32768
  fmpz_cdiv_q(Delta, q, t);
  fmpz_neg(neg_Delta, Delta);
  fmpz_cdiv_q_ui(Delta_2, Delta, 2);
  /* fmpz_sqrt(T2, q); // pow(2, 64) = 18446744073709551616 */
  /* fmpz_sqrt(T, T2); // pow(2, 64) = 18446744073709551616 */
  fmpz_root(T, q, 4);
  fmpz_pow_ui(T2, T, 2);
  fmpz_pow_ui(T3, T, 3);
  l = fmpz_flog(q, T);
  /* printf("l is %d\n",  l); */

  fmpz_poly_init(f);
  fmpz_poly_set_coeff_si(f, 0, 1);
  /* for (int i = d - 1; i >= 1; i--) { */
  /*   fmpz_poly_set_coeff_si(f, i, 0); */
  /* } */
  fmpz_poly_set_coeff_si(f, d, 1);
}

// D_{Z, \sigma}
long SampleZ(int seed) {
  long xx;
  long z;
  double ro;
  long zz;
  int b;

  srand(time(NULL) + seed);

  do {
    xx = rand() % B;
    z = rand() % Max;
    ro = exp(-(PI * xx * xx) / (sigma * sigma));
  } while (!(z < ro * Max));

  b = rand() % 2;
  if (b == 0)
    b = -1;
  xx = xx * b;

  return xx;
}

// 生成（近似）离散高斯分布
void SampleD(fmpz_poly_t v, int seed) {
  int i;
  long x;

  for (i = d - 1; i >= 0; i--) {
    x = SampleZ(i + seed);
    fmpz_poly_set_coeff_si(v, i, x);
  }
}

// s和u都是R2上的多项式
void Gen_R2(fmpz_poly_t v, int seed) {
  int i;
  srand(time(NULL) + seed);
  int b;
  for (i = d - 1; i >= 0; i--) {
    b = rand() % 2;
    fmpz_poly_set_coeff_ui(v, i, b);
  }
}

void Gen_SR2(fmpz_poly_t v, int seed) {
  int i;
  srand(time(NULL) + seed);
  int b;
  for (i = d - 1; i >= 0; i--) {
    b = rand() % 3;
    if (b == 2)
      b = -1;
    fmpz_poly_set_coeff_ui(v, i, b);
  }
}

void Gen_Rq_div_2(fmpz_t t, int seed) {
  int i;
  char *str = (char *)malloc(sizeof(char) * (bit_num - 1));
  int b;
  srand(time(NULL) + seed);

  for (i = bit_num - 2; i >= 0; i--)
    str[i] = rand() % 2 + 48;
  fmpz_set_str(t, str, 2);
  b = rand() % 2;

  if (b == 0)
    b = -1;
  fmpz_mul_si(t, t, b);
}

void Gen_Rq_div_pm2(fmpz_t t, int seed) {
  int i;
  char *str = (char *)malloc(sizeof(char) * (bit_num - 1));
  int b;
  srand(time(NULL) + seed);

  for (i = bit_num - 2; i >= 0; i--)
    str[i] = rand() % 2 + 48;
  fmpz_set_str(t, str, 2);
  b = rand() % 3;

  if (b == 3)
    b = -1;
  fmpz_mul_si(t, t, b);
}

/*
void Gen_Rq(fmpz_poly_t v)
{
  int i;
  fmpz_t x;
  fmpz_init(x);

  flint_rand_t state;
  flint_randinit(state);

  for(i = d-1; i >= 0; i--)
     {
        fmpz_randtest_mod_signed(x, state, q);
        fmpz_poly_set_coeff_fmpz(v, i, x);
     }

  fmpz_clear(x);
  flint_randclear(state);

}
*/

// 这个的问题是好像没有小一点的数。。计算上应该效率要差一些，但是这个正常的吗？不是很清楚
void Gen_Rq(fmpz_poly_t v) {
  int i;
  fmpz_t t;
  fmpz_init(t);
  for (i = d - 1; i >= 0; i--) {
    Gen_Rq_div_2(t, i);
    fmpz_poly_set_coeff_fmpz(v, i, t);
  }

  fmpz_clear(t);
}

void secret_key_Gen(fmpz_poly_t sk) { Gen_R2(sk, 1); }

void public_key_Gen(fmpz_poly_t sk, fmpz_poly_t pk_p0, fmpz_poly_t pk_p1) {
  // pk = ( [-(as+e)]_q, a ), a \in R_q, s \in R_2,  e \in D_{Z, sigma}^d
  fmpz_poly_t e;
  fmpz_poly_t temp;
  fmpz_poly_init(e);
  fmpz_poly_init(temp);

  SampleD(e, 1);
  Gen_Rq(pk_p1);

  fmpz_poly_mul(temp, pk_p1, sk);
  fmpz_poly_add(temp, temp, e);
  fmpz_poly_neg(temp, temp);
  fmpz_poly_scalar_smod_fmpz(temp, temp, q);
  fmpz_poly_set(pk_p0, temp);
  fmpz_poly_clear(e);
  fmpz_poly_clear(temp);
}

void SH_Keygen(fmpz_poly_t sk, fmpz_poly_t pk_p0, fmpz_poly_t pk_p1,
               fmpz_poly_t rlk_r0, fmpz_poly_t rlk_r1) {
  secret_key_Gen(sk);               // 生成私钥
  public_key_Gen(sk, pk_p0, pk_p1); // 生成公钥
  ReKey_Gen(sk, rlk_r0, rlk_r1);    // 生成Relinear key
}

void SH_Keygen_ntru(fmpz_poly_t sk, fmpz_poly_t skinvs, fmpz_poly_t pk_p0,
                    fmpz_poly_t pk_p1, fmpz_poly_t rlk_r0) {
  // We need generate a sk such that skinvs exists
  //
  secret_key_Gen(sk); // 生成私钥
                      // Use fmpz_mod_poly_invmod function to
                      //
                      // calculate skinvs. We need keep sk and f
                      /* nmod_poly_t temp; */
                      /**/
                      /* nmod_poly_init_mod(temp, q); */

  fmpz_poly_t temp;
  fmpz_poly_init(temp);

  // find invertible denominator g
  /* this->rand_poly(this->den_); */
  /* nmod_poly_gcd(temp, this->den_, this->modulus); */
  fmpz_poly_gcd(temp, sk, f);

  /* while (!nmod_poly_is_one(temp)) { */
  while (!fmpz_poly_is_one(temp)) {
    /* this->rand_poly(this->den_); */
    secret_key_Gen(sk); // 生成私钥
    fmpz_poly_gcd(temp, sk, f);

    /* nmod_poly_gcd(temp, this->den_, this->modulus); */
  }

  // std::cout << "# g = ";
  // nmod_poly_print_pretty(this->den_, var);
  // std::cout << '\n';
  /* nmod_poly_t g_inv; */

  /* fmpz_mod_poly_t skinvs; */
  // Below code is to deal with Rq where q is not a prime number over modulus
  // x^n + 1 The similar sage code is available at
  // https://github.com/kpatsakis/NTRU_Sage/blob/5eed52be90e1c9f5e0ec44831a2554558865192e/ntru.sage#L122
  // We need to calculate the inverse of poly g over R_{q^r}[x]/x^n + 1
  // We can first calculate the inverse poly f for poly g over R_{q}[x] / x^n +
  // 1 And for i = 2 to ceil (log_2 r) Do below code: b = b * ( 2 - poly * b) b
  // = b % (x^n + 1) b = b mod q // Let the coefficients in the [0, q - 1]
  /* fmpz_mod_ctx_t temp_mod; */
  /*   fmpz_t temp_modulus = {2}; */
  /* fmpz_mod_ctx_init(temp_mod, temp_modulus); */
  /* fmpz_mod_poly_invmod(skinvs, sk, f, temp_mod); */
  /* nmod_poly_init_mod(g_inv, q_nmod); */
  /* nmod_poly_invmod(g_inv, this->den_, this->modulus); */

  // Below code calculate the inversion of sk over R_{q}[x] / (x^d + 1) where d
  // is equal to 1024 in our parameters. And q is 0x7f001 is a prime, which I
  // have checked in sage.
  fmpz_mod_ctx_t temp_mod;
  /* fmpz_t temp_modulus = {0x8000000000000000000000000000001d}; */

  fmpz_mod_ctx_init(temp_mod, q);
  fmpz_mod_poly_invmod(skinvs, sk, f, temp_mod);
  fmpz_poly_scalar_mod_fmpz(skinvs, skinvs, q);

  // nmod_poly_gcdinv(res, g_inv, keygen.g, keygen.modulus);

  // check g_inv
  /* nmod_poly_mul(temp, this->den_, g_inv); */
  /* nmod_poly_rem(temp, temp, this->modulus); */
  // check skinvs
  fmpz_poly_mul(temp, sk, skinvs);
  fmpz_poly_rem(temp, temp, f);

  fmpz_poly_scalar_mod_fmpz(temp, temp, q);
  /* if (!nmod_poly_is_one(temp)) { */
  /*   debug("Something went wrong in invmod!\n"); */
  /*   // nmod_poly_print_pretty(temp, var); */
  /*   // std::cout << '\n'; */
  /*   return 0; */
  /* } */

  // Check
  if (!fmpz_poly_is_one(temp)) {
    printf("There is something wrong when calculating skinvs\n");
  }

  /* while (1) { */
  /*   secret_key_Gen(sk); // 生成私钥 */
  /* } */

  public_key_Gen(sk, pk_p0, pk_p1);   // 生成公钥
  ReKey_Gen_ntru(skinvs, sk, rlk_r0); // 生成Relinear key with NTRU
}

void Encode(int m, fmpz_poly_t pm) {
  int sgn;
  if (m < 0)
    sgn = -1;
  else
    sgn = 1;
  int num;
  m = abs(m);
  double bitnum = log2(m);
  int i;
  int bit;
  if (bitnum == (int)bitnum)
    num = (int)bitnum + 1;
  else
    num = (int)bitnum + 1;

  for (i = 0; i < num; i++) {
    bit = m % 2;
    fmpz_poly_set_coeff_ui(pm, i, bit);
    m = m / 2;
  }
  fmpz_poly_scalar_mul_si(pm, pm, sgn);
}

void SH_Encrypt(fmpz_poly_t m, fmpz_poly_t pk_p0, fmpz_poly_t pk_p1,
                fmpz_poly_t c0, fmpz_poly_t c1) {
  fmpz_poly_t u;
  fmpz_poly_t e1;
  fmpz_poly_t e2;
  fmpz_poly_t temp;
  fmpz_poly_t temp1;
  fmpz_poly_init(u);
  fmpz_poly_init(e1);
  fmpz_poly_init(e2);
  fmpz_poly_init(temp);
  fmpz_poly_init(temp1);

  Gen_R2(u, 2);

  SampleD(e1, 2);
  SampleD(e2, 3);

  fmpz_poly_mul(temp, pk_p0, u);

  fmpz_poly_add(temp, temp, e1);

  fmpz_poly_scalar_mul_fmpz(temp1, m, Delta);

  fmpz_poly_add(temp, temp, temp1);

  fmpz_poly_scalar_smod_fmpz(temp, temp, q);

  fmpz_poly_set(c0, temp);

  fmpz_poly_mul(temp, pk_p1, u);

  fmpz_poly_add(temp, temp, e2);

  fmpz_poly_scalar_smod_fmpz(c1, temp, q);

  fmpz_poly_clear(u);
  fmpz_poly_clear(e1);
  fmpz_poly_clear(e2);
  fmpz_poly_clear(temp);
  fmpz_poly_clear(temp1);
}

void fmpz_poly_nearest_fmpz(fmpz_poly_t v) {
  long deg = fmpz_poly_degree(v);
  int i;
  fmpz_t temp;
  fmpz_t temp1;
  fmpz_init(temp);
  fmpz_init(temp1);

  for (i = 0; i <= deg; i++) {
    fmpz_poly_get_coeff_fmpz(temp, v, i);

    fmpz_abs(temp1, temp);

    fmpz_fdiv_r(temp1, temp1, Delta);

    if (fmpz_equal(temp, Delta) == 1) {
      fmpz_set_ui(temp, 1);
    } else if (fmpz_equal(temp, neg_Delta) == 1) {
      fmpz_set_si(temp, -1);
    } else if (fmpz_is_zero(temp) == 1) {
      fmpz_set_ui(temp, 0);
    }

    else if (fmpz_sgn(temp) == -1 && fmpz_cmp(temp1, Delta_2) > 0) {
      fmpz_fdiv_q(temp, temp, Delta);
    } else if (fmpz_sgn(temp) == -1 && fmpz_cmp(temp1, Delta_2) <= 0) {
      fmpz_cdiv_q(temp, temp, Delta);
    } else if (fmpz_sgn(temp) == 1 && fmpz_cmp(temp1, Delta_2) >= 0) {
      fmpz_cdiv_q(temp, temp, Delta);
    } else {
      fmpz_fdiv_q(temp, temp, Delta);
    }

    fmpz_poly_set_coeff_fmpz(v, i, temp);
  }

  fmpz_clear(temp);
  fmpz_clear(temp1);
}

// c0 = (c_00, c_01), c1 = (c_10, c_11), c = c0+c1

void SH_Add(fmpz_poly_t c_00, fmpz_poly_t c_01, fmpz_poly_t c_10,
            fmpz_poly_t c_11, fmpz_poly_t c0, fmpz_poly_t c1) {
  fmpz_poly_t temp;
  fmpz_poly_init(temp);
  fmpz_poly_add(temp, c_00, c_10);
  fmpz_poly_scalar_smod_fmpz(c0, temp, q);

  fmpz_poly_add(temp, c_01, c_11);
  fmpz_poly_scalar_smod_fmpz(c1, temp, q);
  fmpz_poly_clear(temp);
}

void SH_Mul(fmpz_poly_t c_00, fmpz_poly_t c_01, fmpz_poly_t c_10,
            fmpz_poly_t c_11, fmpz_poly_t c0, fmpz_poly_t c1, fmpz_poly_t c2) {
  fmpz_poly_t temp1;
  fmpz_poly_t temp2;

  fmpz_poly_init(temp1);
  fmpz_poly_init(temp2);

  fmpz_poly_mul(c0, c_00, c_10);
  fmpz_poly_nearest_fmpz(c0);

  fmpz_poly_mul(temp1, c_00, c_11);
  fmpz_poly_mul(temp2, c_01, c_10);
  fmpz_poly_add(c1, temp1, temp2);
  fmpz_poly_nearest_fmpz(c1);

  fmpz_poly_mul(c2, c_01, c_11);
  fmpz_poly_nearest_fmpz(c2);

  fmpz_poly_clear(temp1);
  fmpz_poly_clear(temp2);
}

void SH_DecMul(fmpz_poly_t sk, fmpz_poly_t c0, fmpz_poly_t c1, fmpz_poly_t c2) {
  fmpz_poly_t temp;
  fmpz_poly_t temp2;
  fmpz_poly_init(temp);
  fmpz_poly_init(temp2);

  fmpz_poly_mul(temp, c1, sk);

  fmpz_poly_mul(temp2, c2, sk);
  fmpz_poly_mul(temp2, temp2, sk);

  fmpz_poly_add(temp, c0, temp);
  fmpz_poly_add(temp, temp, temp2);

  fmpz_poly_scalar_smod_fmpz(temp, temp, q);

  fmpz_poly_nearest_fmpz(temp);

  fmpz_poly_scalar_smod_fmpz(temp, temp, t);
  printf("\ntemp*t/q is \n");
  fmpz_poly_print(temp);

  fmpz_poly_clear(temp);
  fmpz_poly_clear(temp2);
}

void ReKey_Gen(fmpz_poly_t sk, fmpz_poly_t *rlk_r0, fmpz_poly_t *rlk_r1) {
  fmpz_poly_t e;
  int i;
  fmpz_poly_t temp1;
  fmpz_t temp2;
  fmpz_poly_t temp3;
  fmpz_poly_init(e);
  fmpz_poly_init(temp1);
  fmpz_init(temp2);
  fmpz_poly_init(temp3);

  for (i = 0; i < l; i++) {
    Gen_Rq(rlk_r1[i]);
    SampleD(e, i + 9);
    fmpz_poly_mul(temp1, rlk_r1[i], sk);
    fmpz_poly_add(temp1, temp1, e);
    fmpz_poly_neg(temp1, temp1);
    fmpz_pow_ui(temp2, T, i);
    fmpz_poly_mul(temp3, sk, sk);
    fmpz_poly_scalar_mul_fmpz(temp3, temp3, temp2);
    fmpz_poly_add(temp1, temp1, temp3);
    fmpz_poly_scalar_smod_fmpz(rlk_r0[i], temp1, q);
  }

  fmpz_poly_clear(temp1);
  fmpz_clear(temp2);
  fmpz_poly_clear(temp3);
  fmpz_poly_clear(e);
}

void ReKey_Gen_ntru(fmpz_poly_t skinvs, fmpz_poly_t sk, fmpz_poly_t *rlk_r0) {
  /* fmpz_poly_t e; */
  int i;
  fmpz_poly_t temp1;
  fmpz_t temp2;
  fmpz_poly_t temp3;
  /* fmpz_poly_init(e); */
  fmpz_poly_init(temp1);
  fmpz_init(temp2);
  fmpz_poly_init(temp3);

  // hi = gi / s
  // evki = hi + w^i s
  for (i = 0; i < l; i++) {
    Gen_Rq_div_pm2(rlk_r0[i], i + 9);
    /* Gen_Rq(rlk_r1[i]); */
    /* SampleD(e, i + 9); */
    fmpz_poly_mul(temp1, rlk_r0[i], skinvs);
    fmpz_poly_scalar_mod_fmpz(temp1, temp1, q);
    /* fmpz_poly_add(temp1, temp1, e); */
    /* fmpz_poly_neg(temp1, temp1); */
    fmpz_pow_ui(temp2, T, i);
    /* fmpz_poly_mul(temp3, sk, sk); */
    fmpz_poly_scalar_mul_fmpz(temp3, sk, temp2);
    fmpz_poly_add(temp1, temp1, temp3);
    fmpz_poly_scalar_mod_fmpz(rlk_r0[i], temp1, q);
  }

  fmpz_poly_clear(temp1);
  fmpz_clear(temp2);
  fmpz_poly_clear(temp3);
  /* fmpz_poly_clear(e); */
}

// T=sqrt(q), 所以到了i = 2和以后更大的i就没有必要再算了，因为mod q都是0
void Relinear(fmpz_poly_t *rlk_r0, fmpz_poly_t *rlk_r1, fmpz_poly_t c0,
              fmpz_poly_t c1, fmpz_poly_t c2, fmpz_poly_t _c0,
              fmpz_poly_t _c1) {
  int i;
  int deg = fmpz_poly_degree(c2);

  fmpz_poly_t *c2_i;

  c2_i = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * l);

  fmpz_t temp;
  fmpz_t temp1;
  fmpz_t temp2;
  fmpz_t temp3;
  fmpz_poly_t temp4;
  fmpz_poly_t temp5;
  fmpz_poly_t temp6;
  fmpz_poly_t temp7;

  for (i = 0; i < l; i++) {
    fmpz_poly_init(c2_i[i]);
  }

  fmpz_init(temp);
  fmpz_init(temp1);
  fmpz_init(temp2);
  fmpz_init(temp3);
  fmpz_poly_init(temp4);
  fmpz_poly_init(temp5);
  fmpz_poly_init(temp6);
  fmpz_poly_init(temp7);

  for (i = 0; i <= deg; i++) {
    fmpz_poly_get_coeff_fmpz(temp, c2, i);

    fmpz_abs(temp1, temp);
    if (fmpz_cmp(temp1, T) < 0)
      fmpz_poly_set_coeff_fmpz(c2_i[0], i, temp);
    else {
      // temp * T^2 + temp1
      fmpz_fdiv_qr(temp3, temp2, temp, T3);
      // temp = temp2 * T + temp3
      fmpz_fdiv_qr(temp2, temp1, temp2, T2);
      // temp1 = temp2 * T + temp3
      fmpz_fdiv_qr(temp1, temp, temp1, T);

      fmpz_poly_set_coeff_fmpz(c2_i[0], i, temp);
      fmpz_poly_set_coeff_fmpz(c2_i[1], i, temp1);
      fmpz_poly_set_coeff_fmpz(c2_i[2], i, temp2);
      fmpz_poly_set_coeff_fmpz(c2_i[3], i, temp3);
    }
  }

  fmpz_poly_mul(temp4, rlk_r0[0], c2_i[0]);
  fmpz_poly_mul(temp5, rlk_r0[1], c2_i[1]);
  fmpz_poly_mul(temp6, rlk_r0[2], c2_i[2]);
  fmpz_poly_mul(temp7, rlk_r0[3], c2_i[3]);

  fmpz_poly_add(temp4, temp4, temp5);
  fmpz_poly_add(temp6, temp6, temp7);
  fmpz_poly_add(temp4, temp4, temp6);
  fmpz_poly_add(temp4, temp4, c0);
  fmpz_poly_scalar_mod_fmpz(_c0, temp4, q);

  fmpz_poly_mul(temp4, rlk_r1[0], c2_i[0]);
  fmpz_poly_mul(temp5, rlk_r1[1], c2_i[1]);
  fmpz_poly_mul(temp6, rlk_r1[2], c2_i[2]);
  fmpz_poly_mul(temp7, rlk_r1[3], c2_i[3]);

  fmpz_poly_add(temp4, temp4, temp5);
  fmpz_poly_add(temp6, temp6, temp7);
  fmpz_poly_add(temp4, temp4, temp6);
  fmpz_poly_add(temp4, temp4, c1);
  fmpz_poly_scalar_mod_fmpz(_c1, temp4, q);

  for (i = 0; i < l; i++)
    fmpz_poly_clear(c2_i[i]);

  fmpz_clear(temp);
  fmpz_clear(temp1);
  fmpz_clear(temp3);
  fmpz_clear(temp4);
  fmpz_poly_clear(temp4);
  fmpz_poly_clear(temp5);
  fmpz_poly_clear(temp6);
  fmpz_poly_clear(temp7);
}

void Relinear_ntru(fmpz_poly_t *rlk_r0, fmpz_poly_t c0, fmpz_poly_t c1,
                   fmpz_poly_t c2, fmpz_poly_t _c0, fmpz_poly_t _c1) {
  // We need keep c0 and change c1
  int i;
  int deg = fmpz_poly_degree(c2);

  fmpz_poly_t *c2_i;
  c2_i = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * l);

  fmpz_t temp;
  fmpz_t temp1;
  fmpz_t temp2;
  fmpz_t temp3;
  fmpz_poly_t temp4;
  fmpz_poly_t temp5;
  fmpz_poly_t temp6;
  fmpz_poly_t temp7;

  for (i = 0; i < l; i++)
    fmpz_poly_init(c2_i[i]);

  fmpz_init(temp);
  fmpz_init(temp1);
  fmpz_init(temp2);
  fmpz_init(temp3);
  fmpz_poly_init(temp4);
  fmpz_poly_init(temp5);
  fmpz_poly_init(temp6);
  fmpz_poly_init(temp7);

  for (i = 0; i <= deg; i++) {
    fmpz_poly_get_coeff_fmpz(temp, c2, i);

    fmpz_abs(temp1, temp);
    if (fmpz_cmp(temp1, T) < 0)
      fmpz_poly_set_coeff_fmpz(c2_i[0], i, temp);
    else {
      /* fmpz_fdiv_qr(temp, temp1, temp, T); */
      /* fmpz_poly_set_coeff_fmpz(c2_i[0], i, temp1); */
      /* fmpz_poly_set_coeff_fmpz(c2_i[1], i, temp); */
      // temp1 * T2 + temp = temp
      /* fmpz_fdiv_qr(temp, temp1, temp, T2); */
      /* // temp = temp3 * T + temp2 */
      /* fmpz_fdiv_qr(temp3, temp2, temp, T); */
      /* // temp1 = temp1 * T + temp */
      /* fmpz_fdiv_qr(temp1, temp, temp1, T); */

      fmpz_fdiv_qr(temp3, temp2, temp, T3);
      fmpz_fdiv_qr(temp2, temp1, temp2, T2);
      fmpz_fdiv_qr(temp1, temp, temp1, T);
      fmpz_poly_set_coeff_fmpz(c2_i[0], i, temp);
      fmpz_poly_set_coeff_fmpz(c2_i[1], i, temp1);
      fmpz_poly_set_coeff_fmpz(c2_i[2], i, temp2);
      fmpz_poly_set_coeff_fmpz(c2_i[3], i, temp3);
    }
  }
  // We keep c0
  /* fmpz_poly_mul(temp2, rlk_r0[0], c2_i[0]); */
  /* fmpz_poly_mul(temp3, rlk_r0[1], c2_i[1]); */
  /**/
  /* fmpz_poly_add(temp2, temp2, temp3); */
  /* fmpz_poly_add(temp2, temp2, c0); */
  /* fmpz_poly_scalar_smod_fmpz(_c0, temp2, q); */
  fmpz_poly_set(_c0, c0);

  fmpz_poly_mul(temp4, rlk_r0[0], c2_i[0]);
  fmpz_poly_mul(temp5, rlk_r0[1], c2_i[1]);
  fmpz_poly_mul(temp6, rlk_r0[2], c2_i[2]);
  fmpz_poly_mul(temp7, rlk_r0[3], c2_i[3]);

  fmpz_poly_add(temp4, temp4, temp5);
  fmpz_poly_add(temp6, temp6, temp7);
  fmpz_poly_add(temp4, temp4, temp6);
  fmpz_poly_add(temp4, temp4, c1);
  fmpz_poly_scalar_smod_fmpz(_c1, temp4, q);

  for (i = 0; i < l; i++)
    fmpz_poly_clear(c2_i[i]);

  fmpz_clear(temp);
  fmpz_clear(temp1);
  fmpz_clear(temp2);
  fmpz_clear(temp3);
  fmpz_poly_clear(temp4);
  fmpz_poly_clear(temp5);
  fmpz_poly_clear(temp6);
  fmpz_poly_clear(temp7);
}

void SH_Decrypt(fmpz_poly_t sk, fmpz_poly_t c0, fmpz_poly_t c1) {
  fmpz_poly_t temp;
  fmpz_poly_init(temp);
  fmpz_poly_mul(temp, c1, sk);
  fmpz_poly_add(temp, temp, c0);
  fmpz_poly_scalar_smod_fmpz(temp, temp, q);
  fmpz_poly_nearest_fmpz(temp);
  fmpz_poly_scalar_smod_fmpz(temp, temp, t);

  printf("\ntemp*t/q is \n");
  fmpz_poly_print(temp);

  fmpz_poly_clear(temp);
}

struct timespec diff(struct timespec start, struct timespec end) {
  struct timespec temp;
  if ((end.tv_nsec - start.tv_nsec) < 0) {
    temp.tv_sec = end.tv_sec - start.tv_sec - 1;
    temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec - start.tv_sec;
    temp.tv_nsec = end.tv_nsec - start.tv_nsec;
  }
  return temp;
}

void test_bfv() {

  fmpz_poly_t pk_p0;
  fmpz_poly_t pk_p1;
  fmpz_poly_t sk;
  fmpz_poly_t *rlk_r0;
  fmpz_poly_t *rlk_r1;
  int m;
  int i;
  fmpz_poly_t m0;
  fmpz_poly_t m1;
  fmpz_poly_t c_00;
  fmpz_poly_t c_01;
  fmpz_poly_t c_10;
  fmpz_poly_t c_11;
  fmpz_poly_t c0;
  fmpz_poly_t c1;
  fmpz_poly_t c2;
  fmpz_poly_t _c0;
  fmpz_poly_t _c1;
  struct timespec tpstart_full;
  struct timespec tpend_full;
  struct timespec tp_full;
  struct timespec tpstart_mul;
  struct timespec tpend_mul;
  struct timespec tp_mul;
  struct timespec tpstart_rel;
  struct timespec tpend_rel;
  struct timespec tp_rel;
  long timedif;

  fmpz_init(q);
  fmpz_init(neg_Delta);
  fmpz_init(Delta_2);
  fmpz_init(t);
  fmpz_init(Delta);
  fmpz_init(T);
  fmpz_init(T2);

  fmpz_poly_init(pk_p0);
  fmpz_poly_init(pk_p1);
  fmpz_poly_init(sk);
  fmpz_poly_init(m0);
  fmpz_poly_init(m1);
  fmpz_poly_init(c_00);
  fmpz_poly_init(c_01);
  fmpz_poly_init(c_10);
  fmpz_poly_init(c_11);
  fmpz_poly_init(c0);
  fmpz_poly_init(c1);
  fmpz_poly_init(c2);
  fmpz_poly_init(_c0);
  fmpz_poly_init(_c1);

  System_Param();

  rlk_r0 = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * (l));
  rlk_r1 = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * (l));

  for (i = 0; i < l; i++) {
    fmpz_poly_init(rlk_r0[i]);
    fmpz_poly_init(rlk_r1[i]);
  }

  // test begin
  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  SH_Keygen(sk, pk_p0, pk_p1, rlk_r0, rlk_r1);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nKEYGEN%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  scanf("%d", &m);
  Encode(m, m0);
  scanf("%d", &m);
  Encode(m, m1);

  clock_gettime(CLOCK_MONOTONIC, &tpstart_full);
  SH_Encrypt(m0, pk_p0, pk_p1, c_00, c_01);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nENC%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  SH_Encrypt(m1, pk_p0, pk_p1, c_10, c_11);

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  /* SH_Decrypt(sk, c_00, c_01); */
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nDEC%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  SH_Decrypt(sk, c_10, c_11);

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  /* SH_Add(c_00, c_01, c_10, c_11, c0, c1); */
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nADD%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  /* SH_Decrypt(sk, c0, c1); */

  clock_gettime(CLOCK_MONOTONIC, &tpstart_mul);
  SH_Mul(c_00, c_01, c_10, c_11, c0, c1, c2);

  // SH_DecMul(sk, c0, c1, c2);

  clock_gettime(CLOCK_MONOTONIC, &tpstart_rel);
  Relinear(rlk_r0, rlk_r1, c0, c1, c2, _c0, _c1);
  clock_gettime(CLOCK_MONOTONIC, &tpend_mul);
  clock_gettime(CLOCK_MONOTONIC, &tpend_rel);

  tp_mul = diff(tpstart_mul, tpend_mul);
  printf("\nMUL%ds:%ldns \n", tp_mul.tv_sec, tp_mul.tv_nsec);

  tp_rel = diff(tpstart_rel, tpend_rel);
  printf("\nREL%ds:%ldns \n", tp_rel.tv_sec, tp_rel.tv_nsec);

  SH_Decrypt(sk, _c0, _c1);
  clock_gettime(CLOCK_MONOTONIC, &tpend_full);
  tp_full = diff(tpstart_full, tpend_full);
  printf("\nfull%ds:%ldns \n", tp_full.tv_sec, tp_full.tv_nsec);

  fmpz_poly_clear(sk);
  fmpz_poly_clear(pk_p0);
  fmpz_poly_clear(pk_p1);
  fmpz_clear(q);
  fmpz_clear(neg_Delta);
  fmpz_clear(Delta_2);
  fmpz_clear(t);
  fmpz_clear(Delta);
  fmpz_clear(T);
  fmpz_poly_clear(c_00);
  fmpz_poly_clear(c_01);
  fmpz_poly_clear(c_10);
  fmpz_poly_clear(c_11);
  fmpz_poly_clear(c0);
  fmpz_poly_clear(c1);
  fmpz_poly_clear(c2);
  fmpz_poly_clear(_c0);
  fmpz_poly_clear(_c1);
  fmpz_poly_clear(m0);
  fmpz_poly_clear(m1);

  for (i = 0; i < l; i++) {
    fmpz_poly_clear(rlk_r0[i]);
    fmpz_poly_clear(rlk_r1[i]);
  }
}

void test_full_bfv(int m, int i) {

  fmpz_poly_t pk_p0;
  fmpz_poly_t pk_p1;
  fmpz_poly_t sk;
  fmpz_poly_t *rlk_r0;
  fmpz_poly_t *rlk_r1;
  fmpz_poly_t m0;
  fmpz_poly_t m1;
  fmpz_poly_t c_00;
  fmpz_poly_t c_01;
  fmpz_poly_t c_10;
  fmpz_poly_t c_11;
  fmpz_poly_t c0;
  fmpz_poly_t c1;
  fmpz_poly_t c2;
  fmpz_poly_t _c0;
  fmpz_poly_t _c1;
  long timedif;

  fmpz_init(q);
  fmpz_init(neg_Delta);
  fmpz_init(Delta_2);
  fmpz_init(t);
  fmpz_init(Delta);
  fmpz_init(T);
  fmpz_init(T2);

  fmpz_poly_init(pk_p0);
  fmpz_poly_init(pk_p1);
  fmpz_poly_init(sk);
  fmpz_poly_init(m0);
  fmpz_poly_init(m1);
  fmpz_poly_init(c_00);
  fmpz_poly_init(c_01);
  fmpz_poly_init(c_10);
  fmpz_poly_init(c_11);
  fmpz_poly_init(c0);
  fmpz_poly_init(c1);
  fmpz_poly_init(c2);
  fmpz_poly_init(_c0);
  fmpz_poly_init(_c1);

  System_Param();

  rlk_r0 = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * (l));
  rlk_r1 = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * (l));

  for (i = 0; i < l; i++) {
    fmpz_poly_init(rlk_r0[i]);
    fmpz_poly_init(rlk_r1[i]);
  }

  // test begin
  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  SH_Keygen(sk, pk_p0, pk_p1, rlk_r0, rlk_r1);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nKEYGEN%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  /* scanf("%d", &m); */
  Encode(m, m0);
  /* scanf("%d", &m); */
  Encode(m, m1);

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  SH_Encrypt(m0, pk_p0, pk_p1, c_00, c_01);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nENC%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  SH_Encrypt(m1, pk_p0, pk_p1, c_10, c_11);

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  /* SH_Decrypt(sk, c_00, c_01); */
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nDEC%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  /* SH_Decrypt(sk, c_10, c_11); */

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  SH_Add(c_00, c_01, c_10, c_11, c0, c1);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nADD%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  /* SH_Decrypt(sk, c0, c1); */

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  SH_Mul(c_00, c_01, c_10, c_11, c0, c1, c2);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nMUL%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  // SH_DecMul(sk, c0, c1, c2);

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  Relinear(rlk_r0, rlk_r1, c0, c1, c2, _c0, _c1);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nREL%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  SH_Decrypt(sk, _c0, _c1);

  fmpz_poly_clear(sk);
  fmpz_poly_clear(pk_p0);
  fmpz_poly_clear(pk_p1);
  fmpz_clear(q);
  fmpz_clear(neg_Delta);
  fmpz_clear(Delta_2);
  fmpz_clear(t);
  fmpz_clear(Delta);
  fmpz_clear(T);
  fmpz_clear(T2);
  fmpz_poly_clear(c_00);
  fmpz_poly_clear(c_01);
  fmpz_poly_clear(c_10);
  fmpz_poly_clear(c_11);
  fmpz_poly_clear(c0);
  fmpz_poly_clear(c1);
  fmpz_poly_clear(c2);
  fmpz_poly_clear(_c0);
  fmpz_poly_clear(_c1);
  fmpz_poly_clear(m0);
  fmpz_poly_clear(m1);

  for (i = 0; i < l; i++) {
    fmpz_poly_clear(rlk_r0[i]);
    fmpz_poly_clear(rlk_r1[i]);
  }
}

void test_ntru() {

  fmpz_poly_t pk_p0;
  fmpz_poly_t pk_p1;
  fmpz_poly_t sk;
  fmpz_poly_t *rlk_r0;
  fmpz_poly_t skinvs;
  int m;
  int i;
  fmpz_poly_t m0;
  fmpz_poly_t m1;
  fmpz_poly_t c_00;
  fmpz_poly_t c_01;
  fmpz_poly_t c_10;
  fmpz_poly_t c_11;
  fmpz_poly_t c0;
  fmpz_poly_t c1;
  fmpz_poly_t c2;
  fmpz_poly_t _c0;
  fmpz_poly_t _c1;
  struct timespec tpstart_full;
  struct timespec tpend_full;
  struct timespec tp_full;
  struct timespec tpstart_mul;
  struct timespec tpend_mul;
  struct timespec tp_mul;
  struct timespec tpstart_rel;
  struct timespec tpend_rel;
  struct timespec tp_rel;
  long timedif;

  fmpz_init(q);
  fmpz_init(neg_Delta);
  fmpz_init(Delta_2);
  fmpz_init(t);
  fmpz_init(Delta);
  fmpz_init(T);
  fmpz_init(T);

  fmpz_poly_init(pk_p0);
  fmpz_poly_init(pk_p1);
  fmpz_poly_init(skinvs);
  fmpz_poly_init(sk);
  fmpz_poly_init(m0);
  fmpz_poly_init(m1);
  fmpz_poly_init(c_00);
  fmpz_poly_init(c_01);
  fmpz_poly_init(c_10);
  fmpz_poly_init(c_11);
  fmpz_poly_init(c0);
  fmpz_poly_init(c1);
  fmpz_poly_init(c2);
  fmpz_poly_init(_c0);
  fmpz_poly_init(_c1);

  System_Param();

  rlk_r0 = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * (l));

  for (i = 0; i < l; i++) {
    fmpz_poly_init(rlk_r0[i]);
  }

  // test begin
  /* clock_gettime(CLOCK_MONOTONIC, &tpstart_full); */
  SH_Keygen_ntru(sk, skinvs, pk_p0, pk_p1, rlk_r0);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nKEYGEN%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  scanf("%d", &m);
  Encode(m, m0);
  scanf("%d", &m);
  Encode(m, m1);

  clock_gettime(CLOCK_MONOTONIC, &tpstart_full);
  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  SH_Encrypt(m0, pk_p0, pk_p1, c_00, c_01);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nENC%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  SH_Encrypt(m1, pk_p0, pk_p1, c_10, c_11);

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  /* SH_Decrypt(sk, c_00, c_01); */
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nDEC%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  /* SH_Decrypt(sk, c_10, c_11); */

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  /* SH_Add(c_00, c_01, c_10, c_11, c0, c1); */
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nADD%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  /* SH_Decrypt(sk, c0, c1); */

  clock_gettime(CLOCK_MONOTONIC, &tpstart_mul);
  SH_Mul(c_00, c_01, c_10, c_11, c0, c1, c2);

  // SH_DecMul(sk, c0, c1, c2);

  clock_gettime(CLOCK_MONOTONIC, &tpstart_rel);
  /* Relinear(rlk_r0, rlk_r1, c0, c1, c2, _c0, _c1); */
  Relinear_ntru(rlk_r0, c0, c1, c2, _c0, _c1);
  clock_gettime(CLOCK_MONOTONIC, &tpend_rel);
  clock_gettime(CLOCK_MONOTONIC, &tpend_mul);

  tp_mul = diff(tpstart_mul, tpend_mul);
  printf("\nMUL%ds:%ldns \n", tp_mul.tv_sec, tp_mul.tv_nsec);

  tp_rel = diff(tpstart_rel, tpend_rel);
  printf("\nREL%ds:%ldns \n", tp_rel.tv_sec, tp_rel.tv_nsec);

  SH_Decrypt(sk, _c0, _c1);
  clock_gettime(CLOCK_MONOTONIC, &tpend_full);
  tp_full = diff(tpstart_full, tpend_full);
  printf("\nfull%ds:%ldns \n", tp_full.tv_sec, tp_full.tv_nsec);

  fmpz_poly_clear(sk);
  fmpz_poly_clear(skinvs);
  fmpz_poly_clear(pk_p0);
  fmpz_poly_clear(pk_p1);
  fmpz_clear(q);
  fmpz_clear(neg_Delta);
  fmpz_clear(Delta_2);
  fmpz_clear(t);
  fmpz_clear(Delta);
  fmpz_clear(T);
  fmpz_poly_clear(c_00);
  fmpz_poly_clear(c_01);
  fmpz_poly_clear(c_10);
  fmpz_poly_clear(c_11);
  fmpz_poly_clear(c0);
  fmpz_poly_clear(c1);
  fmpz_poly_clear(c2);
  fmpz_poly_clear(_c0);
  fmpz_poly_clear(_c1);
  fmpz_poly_clear(m0);
  fmpz_poly_clear(m1);

  for (i = 0; i < l; i++) {
    fmpz_poly_clear(rlk_r0[i]);
    /* fmpz_poly_clear(rlk_r1[i]); */
  }
}

void test_full_ntru(int m, int i) {

  fmpz_poly_t pk_p0;
  fmpz_poly_t pk_p1;
  fmpz_poly_t sk;
  fmpz_poly_t *rlk_r0;
  fmpz_poly_t skinvs;
  fmpz_poly_t m0;
  fmpz_poly_t m1;
  fmpz_poly_t c_00;
  fmpz_poly_t c_01;
  fmpz_poly_t c_10;
  fmpz_poly_t c_11;
  fmpz_poly_t c0;
  fmpz_poly_t c1;
  fmpz_poly_t c2;
  fmpz_poly_t _c0;
  fmpz_poly_t _c1;
  long timedif;

  fmpz_init(q);
  fmpz_init(neg_Delta);
  fmpz_init(Delta_2);
  fmpz_init(t);
  fmpz_init(Delta);
  fmpz_init(T);

  fmpz_poly_init(pk_p0);
  fmpz_poly_init(pk_p1);
  fmpz_poly_init(skinvs);
  fmpz_poly_init(sk);
  fmpz_poly_init(m0);
  fmpz_poly_init(m1);
  fmpz_poly_init(c_00);
  fmpz_poly_init(c_01);
  fmpz_poly_init(c_10);
  fmpz_poly_init(c_11);
  fmpz_poly_init(c0);
  fmpz_poly_init(c1);
  fmpz_poly_init(c2);
  fmpz_poly_init(_c0);
  fmpz_poly_init(_c1);

  System_Param();

  rlk_r0 = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * (l));

  for (i = 0; i < l; i++) {
    fmpz_poly_init(rlk_r0[i]);
  }

  // test begin
  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  SH_Keygen_ntru(sk, skinvs, pk_p0, pk_p1, rlk_r0);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nKEYGEN%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  Encode(m, m0);
  Encode(m, m1);

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  SH_Encrypt(m0, pk_p0, pk_p1, c_00, c_01);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nENC%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  SH_Encrypt(m1, pk_p0, pk_p1, c_10, c_11);

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  SH_Decrypt(sk, c_00, c_01);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nDEC%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  SH_Decrypt(sk, c_10, c_11);

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  SH_Add(c_00, c_01, c_10, c_11, c0, c1);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nADD%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  SH_Decrypt(sk, c0, c1);

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  SH_Mul(c_00, c_01, c_10, c_11, c0, c1, c2);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nMUL%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  // SH_DecMul(sk, c0, c1, c2);

  /* clock_gettime(CLOCK_MONOTONIC, &tpstart); */
  /* Relinear(rlk_r0, rlk_r1, c0, c1, c2, _c0, _c1); */
  Relinear_ntru(rlk_r0, c0, c1, c2, _c0, _c1);
  /* clock_gettime(CLOCK_MONOTONIC, &tpend); */
  /* tp = diff(tpstart, tpend); */
  /* printf("\nREL%ds:%ldns \n", tp.tv_sec, tp.tv_nsec); */

  SH_Decrypt(sk, _c0, _c1);

  fmpz_poly_clear(sk);
  fmpz_poly_clear(skinvs);
  fmpz_poly_clear(pk_p0);
  fmpz_poly_clear(pk_p1);
  fmpz_clear(q);
  fmpz_clear(neg_Delta);
  fmpz_clear(Delta_2);
  fmpz_clear(t);
  fmpz_clear(Delta);
  fmpz_clear(T);
  fmpz_poly_clear(c_00);
  fmpz_poly_clear(c_01);
  fmpz_poly_clear(c_10);
  fmpz_poly_clear(c_11);
  fmpz_poly_clear(c0);
  fmpz_poly_clear(c1);
  fmpz_poly_clear(c2);
  fmpz_poly_clear(_c0);
  fmpz_poly_clear(_c1);
  fmpz_poly_clear(m0);
  fmpz_poly_clear(m1);

  for (i = 0; i < l; i++) {
    fmpz_poly_clear(rlk_r0[i]);
    /* fmpz_poly_clear(rlk_r1[i]); */
  }
}

test_ntru_and_bfv() {

  fmpz_poly_t pk_p0;
  fmpz_poly_t pk_p1;
  fmpz_poly_t sk;
  fmpz_poly_t *rlk_r0_ntru;
  fmpz_poly_t *rlk_r0;
  fmpz_poly_t *rlk_r1;
  fmpz_poly_t skinvs;
  int m;
  int i;
  fmpz_poly_t m0;
  fmpz_poly_t m1;
  fmpz_poly_t c_00;
  fmpz_poly_t c_01;
  fmpz_poly_t c_10;
  fmpz_poly_t c_11;
  fmpz_poly_t c0;
  fmpz_poly_t c1;
  fmpz_poly_t c2;
  fmpz_poly_t _c0;
  fmpz_poly_t _c1;
  struct timespec tpstart;
  struct timespec tpend;
  struct timespec tp;
  long timedif;

  fmpz_init(q);
  fmpz_init(neg_Delta);
  fmpz_init(Delta_2);
  fmpz_init(t);
  fmpz_init(Delta);
  fmpz_init(T);

  fmpz_poly_init(pk_p0);
  fmpz_poly_init(pk_p1);
  fmpz_poly_init(skinvs);
  fmpz_poly_init(sk);
  fmpz_poly_init(m0);
  fmpz_poly_init(m1);
  fmpz_poly_init(c_00);
  fmpz_poly_init(c_01);
  fmpz_poly_init(c_10);
  fmpz_poly_init(c_11);
  fmpz_poly_init(c0);
  fmpz_poly_init(c1);
  fmpz_poly_init(c2);
  fmpz_poly_init(_c0);
  fmpz_poly_init(_c1);

  System_Param();

  rlk_r0 = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * (l));
  rlk_r1 = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * (l));
  rlk_r0_ntru = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * (l));

  for (i = 0; i < l; i++) {
    fmpz_poly_init(rlk_r0[i]);
    fmpz_poly_init(rlk_r0_ntru[i]);
    fmpz_poly_init(rlk_r1[i]);
  }

  // test begin
  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  SH_Keygen_ntru(sk, skinvs, pk_p0, pk_p1, rlk_r0_ntru);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart, tpend);
  printf("\nKEYGEN%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);

  scanf("%d", &m);
  Encode(m, m0);
  scanf("%d", &m);
  Encode(m, m1);

  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  SH_Encrypt(m0, pk_p0, pk_p1, c_00, c_01);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart, tpend);
  printf("\nENC%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);

  SH_Encrypt(m1, pk_p0, pk_p1, c_10, c_11);

  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  SH_Decrypt(sk, c_00, c_01);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart, tpend);
  printf("\nDEC%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);

  SH_Decrypt(sk, c_10, c_11);

  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  SH_Add(c_00, c_01, c_10, c_11, c0, c1);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart, tpend);
  printf("\nADD%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);

  SH_Decrypt(sk, c0, c1);

  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  SH_Mul(c_00, c_01, c_10, c_11, c0, c1, c2);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart, tpend);
  printf("\nMUL%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);

  // SH_DecMul(sk, c0, c1, c2);

  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  /* Relinear(rlk_r0, rlk_r1, c0, c1, c2, _c0, _c1); */
  Relinear_ntru(rlk_r0_ntru, c0, c1, c2, _c0, _c1);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart, tpend);
  printf("\nREL%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);

  SH_Decrypt(sk, _c0, _c1);

  fmpz_poly_clear(sk);
  fmpz_poly_clear(skinvs);
  fmpz_poly_clear(pk_p0);
  fmpz_poly_clear(pk_p1);
  fmpz_clear(q);
  fmpz_clear(neg_Delta);
  fmpz_clear(Delta_2);
  fmpz_clear(t);
  fmpz_clear(Delta);
  fmpz_clear(T);
  fmpz_clear(T2);
  fmpz_poly_clear(c_00);
  fmpz_poly_clear(c_01);
  fmpz_poly_clear(c_10);
  fmpz_poly_clear(c_11);
  fmpz_poly_clear(c0);
  fmpz_poly_clear(c1);
  fmpz_poly_clear(c2);
  fmpz_poly_clear(_c0);
  fmpz_poly_clear(_c1);
  fmpz_poly_clear(m0);
  fmpz_poly_clear(m1);

  for (i = 0; i < l; i++) {
    fmpz_poly_clear(rlk_r0[i]);
    fmpz_poly_clear(rlk_r1[i]);
    fmpz_poly_clear(rlk_r0_ntru[i]);
    /* fmpz_poly_clear(rlk_r1[i]); */
  }
}

void test_full_speed(int num) {
  struct timespec tpstart;
  struct timespec tpend;
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  test_full_bfv(1, 1);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart, tpend);
  printf("\nFUll_BFV%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);

  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  test_full_ntru(1, 1);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart, tpend);
  printf("\nFULL_BFV_WITH_NTRU:%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);
}

int main() {
  flint_set_num_threads(20);
  test_bfv();
  test_ntru();
  /* test_full_speed(1); */
  /* test_ntru_and_bfv(); */
  return 0;
}
