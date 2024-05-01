#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <x86intrin.h>

//// TO BE COMPILE WITH : g++ -mavx512f -O3 11_nbody.cpp
int main() {
  
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
 
  for(int i=0; i<N; i++) {
    __m512 mvec = _mm512_load_ps(m);

    //float rx
    __m512 xjvec = _mm512_load_ps(x);
    __m512 xivec = _mm512_set1_ps(x[i]);
    //float ry
    __m512 yjvec = _mm512_load_ps(y);
    __m512 yivec = _mm512_set1_ps(y[i]);
    //rx = x[i] - x[j];
    __m512 rxvec = _mm512_sub_ps(xivec, xjvec);
    //ry = y[i] - y[j];
    __m512 ryvec = _mm512_sub_ps(yivec, yjvec);

    //(rx * rx + ry * ry)
    __m512 rxvec2 = _mm512_mul_ps(rxvec, rxvec);
    __m512 ryvec2 = _mm512_mul_ps(ryvec, ryvec);
    __m512 rvec2 = _mm512_add_ps(rxvec2, ryvec2);

    //float r = std::sqrt(rx * rx + ry * ry);
    __m512 r_vec = _mm512_rsqrt14_ps(rvec2);

    //(r * r * r)
    __m512 r2_vec = _mm512_mul_ps(r_vec, r_vec);
    __m512 r3_vec = _mm512_mul_ps(r2_vec, r_vec);
    
    //rx * m[j] / (r * r * r);
    __m512 fx_sub_vec = _mm512_mul_ps(rxvec, r3_vec);
    __m512 fx_vec = _mm512_mul_ps(fx_sub_vec, mvec);
    
    //ry * m[j] / (r * r * r);
    __m512 fy_sub_vec = _mm512_mul_ps(ryvec, r3_vec);
    __m512 fy_vec = _mm512_mul_ps(fy_sub_vec, mvec);

    //Mask initialization
    __m512 fx_mask_vec = _mm512_setzero_ps();
    __m512 fy_mask_vec = _mm512_setzero_ps();
    __m512 limit = _mm512_set1_ps(10e+5);
    
    __mmask16 mask_x = _mm512_cmp_ps_mask(fx_vec, limit, _MM_CMPINT_GT);
    fx_mask_vec = _mm512_mask_blend_ps(mask_x, fx_vec, fx_mask_vec);
    __mmask16 mask_y = _mm512_cmp_ps_mask(fy_vec, limit, _MM_CMPINT_GT);
    fy_mask_vec = _mm512_mask_blend_ps(mask_y, fy_vec, fy_mask_vec);
    
    //fx[i] -= rx * m[j] / (r * r * r);
    fx[i] = _mm512_reduce_add_ps(fx_mask_vec);
    //fy[i] -= ry * m[j] / (r * r * r);
    fy[i] = _mm512_reduce_add_ps(fy_mask_vec);

    printf("%d %g %g\n",i,fx[i],fy[i]);
  } 
}
