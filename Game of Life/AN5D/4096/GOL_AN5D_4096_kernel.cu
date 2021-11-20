#include "GOL_AN5D_4096_kernel.hu"
__device__ char __sbref_wrap(char *sb, size_t index) { return sb[index]; }

__global__ void kernel0_4(char *grid, char *rule_table, int c0)
{
#ifndef AN5D_TYPE
#define AN5D_TYPE unsigned
#endif
    const AN5D_TYPE __c0Len = (1023 - 0 + 1);
    const AN5D_TYPE __c0Pad = (0);
    #define __c0 c0
    const AN5D_TYPE __c1Len = (4096 - 1 + 1);
    const AN5D_TYPE __c1Pad = (1);
    #define __c1 c1
    const AN5D_TYPE __c2Len = (4096 - 1 + 1);
    const AN5D_TYPE __c2Pad = (1);
    #define __c2 c2
    const AN5D_TYPE __halo1 = 1;
    const AN5D_TYPE __halo2 = 1;
    const AN5D_TYPE __side0Len = 4;
    const AN5D_TYPE __side1Len = 128;
    const AN5D_TYPE __side2Len = 24;
    const AN5D_TYPE __OlLen1 = (__halo1 * __side0Len);
    const AN5D_TYPE __OlLen2 = (__halo2 * __side0Len);
    const AN5D_TYPE __side1LenOl = (__side1Len + 2 * __OlLen1);
    const AN5D_TYPE __side2LenOl = (__side2Len + 2 * __OlLen2);
    const AN5D_TYPE __blockSize = 1 * __side2LenOl;
    const AN5D_TYPE __side1Num = (__c1Len + __side1Len - 1) / __side1Len;
    const AN5D_TYPE __side2Num = (__c2Len + __side2Len - 1) / __side2Len;
    const AN5D_TYPE __tid = threadIdx.y * blockDim.x + threadIdx.x;
    const AN5D_TYPE __local_c2 = __tid;
    const AN5D_TYPE __c1Id = blockIdx.x / __side2Num;
    const AN5D_TYPE __c2 = (blockIdx.x % __side2Num) * __side2Len + __local_c2 + __c2Pad - __OlLen2;
    char __reg_0_0;
    char __reg_0_1;
    char __reg_0_2;
    char __reg_1_0;
    char __reg_1_1;
    char __reg_1_2;
    char __reg_2_0;
    char __reg_2_1;
    char __reg_2_2;
    char __reg_3_0;
    char __reg_3_1;
    char __reg_3_2;
    __shared__ char __a_sb_double[__blockSize * 2];
    char *__a_sb = __a_sb_double;
    __shared__ char __b_sb_double[__blockSize * 2];
    char *__b_sb = __b_sb_double;
    __shared__ char __c_sb_double[__blockSize * 2];
    char *__c_sb = __c_sb_double;
    const AN5D_TYPE __loadValid = 1 && __c2 >= __c2Pad - __halo2 && __c2 < __c2Pad + __c2Len + __halo2;
    const AN5D_TYPE __updateValid = 1 && __c2 >= __c2Pad && __c2 < __c2Pad + __c2Len;
    const AN5D_TYPE __writeValid1 = __updateValid && __local_c2 >= (__halo2 * 1) && __local_c2 < __side2LenOl - (__halo2 * 1);
    const AN5D_TYPE __writeValid2 = __updateValid && __local_c2 >= (__halo2 * 2) && __local_c2 < __side2LenOl - (__halo2 * 2);
    const AN5D_TYPE __writeValid3 = __updateValid && __local_c2 >= (__halo2 * 3) && __local_c2 < __side2LenOl - (__halo2 * 3);
    const AN5D_TYPE __writeValid4 = __updateValid && __local_c2 >= (__halo2 * 4) && __local_c2 < __side2LenOl - (__halo2 * 4);
    const AN5D_TYPE __storeValid = __writeValid4;
    AN5D_TYPE __c1;
    AN5D_TYPE __h;
    const AN5D_TYPE __c1Pad2 = __c1Pad + __side1Len * __c1Id;
    #define __LOAD(reg, h) do { if (__loadValid) { __c1 = __c1Pad2 - __halo1 + h; reg = grid[((__c0 % 2) * 4098 + __c1) * 4098 + __c2]; }} while (0)
    #define __DEST (grid[(((c0 + 1) % 2) * 4098 + c1) * 4098 + c2])
    #define __REGREF(reg, i2) reg
    #define __SBREF(sb, i2) __sbref_wrap(sb, (int)__tid + i2)
    #define __CALCEXPR(__rn0, __a, __b, __c) do { __rn0 = rule_table[grid[((c0 % 2) * 4098 + c1) * 4098 + c2] * 9 + (grid[((c0 % 2) * 4098 + (c1 - 1)) * 4098 + c2] + grid[((c0 % 2) * 4098 + (c1 + 1)) * 4098 + c2] + grid[((c0 % 2) * 4098 + c1) * 4098 + (c2 - 1)] + grid[((c0 % 2) * 4098 + c1) * 4098 + (c2 + 1)] + grid[((c0 % 2) * 4098 + (c1 - 1)) * 4098 + (c2 - 1)] + grid[((c0 % 2) * 4098 + (c1 - 1)) * 4098 + (c2 + 1)] + grid[((c0 % 2) * 4098 + (c1 + 1)) * 4098 + (c2 - 1)] + grid[((c0 % 2) * 4098 + (c1 + 1)) * 4098 + (c2 + 1)])]; } while (0)
    #define __DB_SWITCH() do { __a_sb = &__a_sb_double[(__a_sb == __a_sb_double) ? __blockSize : 0]; __b_sb = &__b_sb_double[(__b_sb == __b_sb_double) ? __blockSize : 0]; __c_sb = &__c_sb_double[(__c_sb == __c_sb_double) ? __blockSize : 0]; } while (0)
    #define __CALCSETUP(a, b, c) do { __DB_SWITCH(); __a_sb[__tid] = a; __b_sb[__tid] = b; __c_sb[__tid] = c; __syncthreads(); } while (0)
    #define __CALC1(out, reg0, reg1, reg2) do { __CALCSETUP(reg0, reg1, reg2); if (__writeValid1) __CALCEXPR(out, reg0, reg1, reg2); else out = reg1; } while (0)
    #define __CALC2(out, reg0, reg1, reg2) do { __CALCSETUP(reg0, reg1, reg2); if (__writeValid2) __CALCEXPR(out, reg0, reg1, reg2); else out = reg1; } while (0)
    #define __CALC3(out, reg0, reg1, reg2) do { __CALCSETUP(reg0, reg1, reg2); if (__writeValid3) __CALCEXPR(out, reg0, reg1, reg2); else out = reg1; } while (0)
    #define __STORE(h, reg0, reg1, reg2) do { __CALCSETUP(reg0, reg1, reg2); if (__storeValid) { __c1 = __c1Pad2 - __halo1 + h; __CALCEXPR(__DEST, reg0, reg1, reg2); } } while (0)
    if (__c1Id == 0)
    {
      __LOAD(__reg_3_0, 0);
      __LOAD(__reg_0_1, 1);
      __LOAD(__reg_0_2, 2);
      __CALC1(__reg_1_1, __reg_3_0, __reg_0_1, __reg_0_2);
      __LOAD(__reg_0_0, 3);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __CALC2(__reg_2_1, __reg_3_0, __reg_1_1, __reg_1_2);
      __LOAD(__reg_0_1, 4);
      __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
      __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
      __CALC3(__reg_3_1, __reg_3_0, __reg_2_1, __reg_2_2);
      __LOAD(__reg_0_2, 5);
      __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
      __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
      __CALC3(__reg_3_2, __reg_2_1, __reg_2_2, __reg_2_0);
      __STORE(1, __reg_3_0, __reg_3_1, __reg_3_2);
      __LOAD(__reg_0_0, 6);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_1_2);
      __CALC3(__reg_3_0, __reg_2_2, __reg_2_0, __reg_2_1);
      __STORE(2, __reg_3_1, __reg_3_2, __reg_3_0);
      __LOAD(__reg_0_1, 7);
      __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
      __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
      __CALC3(__reg_3_1, __reg_2_0, __reg_2_1, __reg_2_2);
      __STORE(3, __reg_3_2, __reg_3_0, __reg_3_1);
      __LOAD(__reg_0_2, 8);
      __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
      __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
      __CALC3(__reg_3_2, __reg_2_1, __reg_2_2, __reg_2_0);
      __STORE(4, __reg_3_0, __reg_3_1, __reg_3_2);
    }
    else
    {
      __LOAD(__reg_0_0, 0);
      __LOAD(__reg_0_1, 1);
      __LOAD(__reg_0_2, 2);
      __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
      __LOAD(__reg_0_0, 3);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __LOAD(__reg_0_1, 4);
      __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
      __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
      __LOAD(__reg_0_2, 5);
      __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
      __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
      __LOAD(__reg_0_0, 6);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_1_2);
      __CALC3(__reg_3_0, __reg_2_2, __reg_2_0, __reg_2_1);
      __LOAD(__reg_0_1, 7);
      __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
      __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
      __CALC3(__reg_3_1, __reg_2_0, __reg_2_1, __reg_2_2);
      __LOAD(__reg_0_2, 8);
      __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
      __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
      __CALC3(__reg_3_2, __reg_2_1, __reg_2_2, __reg_2_0);
      __STORE(4, __reg_3_0, __reg_3_1, __reg_3_2);
    }
    __a_sb = __a_sb_double + __blockSize * 0;
    __b_sb = __b_sb_double + __blockSize * 0;
    __c_sb = __c_sb_double + __blockSize * 0;
    if (__c1Id == __side1Num - 1)
    {
      for (__h = 9; __h <= __c1Len - __side1Len * __c1Id + __halo1 * 2 - 3;)
      {
        __LOAD(__reg_0_0, __h);
        __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
        __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_1_2);
        __CALC3(__reg_3_0, __reg_2_2, __reg_2_0, __reg_2_1);
        __STORE(__h - 4, __reg_3_1, __reg_3_2, __reg_3_0);
        __h++;
        __LOAD(__reg_0_1, __h);
        __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
        __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
        __CALC3(__reg_3_1, __reg_2_0, __reg_2_1, __reg_2_2);
        __STORE(__h - 4, __reg_3_2, __reg_3_0, __reg_3_1);
        __h++;
        __LOAD(__reg_0_2, __h);
        __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
        __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
        __CALC3(__reg_3_2, __reg_2_1, __reg_2_2, __reg_2_0);
        __STORE(__h - 4, __reg_3_0, __reg_3_1, __reg_3_2);
        __h++;
      }
      if (0) {}
      else if (__h + 0 == __c1Len - __side1Len * __c1Id + __halo1 * 2)
      {
        __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_0_2);
        __CALC3(__reg_3_0, __reg_2_2, __reg_2_0, __reg_2_1);
        __STORE(__h - 4, __reg_3_1, __reg_3_2, __reg_3_0);
        __CALC3(__reg_3_1, __reg_2_0, __reg_2_1, __reg_0_2);
        __STORE(__h - 3, __reg_3_2, __reg_3_0, __reg_3_1);
        __STORE(__h - 2, __reg_3_0, __reg_3_1, __reg_0_2);
      }
      else if (__h + 1 == __c1Len - __side1Len * __c1Id + __halo1 * 2)
      {
        __LOAD(__reg_0_0, __h + 0);
        __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
        __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_1_2);
        __CALC3(__reg_3_0, __reg_2_2, __reg_2_0, __reg_2_1);
        __STORE(__h - 4, __reg_3_1, __reg_3_2, __reg_3_0);
        __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_0_0);
        __CALC3(__reg_3_1, __reg_2_0, __reg_2_1, __reg_2_2);
        __STORE(__h - 3, __reg_3_2, __reg_3_0, __reg_3_1);
        __CALC3(__reg_3_2, __reg_2_1, __reg_2_2, __reg_0_0);
        __STORE(__h - 2, __reg_3_0, __reg_3_1, __reg_3_2);
        __STORE(__h - 1, __reg_3_1, __reg_3_2, __reg_0_0);
      }
      else if (__h + 2 == __c1Len - __side1Len * __c1Id + __halo1 * 2)
      {
        __LOAD(__reg_0_0, __h + 0);
        __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
        __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_1_2);
        __CALC3(__reg_3_0, __reg_2_2, __reg_2_0, __reg_2_1);
        __STORE(__h - 4, __reg_3_1, __reg_3_2, __reg_3_0);
        __LOAD(__reg_0_1, __h + 1);
        __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
        __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
        __CALC3(__reg_3_1, __reg_2_0, __reg_2_1, __reg_2_2);
        __STORE(__h - 3, __reg_3_2, __reg_3_0, __reg_3_1);
        __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_0_1);
        __CALC3(__reg_3_2, __reg_2_1, __reg_2_2, __reg_2_0);
        __STORE(__h - 2, __reg_3_0, __reg_3_1, __reg_3_2);
        __CALC3(__reg_3_0, __reg_2_2, __reg_2_0, __reg_0_1);
        __STORE(__h - 1, __reg_3_1, __reg_3_2, __reg_3_0);
        __STORE(__h + 0, __reg_3_2, __reg_3_0, __reg_0_1);
      }
    }
    else
    {
      for (__h = 9; __h <= __side1LenOl - 3;)
      {
        __LOAD(__reg_0_0, __h);
        __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
        __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_1_2);
        __CALC3(__reg_3_0, __reg_2_2, __reg_2_0, __reg_2_1);
        __STORE(__h - 4, __reg_3_1, __reg_3_2, __reg_3_0);
        __h++;
        __LOAD(__reg_0_1, __h);
        __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
        __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
        __CALC3(__reg_3_1, __reg_2_0, __reg_2_1, __reg_2_2);
        __STORE(__h - 4, __reg_3_2, __reg_3_0, __reg_3_1);
        __h++;
        __LOAD(__reg_0_2, __h);
        __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
        __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
        __CALC3(__reg_3_2, __reg_2_1, __reg_2_2, __reg_2_0);
        __STORE(__h - 4, __reg_3_0, __reg_3_1, __reg_3_2);
        __h++;
      }
      if (__h == __side1LenOl) return;
      __LOAD(__reg_0_0, __h);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_1_2);
      __CALC3(__reg_3_0, __reg_2_2, __reg_2_0, __reg_2_1);
      __STORE(__h - 4, __reg_3_1, __reg_3_2, __reg_3_0);
      __h++;
      if (__h == __side1LenOl) return;
      __LOAD(__reg_0_1, __h);
      __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
      __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
      __CALC3(__reg_3_1, __reg_2_0, __reg_2_1, __reg_2_2);
      __STORE(__h - 4, __reg_3_2, __reg_3_0, __reg_3_1);
      __h++;
      if (__h == __side1LenOl) return;
      __LOAD(__reg_0_2, __h);
      __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
      __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
      __CALC3(__reg_3_2, __reg_2_1, __reg_2_2, __reg_2_0);
      __STORE(__h - 4, __reg_3_0, __reg_3_1, __reg_3_2);
      __h++;
    }
}
__global__ void kernel0_3(char *grid, char *rule_table, int c0)
{
#ifndef AN5D_TYPE
#define AN5D_TYPE unsigned
#endif
    const AN5D_TYPE __c0Len = (1023 - 0 + 1);
    const AN5D_TYPE __c0Pad = (0);
    #define __c0 c0
    const AN5D_TYPE __c1Len = (4096 - 1 + 1);
    const AN5D_TYPE __c1Pad = (1);
    #define __c1 c1
    const AN5D_TYPE __c2Len = (4096 - 1 + 1);
    const AN5D_TYPE __c2Pad = (1);
    #define __c2 c2
    const AN5D_TYPE __halo1 = 1;
    const AN5D_TYPE __halo2 = 1;
    const AN5D_TYPE __side0Len = 3;
    const AN5D_TYPE __side1Len = 128;
    const AN5D_TYPE __side2Len = 26;
    const AN5D_TYPE __OlLen1 = (__halo1 * __side0Len);
    const AN5D_TYPE __OlLen2 = (__halo2 * __side0Len);
    const AN5D_TYPE __side1LenOl = (__side1Len + 2 * __OlLen1);
    const AN5D_TYPE __side2LenOl = (__side2Len + 2 * __OlLen2);
    const AN5D_TYPE __blockSize = 1 * __side2LenOl;
    const AN5D_TYPE __side1Num = (__c1Len + __side1Len - 1) / __side1Len;
    const AN5D_TYPE __side2Num = (__c2Len + __side2Len - 1) / __side2Len;
    const AN5D_TYPE __tid = threadIdx.y * blockDim.x + threadIdx.x;
    const AN5D_TYPE __local_c2 = __tid;
    const AN5D_TYPE __c1Id = blockIdx.x / __side2Num;
    const AN5D_TYPE __c2 = (blockIdx.x % __side2Num) * __side2Len + __local_c2 + __c2Pad - __OlLen2;
    char __reg_0_0;
    char __reg_0_1;
    char __reg_0_2;
    char __reg_1_0;
    char __reg_1_1;
    char __reg_1_2;
    char __reg_2_0;
    char __reg_2_1;
    char __reg_2_2;
    __shared__ char __a_sb_double[__blockSize * 2];
    char *__a_sb = __a_sb_double;
    __shared__ char __b_sb_double[__blockSize * 2];
    char *__b_sb = __b_sb_double;
    __shared__ char __c_sb_double[__blockSize * 2];
    char *__c_sb = __c_sb_double;
    const AN5D_TYPE __loadValid = 1 && __c2 >= __c2Pad - __halo2 && __c2 < __c2Pad + __c2Len + __halo2;
    const AN5D_TYPE __updateValid = 1 && __c2 >= __c2Pad && __c2 < __c2Pad + __c2Len;
    const AN5D_TYPE __writeValid1 = __updateValid && __local_c2 >= (__halo2 * 1) && __local_c2 < __side2LenOl - (__halo2 * 1);
    const AN5D_TYPE __writeValid2 = __updateValid && __local_c2 >= (__halo2 * 2) && __local_c2 < __side2LenOl - (__halo2 * 2);
    const AN5D_TYPE __writeValid3 = __updateValid && __local_c2 >= (__halo2 * 3) && __local_c2 < __side2LenOl - (__halo2 * 3);
    const AN5D_TYPE __storeValid = __writeValid3;
    AN5D_TYPE __c1;
    AN5D_TYPE __h;
    const AN5D_TYPE __c1Pad2 = __c1Pad + __side1Len * __c1Id;
    #define __LOAD(reg, h) do { if (__loadValid) { __c1 = __c1Pad2 - __halo1 + h; reg = grid[((__c0 % 2) * 4098 + __c1) * 4098 + __c2]; }} while (0)
    #define __DEST (grid[(((c0 + 1) % 2) * 4098 + c1) * 4098 + c2])
    #define __REGREF(reg, i2) reg
    #define __SBREF(sb, i2) __sbref_wrap(sb, (int)__tid + i2)
    #define __CALCEXPR(__rn0, __a, __b, __c) do { __rn0 = rule_table[grid[((c0 % 2) * 4098 + c1) * 4098 + c2] * 9 + (grid[((c0 % 2) * 4098 + (c1 - 1)) * 4098 + c2] + grid[((c0 % 2) * 4098 + (c1 + 1)) * 4098 + c2] + grid[((c0 % 2) * 4098 + c1) * 4098 + (c2 - 1)] + grid[((c0 % 2) * 4098 + c1) * 4098 + (c2 + 1)] + grid[((c0 % 2) * 4098 + (c1 - 1)) * 4098 + (c2 - 1)] + grid[((c0 % 2) * 4098 + (c1 - 1)) * 4098 + (c2 + 1)] + grid[((c0 % 2) * 4098 + (c1 + 1)) * 4098 + (c2 - 1)] + grid[((c0 % 2) * 4098 + (c1 + 1)) * 4098 + (c2 + 1)])]; } while (0)
    #define __DB_SWITCH() do { __a_sb = &__a_sb_double[(__a_sb == __a_sb_double) ? __blockSize : 0]; __b_sb = &__b_sb_double[(__b_sb == __b_sb_double) ? __blockSize : 0]; __c_sb = &__c_sb_double[(__c_sb == __c_sb_double) ? __blockSize : 0]; } while (0)
    #define __CALCSETUP(a, b, c) do { __DB_SWITCH(); __a_sb[__tid] = a; __b_sb[__tid] = b; __c_sb[__tid] = c; __syncthreads(); } while (0)
    #define __CALC1(out, reg0, reg1, reg2) do { __CALCSETUP(reg0, reg1, reg2); if (__writeValid1) __CALCEXPR(out, reg0, reg1, reg2); else out = reg1; } while (0)
    #define __CALC2(out, reg0, reg1, reg2) do { __CALCSETUP(reg0, reg1, reg2); if (__writeValid2) __CALCEXPR(out, reg0, reg1, reg2); else out = reg1; } while (0)
    #define __STORE(h, reg0, reg1, reg2) do { __CALCSETUP(reg0, reg1, reg2); if (__storeValid) { __c1 = __c1Pad2 - __halo1 + h; __CALCEXPR(__DEST, reg0, reg1, reg2); } } while (0)
    if (__c1Id == 0)
    {
      __LOAD(__reg_2_0, 0);
      __LOAD(__reg_0_1, 1);
      __LOAD(__reg_0_2, 2);
      __CALC1(__reg_1_1, __reg_2_0, __reg_0_1, __reg_0_2);
      __LOAD(__reg_0_0, 3);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __CALC2(__reg_2_1, __reg_2_0, __reg_1_1, __reg_1_2);
      __LOAD(__reg_0_1, 4);
      __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
      __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
      __STORE(1, __reg_2_0, __reg_2_1, __reg_2_2);
      __LOAD(__reg_0_2, 5);
      __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
      __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
      __STORE(2, __reg_2_1, __reg_2_2, __reg_2_0);
      __LOAD(__reg_0_0, 6);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_1_2);
      __STORE(3, __reg_2_2, __reg_2_0, __reg_2_1);
    }
    else
    {
      __LOAD(__reg_0_0, 0);
      __LOAD(__reg_0_1, 1);
      __LOAD(__reg_0_2, 2);
      __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
      __LOAD(__reg_0_0, 3);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __LOAD(__reg_0_1, 4);
      __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
      __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
      __LOAD(__reg_0_2, 5);
      __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
      __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
      __LOAD(__reg_0_0, 6);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_1_2);
      __STORE(3, __reg_2_2, __reg_2_0, __reg_2_1);
      __DB_SWITCH(); __syncthreads();
    }
    __a_sb = __a_sb_double + __blockSize * 0;
    __b_sb = __b_sb_double + __blockSize * 0;
    __c_sb = __c_sb_double + __blockSize * 0;
    if (__c1Id == __side1Num - 1)
    {
      for (__h = 7; __h <= __c1Len - __side1Len * __c1Id + __halo1 * 2 - 3;)
      {
        __LOAD(__reg_0_1, __h);
        __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
        __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
        __STORE(__h - 3, __reg_2_0, __reg_2_1, __reg_2_2);
        __h++;
        __LOAD(__reg_0_2, __h);
        __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
        __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
        __STORE(__h - 3, __reg_2_1, __reg_2_2, __reg_2_0);
        __h++;
        __LOAD(__reg_0_0, __h);
        __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
        __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_1_2);
        __STORE(__h - 3, __reg_2_2, __reg_2_0, __reg_2_1);
        __h++;
        __DB_SWITCH(); __syncthreads();
      }
      if (0) {}
      else if (__h + 0 == __c1Len - __side1Len * __c1Id + __halo1 * 2)
      {
        __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_0_0);
        __STORE(__h - 3, __reg_2_0, __reg_2_1, __reg_2_2);
        __STORE(__h - 2, __reg_2_1, __reg_2_2, __reg_0_0);
      }
      else if (__h + 1 == __c1Len - __side1Len * __c1Id + __halo1 * 2)
      {
        __LOAD(__reg_0_1, __h + 0);
        __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
        __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
        __STORE(__h - 3, __reg_2_0, __reg_2_1, __reg_2_2);
        __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_0_1);
        __STORE(__h - 2, __reg_2_1, __reg_2_2, __reg_2_0);
        __STORE(__h - 1, __reg_2_2, __reg_2_0, __reg_0_1);
      }
      else if (__h + 2 == __c1Len - __side1Len * __c1Id + __halo1 * 2)
      {
        __LOAD(__reg_0_1, __h + 0);
        __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
        __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
        __STORE(__h - 3, __reg_2_0, __reg_2_1, __reg_2_2);
        __LOAD(__reg_0_2, __h + 1);
        __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
        __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
        __STORE(__h - 2, __reg_2_1, __reg_2_2, __reg_2_0);
        __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_0_2);
        __STORE(__h - 1, __reg_2_2, __reg_2_0, __reg_2_1);
        __STORE(__h + 0, __reg_2_0, __reg_2_1, __reg_0_2);
      }
    }
    else
    {
      for (__h = 7; __h <= __side1LenOl - 3;)
      {
        __LOAD(__reg_0_1, __h);
        __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
        __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
        __STORE(__h - 3, __reg_2_0, __reg_2_1, __reg_2_2);
        __h++;
        __LOAD(__reg_0_2, __h);
        __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
        __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
        __STORE(__h - 3, __reg_2_1, __reg_2_2, __reg_2_0);
        __h++;
        __LOAD(__reg_0_0, __h);
        __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
        __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_1_2);
        __STORE(__h - 3, __reg_2_2, __reg_2_0, __reg_2_1);
        __h++;
        __DB_SWITCH();  __syncthreads();
      }
      if (__h == __side1LenOl) return;
      __LOAD(__reg_0_1, __h);
      __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
      __CALC2(__reg_2_2, __reg_1_1, __reg_1_2, __reg_1_0);
      __STORE(__h - 3, __reg_2_0, __reg_2_1, __reg_2_2);
      __h++;
      if (__h == __side1LenOl) return;
      __LOAD(__reg_0_2, __h);
      __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
      __CALC2(__reg_2_0, __reg_1_2, __reg_1_0, __reg_1_1);
      __STORE(__h - 3, __reg_2_1, __reg_2_2, __reg_2_0);
      __h++;
      if (__h == __side1LenOl) return;
      __LOAD(__reg_0_0, __h);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __CALC2(__reg_2_1, __reg_1_0, __reg_1_1, __reg_1_2);
      __STORE(__h - 3, __reg_2_2, __reg_2_0, __reg_2_1);
      __h++;
    }
}
__global__ void kernel0_2(char *grid, char *rule_table, int c0)
{
#ifndef AN5D_TYPE
#define AN5D_TYPE unsigned
#endif
    const AN5D_TYPE __c0Len = (1023 - 0 + 1);
    const AN5D_TYPE __c0Pad = (0);
    #define __c0 c0
    const AN5D_TYPE __c1Len = (4096 - 1 + 1);
    const AN5D_TYPE __c1Pad = (1);
    #define __c1 c1
    const AN5D_TYPE __c2Len = (4096 - 1 + 1);
    const AN5D_TYPE __c2Pad = (1);
    #define __c2 c2
    const AN5D_TYPE __halo1 = 1;
    const AN5D_TYPE __halo2 = 1;
    const AN5D_TYPE __side0Len = 2;
    const AN5D_TYPE __side1Len = 128;
    const AN5D_TYPE __side2Len = 28;
    const AN5D_TYPE __OlLen1 = (__halo1 * __side0Len);
    const AN5D_TYPE __OlLen2 = (__halo2 * __side0Len);
    const AN5D_TYPE __side1LenOl = (__side1Len + 2 * __OlLen1);
    const AN5D_TYPE __side2LenOl = (__side2Len + 2 * __OlLen2);
    const AN5D_TYPE __blockSize = 1 * __side2LenOl;
    const AN5D_TYPE __side1Num = (__c1Len + __side1Len - 1) / __side1Len;
    const AN5D_TYPE __side2Num = (__c2Len + __side2Len - 1) / __side2Len;
    const AN5D_TYPE __tid = threadIdx.y * blockDim.x + threadIdx.x;
    const AN5D_TYPE __local_c2 = __tid;
    const AN5D_TYPE __c1Id = blockIdx.x / __side2Num;
    const AN5D_TYPE __c2 = (blockIdx.x % __side2Num) * __side2Len + __local_c2 + __c2Pad - __OlLen2;
    char __reg_0_0;
    char __reg_0_1;
    char __reg_0_2;
    char __reg_1_0;
    char __reg_1_1;
    char __reg_1_2;
    __shared__ char __a_sb_double[__blockSize * 2];
    char *__a_sb = __a_sb_double;
    __shared__ char __b_sb_double[__blockSize * 2];
    char *__b_sb = __b_sb_double;
    __shared__ char __c_sb_double[__blockSize * 2];
    char *__c_sb = __c_sb_double;
    const AN5D_TYPE __loadValid = 1 && __c2 >= __c2Pad - __halo2 && __c2 < __c2Pad + __c2Len + __halo2;
    const AN5D_TYPE __updateValid = 1 && __c2 >= __c2Pad && __c2 < __c2Pad + __c2Len;
    const AN5D_TYPE __writeValid1 = __updateValid && __local_c2 >= (__halo2 * 1) && __local_c2 < __side2LenOl - (__halo2 * 1);
    const AN5D_TYPE __writeValid2 = __updateValid && __local_c2 >= (__halo2 * 2) && __local_c2 < __side2LenOl - (__halo2 * 2);
    const AN5D_TYPE __storeValid = __writeValid2;
    AN5D_TYPE __c1;
    AN5D_TYPE __h;
    const AN5D_TYPE __c1Pad2 = __c1Pad + __side1Len * __c1Id;
    #define __LOAD(reg, h) do { if (__loadValid) { __c1 = __c1Pad2 - __halo1 + h; reg = grid[((__c0 % 2) * 4098 + __c1) * 4098 + __c2]; }} while (0)
    #define __DEST (grid[(((c0 + 1) % 2) * 4098 + c1) * 4098 + c2])
    #define __REGREF(reg, i2) reg
    #define __SBREF(sb, i2) __sbref_wrap(sb, (int)__tid + i2)
    #define __CALCEXPR(__rn0, __a, __b, __c) do { __rn0 = rule_table[grid[((c0 % 2) * 4098 + c1) * 4098 + c2] * 9 + (grid[((c0 % 2) * 4098 + (c1 - 1)) * 4098 + c2] + grid[((c0 % 2) * 4098 + (c1 + 1)) * 4098 + c2] + grid[((c0 % 2) * 4098 + c1) * 4098 + (c2 - 1)] + grid[((c0 % 2) * 4098 + c1) * 4098 + (c2 + 1)] + grid[((c0 % 2) * 4098 + (c1 - 1)) * 4098 + (c2 - 1)] + grid[((c0 % 2) * 4098 + (c1 - 1)) * 4098 + (c2 + 1)] + grid[((c0 % 2) * 4098 + (c1 + 1)) * 4098 + (c2 - 1)] + grid[((c0 % 2) * 4098 + (c1 + 1)) * 4098 + (c2 + 1)])]; } while (0)
    #define __DB_SWITCH() do { __a_sb = &__a_sb_double[(__a_sb == __a_sb_double) ? __blockSize : 0]; __b_sb = &__b_sb_double[(__b_sb == __b_sb_double) ? __blockSize : 0]; __c_sb = &__c_sb_double[(__c_sb == __c_sb_double) ? __blockSize : 0]; } while (0)
    #define __CALCSETUP(a, b, c) do { __DB_SWITCH(); __a_sb[__tid] = a; __b_sb[__tid] = b; __c_sb[__tid] = c; __syncthreads(); } while (0)
    #define __CALC1(out, reg0, reg1, reg2) do { __CALCSETUP(reg0, reg1, reg2); if (__writeValid1) __CALCEXPR(out, reg0, reg1, reg2); else out = reg1; } while (0)
    #define __STORE(h, reg0, reg1, reg2) do { __CALCSETUP(reg0, reg1, reg2); if (__storeValid) { __c1 = __c1Pad2 - __halo1 + h; __CALCEXPR(__DEST, reg0, reg1, reg2); } } while (0)
    if (__c1Id == 0)
    {
      __LOAD(__reg_1_0, 0);
      __LOAD(__reg_0_1, 1);
      __LOAD(__reg_0_2, 2);
      __CALC1(__reg_1_1, __reg_1_0, __reg_0_1, __reg_0_2);
      __LOAD(__reg_0_0, 3);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __STORE(1, __reg_1_0, __reg_1_1, __reg_1_2);
      __LOAD(__reg_0_1, 4);
      __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
      __STORE(2, __reg_1_1, __reg_1_2, __reg_1_0);
    }
    else
    {
      __LOAD(__reg_0_0, 0);
      __LOAD(__reg_0_1, 1);
      __LOAD(__reg_0_2, 2);
      __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
      __LOAD(__reg_0_0, 3);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __LOAD(__reg_0_1, 4);
      __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
      __STORE(2, __reg_1_1, __reg_1_2, __reg_1_0);
      __DB_SWITCH(); __syncthreads();
    }
    __a_sb = __a_sb_double + __blockSize * 1;
    __b_sb = __b_sb_double + __blockSize * 1;
    __c_sb = __c_sb_double + __blockSize * 1;
    if (__c1Id == __side1Num - 1)
    {
      for (__h = 5; __h <= __c1Len - __side1Len * __c1Id + __halo1 * 2 - 3;)
      {
        __LOAD(__reg_0_2, __h);
        __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
        __STORE(__h - 2, __reg_1_2, __reg_1_0, __reg_1_1);
        __h++;
        __LOAD(__reg_0_0, __h);
        __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
        __STORE(__h - 2, __reg_1_0, __reg_1_1, __reg_1_2);
        __h++;
        __LOAD(__reg_0_1, __h);
        __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
        __STORE(__h - 2, __reg_1_1, __reg_1_2, __reg_1_0);
        __h++;
      }
      if (0) {}
      else if (__h + 0 == __c1Len - __side1Len * __c1Id + __halo1 * 2)
      {
        __STORE(__h - 2, __reg_1_2, __reg_1_0, __reg_0_1);
      }
      else if (__h + 1 == __c1Len - __side1Len * __c1Id + __halo1 * 2)
      {
        __LOAD(__reg_0_2, __h + 0);
        __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
        __STORE(__h - 2, __reg_1_2, __reg_1_0, __reg_1_1);
        __STORE(__h - 1, __reg_1_0, __reg_1_1, __reg_0_2);
      }
      else if (__h + 2 == __c1Len - __side1Len * __c1Id + __halo1 * 2)
      {
        __LOAD(__reg_0_2, __h + 0);
        __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
        __STORE(__h - 2, __reg_1_2, __reg_1_0, __reg_1_1);
        __LOAD(__reg_0_0, __h + 1);
        __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
        __STORE(__h - 1, __reg_1_0, __reg_1_1, __reg_1_2);
        __STORE(__h + 0, __reg_1_1, __reg_1_2, __reg_0_0);
      }
    }
    else
    {
      for (__h = 5; __h <= __side1LenOl - 3;)
      {
        __LOAD(__reg_0_2, __h);
        __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
        __STORE(__h - 2, __reg_1_2, __reg_1_0, __reg_1_1);
        __h++;
        __LOAD(__reg_0_0, __h);
        __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
        __STORE(__h - 2, __reg_1_0, __reg_1_1, __reg_1_2);
        __h++;
        __LOAD(__reg_0_1, __h);
        __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
        __STORE(__h - 2, __reg_1_1, __reg_1_2, __reg_1_0);
        __h++;
      }
      if (__h == __side1LenOl) return;
      __LOAD(__reg_0_2, __h);
      __CALC1(__reg_1_1, __reg_0_0, __reg_0_1, __reg_0_2);
      __STORE(__h - 2, __reg_1_2, __reg_1_0, __reg_1_1);
      __h++;
      if (__h == __side1LenOl) return;
      __LOAD(__reg_0_0, __h);
      __CALC1(__reg_1_2, __reg_0_1, __reg_0_2, __reg_0_0);
      __STORE(__h - 2, __reg_1_0, __reg_1_1, __reg_1_2);
      __h++;
      if (__h == __side1LenOl) return;
      __LOAD(__reg_0_1, __h);
      __CALC1(__reg_1_0, __reg_0_2, __reg_0_0, __reg_0_1);
      __STORE(__h - 2, __reg_1_1, __reg_1_2, __reg_1_0);
      __h++;
    }
}
__global__ void kernel0_1(char *grid, char *rule_table, int c0)
{
#ifndef AN5D_TYPE
#define AN5D_TYPE unsigned
#endif
    const AN5D_TYPE __c0Len = (1023 - 0 + 1);
    const AN5D_TYPE __c0Pad = (0);
    #define __c0 c0
    const AN5D_TYPE __c1Len = (4096 - 1 + 1);
    const AN5D_TYPE __c1Pad = (1);
    #define __c1 c1
    const AN5D_TYPE __c2Len = (4096 - 1 + 1);
    const AN5D_TYPE __c2Pad = (1);
    #define __c2 c2
    const AN5D_TYPE __halo1 = 1;
    const AN5D_TYPE __halo2 = 1;
    const AN5D_TYPE __side0Len = 1;
    const AN5D_TYPE __side1Len = 128;
    const AN5D_TYPE __side2Len = 30;
    const AN5D_TYPE __OlLen1 = (__halo1 * __side0Len);
    const AN5D_TYPE __OlLen2 = (__halo2 * __side0Len);
    const AN5D_TYPE __side1LenOl = (__side1Len + 2 * __OlLen1);
    const AN5D_TYPE __side2LenOl = (__side2Len + 2 * __OlLen2);
    const AN5D_TYPE __blockSize = 1 * __side2LenOl;
    const AN5D_TYPE __side1Num = (__c1Len + __side1Len - 1) / __side1Len;
    const AN5D_TYPE __side2Num = (__c2Len + __side2Len - 1) / __side2Len;
    const AN5D_TYPE __tid = threadIdx.y * blockDim.x + threadIdx.x;
    const AN5D_TYPE __local_c2 = __tid;
    const AN5D_TYPE __c1Id = blockIdx.x / __side2Num;
    const AN5D_TYPE __c2 = (blockIdx.x % __side2Num) * __side2Len + __local_c2 + __c2Pad - __OlLen2;
    char __reg_0_0;
    char __reg_0_1;
    char __reg_0_2;
    __shared__ char __a_sb_double[__blockSize * 2];
    char *__a_sb = __a_sb_double;
    __shared__ char __b_sb_double[__blockSize * 2];
    char *__b_sb = __b_sb_double;
    __shared__ char __c_sb_double[__blockSize * 2];
    char *__c_sb = __c_sb_double;
    const AN5D_TYPE __loadValid = 1 && __c2 >= __c2Pad - __halo2 && __c2 < __c2Pad + __c2Len + __halo2;
    const AN5D_TYPE __updateValid = 1 && __c2 >= __c2Pad && __c2 < __c2Pad + __c2Len;
    const AN5D_TYPE __writeValid1 = __updateValid && __local_c2 >= (__halo2 * 1) && __local_c2 < __side2LenOl - (__halo2 * 1);
    const AN5D_TYPE __storeValid = __writeValid1;
    AN5D_TYPE __c1;
    AN5D_TYPE __h;
    const AN5D_TYPE __c1Pad2 = __c1Pad + __side1Len * __c1Id;
    #define __LOAD(reg, h) do { if (__loadValid) { __c1 = __c1Pad2 - __halo1 + h; reg = grid[((__c0 % 2) * 4098 + __c1) * 4098 + __c2]; }} while (0)
    #define __DEST (grid[(((c0 + 1) % 2) * 4098 + c1) * 4098 + c2])
    #define __REGREF(reg, i2) reg
    #define __SBREF(sb, i2) __sbref_wrap(sb, (int)__tid + i2)
    #define __CALCEXPR(__rn0, __a, __b, __c) do { __rn0 = rule_table[grid[((c0 % 2) * 4098 + c1) * 4098 + c2] * 9 + (grid[((c0 % 2) * 4098 + (c1 - 1)) * 4098 + c2] + grid[((c0 % 2) * 4098 + (c1 + 1)) * 4098 + c2] + grid[((c0 % 2) * 4098 + c1) * 4098 + (c2 - 1)] + grid[((c0 % 2) * 4098 + c1) * 4098 + (c2 + 1)] + grid[((c0 % 2) * 4098 + (c1 - 1)) * 4098 + (c2 - 1)] + grid[((c0 % 2) * 4098 + (c1 - 1)) * 4098 + (c2 + 1)] + grid[((c0 % 2) * 4098 + (c1 + 1)) * 4098 + (c2 - 1)] + grid[((c0 % 2) * 4098 + (c1 + 1)) * 4098 + (c2 + 1)])]; } while (0)
    #define __DB_SWITCH() do { __a_sb = &__a_sb_double[(__a_sb == __a_sb_double) ? __blockSize : 0]; __b_sb = &__b_sb_double[(__b_sb == __b_sb_double) ? __blockSize : 0]; __c_sb = &__c_sb_double[(__c_sb == __c_sb_double) ? __blockSize : 0]; } while (0)
    #define __CALCSETUP(a, b, c) do { __DB_SWITCH(); __a_sb[__tid] = a; __b_sb[__tid] = b; __c_sb[__tid] = c; __syncthreads(); } while (0)
    #define __STORE(h, reg0, reg1, reg2) do { __CALCSETUP(reg0, reg1, reg2); if (__storeValid) { __c1 = __c1Pad2 - __halo1 + h; __CALCEXPR(__DEST, reg0, reg1, reg2); } } while (0)
    if (__c1Id == 0)
    {
      __LOAD(__reg_0_0, 0);
      __LOAD(__reg_0_1, 1);
      __LOAD(__reg_0_2, 2);
      __STORE(1, __reg_0_0, __reg_0_1, __reg_0_2);
    }
    else
    {
      __LOAD(__reg_0_0, 0);
      __LOAD(__reg_0_1, 1);
      __LOAD(__reg_0_2, 2);
      __STORE(1, __reg_0_0, __reg_0_1, __reg_0_2);
    }
    __a_sb = __a_sb_double + __blockSize * 1;
    __b_sb = __b_sb_double + __blockSize * 1;
    __c_sb = __c_sb_double + __blockSize * 1;
    if (__c1Id == __side1Num - 1)
    {
      for (__h = 3; __h <= __c1Len - __side1Len * __c1Id + __halo1 * 2 - 3;)
      {
        __LOAD(__reg_0_0, __h);
        __STORE(__h - 1, __reg_0_1, __reg_0_2, __reg_0_0);
        __h++;
        __LOAD(__reg_0_1, __h);
        __STORE(__h - 1, __reg_0_2, __reg_0_0, __reg_0_1);
        __h++;
        __LOAD(__reg_0_2, __h);
        __STORE(__h - 1, __reg_0_0, __reg_0_1, __reg_0_2);
        __h++;
        __DB_SWITCH(); __syncthreads();
      }
      if (0) {}
      else if (__h + 0 == __c1Len - __side1Len * __c1Id + __halo1 * 2)
      {
      }
      else if (__h + 1 == __c1Len - __side1Len * __c1Id + __halo1 * 2)
      {
        __LOAD(__reg_0_0, __h + 0);
        __STORE(__h - 1, __reg_0_1, __reg_0_2, __reg_0_0);
      }
      else if (__h + 2 == __c1Len - __side1Len * __c1Id + __halo1 * 2)
      {
        __LOAD(__reg_0_0, __h + 0);
        __STORE(__h - 1, __reg_0_1, __reg_0_2, __reg_0_0);
        __LOAD(__reg_0_1, __h + 1);
        __STORE(__h + 0, __reg_0_2, __reg_0_0, __reg_0_1);
      }
    }
    else
    {
      for (__h = 3; __h <= __side1LenOl - 3;)
      {
        __LOAD(__reg_0_0, __h);
        __STORE(__h - 1, __reg_0_1, __reg_0_2, __reg_0_0);
        __h++;
        __LOAD(__reg_0_1, __h);
        __STORE(__h - 1, __reg_0_2, __reg_0_0, __reg_0_1);
        __h++;
        __LOAD(__reg_0_2, __h);
        __STORE(__h - 1, __reg_0_0, __reg_0_1, __reg_0_2);
        __h++;
        __DB_SWITCH();  __syncthreads();
      }
      if (__h == __side1LenOl) return;
      __LOAD(__reg_0_0, __h);
      __STORE(__h - 1, __reg_0_1, __reg_0_2, __reg_0_0);
      __h++;
      if (__h == __side1LenOl) return;
      __LOAD(__reg_0_1, __h);
      __STORE(__h - 1, __reg_0_2, __reg_0_0, __reg_0_1);
      __h++;
      if (__h == __side1LenOl) return;
      __LOAD(__reg_0_2, __h);
      __STORE(__h - 1, __reg_0_0, __reg_0_1, __reg_0_2);
      __h++;
    }
}
