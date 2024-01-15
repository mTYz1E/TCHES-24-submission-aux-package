/*
 * MIT License
 *
 * Copyright (c) 2021-2022: imec-COSIC KU Leuven, 3001 Leuven, Belgium 
 * Author:  Michiel Van Beirendonck <michiel.vanbeirendonck@esat.kuleuven.be>
 *          Jan-Pieter D'Anvers <janpieter.danvers@esat.kuleuven.be>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include "MaskedComparison.h"
#include "randombytes.h"
#include "params.h"
#include "bitmask.h"
#include "hal.h"

#include <stdint.h>
#include <stdio.h>
#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

static void get_rand_share(size_t ncoeffs, uint32_t mod, uint32_t* x)
{
    for (size_t j = 0; j < ncoeffs; j++)
    {
        x[j] = random_uint32() % mod;
    }
}

static void mask(size_t nshares, size_t ncoeffs, uint32_t x_masked[ncoeffs][nshares], uint32_t x[ncoeffs])
{
    for (size_t i = 0; i < ncoeffs; i++)
    {
        x_masked[i][0] = x[i];

        for (size_t j = 1; j < nshares; j++)
        {
            uint32_t R = random_uint32() % Q;
            x_masked[i][0] = (x_masked[i][0] + (Q - R)) % Q;
            x_masked[i][j] = R;
        }
    #ifdef DEBUG

        uint32_t x_unmasked = 0;

        for (size_t j = 0; j < nshares; j++)
        {
            x_unmasked += x_masked[i][j];
        }

        assert((x_unmasked % Q) == x[i]);

    #endif
    }
}

#ifdef SABER
static void compress(size_t ncoeffs, size_t compressfrom, size_t compressto, uint32_t submitted_poly[ncoeffs])
{
    for (size_t i = 0; i < ncoeffs; i++)
    {
        submitted_poly[i] >>= (compressfrom - compressto);
    }
}
#endif

#ifdef KYBER
static void compress(size_t ncoeffs, __attribute__((unused)) size_t compressfrom, size_t compressto, uint32_t submitted_poly[ncoeffs])
{
    for (size_t i = 0; i < ncoeffs; i++)
    {
        submitted_poly[i] = (((submitted_poly[i] << compressto) + Q/2) / Q) % (1 << compressto);
    }
}
#endif

#ifdef KYBER
static void shared_compress(size_t ncoeffs, size_t compressto, uint32_t B[ncoeffs][NSHARES])
{
    for (size_t i = 0; i < ncoeffs; i++)
    {
    #ifdef DEBUG

        uint32_t B_unmasked = 0, BA_unmasked = 0;

        for (size_t j = 0; j < NSHARES; j++)
        {
            B_unmasked = (B_unmasked + B[i][j]) % Q;
        }

        B_unmasked = (((B_unmasked << compressto) + Q/2) / Q) & bit_mask(compressto);

    #endif

        for (size_t j = 0; j < NSHARES; j++)
        {
            uint64_t tmp = (uint64_t)B[i][j] << (compressto + KYBER_FRAC_BITS);

            if (j == 0)
            {
                tmp += (Q << KYBER_FRAC_BITS)/2;
            }

            B[i][j] = (tmp / Q); // ! make sure (uint64_t/constant) is constant-time on your platform
        }

    #ifdef DEBUG

        for (size_t j = 0; j < NSHARES; j++)
        {
            BA_unmasked =  (BA_unmasked + B[i][j]) & bit_mask(compressto + KYBER_FRAC_BITS);
        }

        BA_unmasked >>= KYBER_FRAC_BITS;

        assert(B_unmasked== BA_unmasked);

    #endif
    }
}
#endif

uint32_t* public_B;
uint32_t* public_C;
uint32_t* B;
uint32_t* C;

uint32_t* BC_Bitsliced;
uint32_t* public_B_rand;
uint32_t* public_C_rand;

// for sanity check
uint32_t* g_ct1;
uint32_t* g_ct2;
uint32_t* g_ct_diff;

static inline void
concat_ct(uint32_t* out, const uint32_t* u, const uint32_t* v) {
    size_t offset = sizeof(uint32_t) * NCOEFFS_B;
    for(int32_t i = 0; i < NCOEFFS_B; ++i) out[i] = u[i];
    for(int32_t i = 0; i < NCOEFFS_C; ++i) out[NCOEFFS_B + i] = v[i];
}

static inline void
unmask_ct(uint32_t* out, const uint32_t* MB, const uint32_t* MC) {
    for(int32_t i = 0; i < NCOEFFS_B; ++i) {
        uint32_t sum = 0;
        for(int32_t j = 0; j < NSHARES; ++j) {
            sum += MB[i * NSHARES + j];
        }

        out[i] = sum % Q;
    }

    for(int32_t i = 0; i < NCOEFFS_C; ++i) {
        uint32_t sum = 0;
        for(int32_t j = 0; j < NSHARES; ++j) {
            sum += MC[i * NSHARES + j];
        }

        out[NCOEFFS_B + i] = sum % Q;
    }
}

static inline void
diff_ct(uint32_t* out, const uint32_t* ct1, const uint32_t* ct2, int32_t len) {
    for(int32_t i = 0; i < len; ++i) {
        out[i] = abs(ct1[i] - ct2[i]) % Q;
    }
}

static inline void
print_ct(const uint32_t* ct, int32_t len) {
    assert(len >= 1);
    if(ct[0])
        printf("[ %u", ct[0]);
    else
        printf(", ");

    for(int32_t i = 1; i < len; ++i) {
        if(ct[i])
            printf(", %u", ct[i]);
        else
            printf(", ");
    }
    printf(" ]\n");
}

static inline void
prep_ct_no_dec() {
  // new ciphertext
  get_rand_share(NCOEFFS_B, Q, public_B);
  get_rand_share(NCOEFFS_C, P, public_C);

  // create reencrypted poly as a (correct) sharing of public poly
  mask(NSHARES, NCOEFFS_B, B, public_B);
  mask(NSHARES, NCOEFFS_C, C, public_C);

  public_C[0] += Q / 4;

  // Now compress public poly
  compress(NCOEFFS_B, COMPRESSFROM_B, COMPRESSTO_B, public_B);
  compress(NCOEFFS_C, COMPRESSFROM_C, COMPRESSTO_C, public_C);

  unmask_ct(g_ct1, B, C);
  compress(NCOEFFS_B, COMPRESSFROM_B, COMPRESSTO_B, g_ct1);
  compress(NCOEFFS_C, COMPRESSFROM_C, COMPRESSTO_C, g_ct1 + NCOEFFS_B);

  concat_ct(g_ct2, public_B, public_C);
  diff_ct(g_ct_diff, g_ct1, g_ct2, NCOEFFS_B + NCOEFFS_C);  
  print_ct(g_ct_diff, NCOEFFS_B + NCOEFFS_C);
}

static inline void
prep_ct_dec() {
  // first steps are the same as prep_ct_no_dec
  // new ciphertext
  get_rand_share(NCOEFFS_B, Q, public_B);
  get_rand_share(NCOEFFS_C, P, public_C);

  // NOTE: no need to compute this. B and C are overwritten later
  // create reencrypted poly as a (correct) sharing of public poly
  //mask(NSHARES, NCOEFFS_B, B, public_B);
  //mask(NSHARES, NCOEFFS_C, C, public_C);
  
  public_C[0] += Q / 4;

  // Now compress public poly
  compress(NCOEFFS_B, COMPRESSFROM_B, COMPRESSTO_B, public_B);
  compress(NCOEFFS_C, COMPRESSFROM_C, COMPRESSTO_C, public_C);

  // new random masked ciphertext
  get_rand_share(NCOEFFS_B, Q, public_B_rand);
  get_rand_share(NCOEFFS_C, P, public_C_rand);

  // NOTE: reuse B for B_rand and C for C_rand to keep
  // memeory addresses the same for both cases
  mask(NSHARES, NCOEFFS_B, B, public_B_rand);
  mask(NSHARES, NCOEFFS_C, C, public_C_rand);

  unmask_ct(g_ct1, B, C);
  compress(NCOEFFS_B, COMPRESSFROM_B, COMPRESSTO_B, g_ct1);
  compress(NCOEFFS_C, COMPRESSFROM_C, COMPRESSTO_C, g_ct1 + NCOEFFS_B);

  concat_ct(g_ct2, public_B, public_C);
  diff_ct(g_ct_diff, g_ct1, g_ct2, NCOEFFS_B + NCOEFFS_C);  
  print_ct(g_ct_diff, NCOEFFS_B + NCOEFFS_C);
}

static inline void
run_masked_cmp_gf() {
  uint32_t rc = MaskedComparison_GF(B, C, public_B, public_C);
  assert(rc == 0); // should fail
}

static inline void
dump_BC_Bitsliced(FILE* f, int32_t trace_num) {

    uint32_t size_ss = sizeof(uint32_t) * (SIMPLECOMPBITS) * NSHARES;
    fwrite(&size_ss, sizeof(uint32_t), 1, f); 
    
    // perform the same operations as ChipWhisperer
    for(int32_t i = 0; i < trace_num; ++i) {
        prep_ct_no_dec();
        run_masked_cmp_gf();
        fwrite(BC_Bitsliced, sizeof(uint32_t), (SIMPLECOMPBITS) * NSHARES, f);

        prep_ct_dec();
        run_masked_cmp_gf();
        fwrite(BC_Bitsliced, sizeof(uint32_t), (SIMPLECOMPBITS) * NSHARES, f);
    }
}

// enable TRNG/PRNG
uint8_t en_rand = 1;

int main(void)
{

    hal_setup();
    const char* ofname = "BC_Bitsliced-dump.bin";
    const int32_t trace_num = 5000;

    FILE* ofile = fopen(ofname, "wb");
    if(!ofile) {
        fprintf(stderr, "fail to open file: %s\n", ofname);
        return 1;
    }

    public_B = malloc(sizeof(uint32_t) * NCOEFFS_B);
    public_C = malloc(sizeof(uint32_t) * NCOEFFS_C);
    B = malloc(sizeof(uint32_t) * NCOEFFS_B * NSHARES);
    C = malloc(sizeof(uint32_t) * NCOEFFS_C * NSHARES);
    BC_Bitsliced = malloc(sizeof(uint32_t*) * (SIMPLECOMPBITS) * NSHARES);

    public_B_rand = malloc(sizeof(uint32_t) * NCOEFFS_B);
    public_C_rand = malloc(sizeof(uint32_t) * NCOEFFS_C);

    g_ct1 = malloc(sizeof(uint32_t) * (NCOEFFS_B + NCOEFFS_C));
    g_ct2 = malloc(sizeof(uint32_t) * (NCOEFFS_B + NCOEFFS_C));
    g_ct_diff = malloc(sizeof(uint32_t) * (NCOEFFS_B + NCOEFFS_C));

    if(!public_B || !public_C || !B || !C || !BC_Bitsliced ||
       !public_B_rand || !public_C_rand || !g_ct1 || !g_ct2 || !g_ct_diff) {
        fprintf(stderr, "fail to allocate memory\n");
        goto cleanup;
    }

    dump_BC_Bitsliced(ofile, trace_num);

cleanup:
    free(public_B);
    free(public_C);
    free(B);
    free(C);
    free(BC_Bitsliced);
    free(public_B_rand);
    free(public_C_rand);

    free(g_ct1);
    free(g_ct2);
    free(g_ct_diff);

    fclose(ofile);

    return 0;
}
