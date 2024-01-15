/*
    This file is part of the ChipWhisperer Example Targets
    Copyright (C) 2012-2017 NewAE Technology Inc.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
#include <stddef.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // memset

#include "MaskedComparison.h"
#include "randombytes.h"
#include "params.h"
#include "profile_util.h" // the original common/hal.h

// ignore hal_send_str
#define hal_send_str(str)

// Enable/disable randomness for masks ON/OFF
uint8_t en_rand = 1;

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
    const size_t shift_bit_num = compressfrom - compressto;
    for (size_t i = 0; i < ncoeffs; i++)
    {
        submitted_poly[i] >>= shift_bit_num;
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

uint32_t* public_B;
uint32_t* public_C;
uint32_t* B;
uint32_t* C;

uint32_t* public_B_rand;
uint32_t* public_C_rand;

uint8_t reset(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t* buf) {
  memset(public_B, 0x0, sizeof(uint32_t) * NCOEFFS_B);
  memset(public_C, 0x0, sizeof(uint32_t) * NCOEFFS_C);
  memset(B, 0x0, sizeof(uint32_t) * NCOEFFS_B * NSHARES);
  memset(C, 0x0, sizeof(uint32_t) * NCOEFFS_C * NSHARES);
  memset(public_B_rand, 0x0, sizeof(uint32_t) * NCOEFFS_B);
  memset(public_C_rand, 0x0, sizeof(uint32_t) * NCOEFFS_C);
  return 0x0;
}

uint8_t prep_ct_no_dec(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t* buf) {
  // new ciphertext
  get_rand_share(NCOEFFS_B, Q, public_B);
  get_rand_share(NCOEFFS_C, P, public_C);

  // create reencrypted poly as a (correct) sharing of public poly
  mask(NSHARES, NCOEFFS_B, B, public_B);
  mask(NSHARES, NCOEFFS_C, C, public_C);

  // Now compress public poly
  compress(NCOEFFS_B, COMPRESSFROM_B, COMPRESSTO_B, public_B);
  compress(NCOEFFS_C, COMPRESSFROM_C, COMPRESSTO_C, public_C);

  // tamper the ciphertext after compression
  public_C[0] += 1;

  return 0x00;
}

uint8_t prep_ct_dec(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t* buf) {
  // first steps are the same as prep_ct_no_dec
  // new ciphertext
  get_rand_share(NCOEFFS_B, Q, public_B);
  get_rand_share(NCOEFFS_C, P, public_C);

  // NOTE: no need to compute this. B and C are overwritten later
  // create reencrypted poly as a (correct) sharing of public poly
  //mask(NSHARES, NCOEFFS_B, B, public_B);
  //mask(NSHARES, NCOEFFS_C, C, public_C);

  // Now compress public poly
  compress(NCOEFFS_B, COMPRESSFROM_B, COMPRESSTO_B, public_B);
  compress(NCOEFFS_C, COMPRESSFROM_C, COMPRESSTO_C, public_C);

  // tamper the ciphertext after compression
  public_C[0] += 1;
  
  // new random masked ciphertext
  get_rand_share(NCOEFFS_B, Q, public_B_rand);
  get_rand_share(NCOEFFS_C, P, public_C_rand);

  // NOTE: reuse B for B_rand and C for C_rand to keep
  // memeory addresses the same for both cases
  mask(NSHARES, NCOEFFS_B, B, public_B_rand); // B_rand (reuse B)
  mask(NSHARES, NCOEFFS_C, C, public_C_rand); // C_rand (reuse C)
  
  return 0x00;
}

uint8_t run_masked_cmp_gf(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t* buf) {
  uint64_t rc = MaskedComparison_GF(B, C, public_B, public_C);
  assert(rc == 0); // should fail
  return 0x00;
}

int main(void) {
    platform_init();
    init_uart();
    trigger_setup();
    simpleserial_init();

    // https://chipwhisperer.readthedocs.io/en/latest/simpleserial.html
    // we use SimpleSerial v2.1
    #if (SS_VER != SS_VER_2_1)
    // TODO: is there a better way to print this?
    char const * const ss_error_msg = "Wrong SimpleSerial version. Need 2.1\n";
    for(int i = 0; i < strlen(ss_error_msg); ++i)
        putch(ss_error_msg[i]);

    // TODO: proper way to terminate?
    return 0;
    #endif

    // allocate buffers
    public_B = malloc(sizeof(uint32_t) * NCOEFFS_B);
    public_C = malloc(sizeof(uint32_t) * NCOEFFS_C);
    B = malloc(sizeof(uint32_t) * NCOEFFS_B * NSHARES);
    C = malloc(sizeof(uint32_t) * NCOEFFS_C * NSHARES);

    public_B_rand = malloc(sizeof(uint32_t) * NCOEFFS_B);
    public_C_rand = malloc(sizeof(uint32_t) * NCOEFFS_C);

    if(!public_B || !public_C || !B || !C
       || !public_B_rand || !public_C_rand 
       ) {
	    putch('m');
	    putch('\n');
	    return 0;
    }

    
    if(simpleserial_addcmd('c', 0, run_masked_cmp_gf)
       || simpleserial_addcmd('x', 0, reset)
       || simpleserial_addcmd('n', 0, prep_ct_no_dec)
       || simpleserial_addcmd('f', 0, prep_ct_dec)
       ) {
      // fail to add all the listeners
      putch('f');
      putch('\n');
      return 2;
    }
    
    while(1)
        simpleserial_get();

    free(public_B);
    free(public_C);
    free(B);
    free(C);
    free(public_B_rand);
    free(public_C_rand);

    return 0;
}
