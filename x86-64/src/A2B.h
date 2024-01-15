/*
 * MIT License
 *
 * Copyright (c) 2021-2022: imec-COSIC KU Leuven, 3001 Leuven, Belgium
 * Authors: Michiel Van Beirendonck <michiel.vanbeirendonck@esat.kuleuven.be>
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

#ifndef A2B_H
#define A2B_H

#include <stdint.h>
#include <stddef.h>
#include "params.h"

#ifdef DEBUG
#include <stdio.h>
#include <assert.h>
#endif

void A2B(size_t nshares, uint64_t B[nshares], const uint64_t A[nshares]);
void A2B32(size_t nshares, uint32_t B[nshares], const uint32_t A[nshares]);
void A2B_bitsliced(size_t nshares, size_t nbits, uint32_t B[32][nshares], const uint32_t A[32][nshares]);
void A2B_keepbitsliced(size_t nshares, uint32_t ncoefsb, uint32_t compressfrom_b, uint32_t compressto_b, uint32_t ncoefsc, uint32_t compressfrom_c, uint32_t compressto_c, uint32_t out[SIMPLECOMPBITS][NSHARES], const uint32_t Bp[ncoefsb][NSHARES], const uint32_t Cp[ncoefsc][NSHARES]);

#endif // A2B_H