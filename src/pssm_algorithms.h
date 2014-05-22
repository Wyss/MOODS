// pssm_algorithms - collection of algorithms for finding PSSM matches from sequences
// Copyright (C) 2007-2009 Pasi Rastas, Janne Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.

#ifndef _PSSM_ALG_H
#define _PSSM_ALG_H

#include <stdint.h>

// for data about position in a sequence
typedef long position_t;

// for bit parallel magic
typedef uint_fast32_t bits_t;

typedef int bool;

// matrix scores
typedef double score_t;
const score_t SCORE_MIN = DBL_MIN;
const score_t PVAL_DP_MULTIPLIER = 1000.0;


// Struct for storing data about matrix matches
typedef struct {
    position_t position;
    score_t score;
} match_data_t;

int expectedDifferences(const score_t **mat, int n, int m, const double* bg, double**ret);

// Output list element for mm AC automaton
struct {
    score_t score;
    int matrix;
    bool full;
} OutputListElementMulti;

void multipleMatrixLookaheadFiltrationDNASetup(const int q,  
    const score_t ***matrices, int *matrices_dims,
    OutputListElementMulti **output, 
    int *window_positions, int *m, int** orders,
    score_t **L,
    const double *bg, const score_t *thresholds);

match_data_t * doScan(const unsigned char *s, 
	const int q, const score_t ***matrices, 
    OutputListElementMulti **output, 
	const int *window_positions, const int *m, int** orders, 
	const score_t *L,
	const score_t *thresholds);
#endif
