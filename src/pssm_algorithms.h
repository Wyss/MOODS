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
#include "kvec.h"

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

// adds a chance for error handling
#define kv_resize_plus(type, v, s, temp, label)  \
(v).m = (s); \
temp = (type *)realloc((v).a, sizeof(type) * (v).m); \
if (temp != NULL) { \
    (v).a = temp; \
} else { \
    free((v).a); \
    goto label; \
}

// adds a chance for error handling
#define kv_push_plus(type, v, x, temp, label) do {                  \
        if ((v).n == (v).m) {                                       \
            (v).m = (v).m ? (v).m << 1 : 2;                         \
            temp = (type*)realloc((v).a, sizeof(type) * (v).m);     \
            if (temp == NULL) {                                     \
                free((v).a);                                        \
                (v).a = NULL;                                       \
                goto label;                                         \
            }                                                       \
            (v).a = temp;                                           \
        }                                                           \
        (v).a[(v).n++] = (x);                                       \
    } while (0)


// define vectors
typedef kvec_t(int) int_vec_t;
kv_resize_init(double_vec_t, double);

typedef kvec_t(int_vec_t) int_matrix_t;
typedef kvec_t(score_t) score_vec_t;
typedef kvec_t(score_vec_t) score_matrix_t;
typedef kvec_t(score_matrix_t) score_matrix_vec_t;
typedef kvec_t(double) double_vec_t;

int expectedDifferences(const score_matrix_t *mat, const double *bg, double **ret);

// Output list element for mm AC automaton
struct {
    score_t score;
    int matrix;
    bool full;
} OutputListElementMulti;

void multipleMatrixLookaheadFiltrationDNASetup(const int q,  
    const score_matrix_vec_t *matrices,
    OutputListElementMulti **output, 
    int_vec_t *window_positions, int_vec_t *m, int_matrix_t *orders,
    score_matrix_t *L,
    const double *bg, const score_vec_t *thresholds);

match_data_t * doScan(const unsigned char *s, 
	const int q, const score_matrix_vec_t *matrices, 
    OutputListElementMulti **output, 
	const int_vec_t *window_positions, const int_vec_t *m, int_matrix_t *orders, 
	const score_matrix_t *L,
	const score_vec_t *thresholds);
#endif
