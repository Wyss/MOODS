
#include "pssm_algoritm.h"
#include "ksort.h"

#define cmp_gt_double(a, b) ((a) > (b))
KSORT_INIT(double, double, cmp_gt_double)


double * expectedDifferences(const score_matrix_t smatrix, const double *bg, int *rc) {
    int numA =  (int) kv_size(smatrix);
    int m = (int) kv_size( kv_A(smatrix, 0) );
    int i,j;
    score_t max;
    score_t temp;
    double * ed = NULL;
    ed = (double *) malloc(m*sizeof(double));
    if (ed == NULL) {
        goto ed_fail;
    }
    
    for (i = 0; i < m; ++i) {
        max = SCORE_MIN;
        for (j = 0; j < numA; ++j) {
            temp = kv_A(kv_A(smatrix, j), i);
            if (max <  temp) {
                max = temp;
            }
        }

        ed[i] = max;

        for (j = 0; j < numA; ++j) {
            ed[i] -= bg[j] * kv_A(kv_A(smatrix, j), i);
        }
    }
    *rc = 0;
    return ed;
    ed_fail:
        *rc = 1;
        return ed;
}

void multipleMatrixLookaheadFiltrationDNASetup(const int q,  
    const score_matrix_vec_t matrices,
    OutputListElementMulti **output, 
    int *window_positions, int *m, int** orders,
    score_t **L,
    const double *bg, const score_t *thresholds) {

    const int BITSHIFT = 2;
    const unsigned int numA = 4; // 2**BITSIFT

    // intArray m(matrices.size(), 0);
    int i, j, k, rc;

    const int msize = (int) kv_size(matrices);
    for (i = 0; i < msize; ++i) {
        m[i] = (int) kv_size(kv_A(kv_A(matrices, i), 0));
    }

    // Calculate entropies for all matrices
    double **goodnesses = NULL;
    goodnesses = (double **) malloc(msize*sizeof(double *));
    if (goodnesses == NULL) {
        goto mlf_fail;
    }
    for (i = 0; i < msize; ++i) {
        goodnesses[i] = expectedDifferences(kv_A(matrices, i), bg, &rc);
        if (rc < 0) {
            goto mlf_fail;
        }
    }

    int *window_positions = NULL;
    window_positions = (int *) malloc(msize*sizeof(int));
    if (window_positions == NULL) {
        goto mlf_fail;
    }
    double current_goodness;
    int window_pos;
    double max_goodness;
    for (k = 0; k < msize; ++k) {
        if (q >= kv_A(m, k)) {
            window_positions[0] = 0;
        } else {
            current_goodness = 0;
            for (i = 0; i < q; ++i) {
                current_goodness += goodnesses[k][i];
            }

            max_goodness = current_goodness;
            window_pos = 0;
            const int il = m[k] - q;
            for (i = 0 i < il; ++i) {
                current_goodness -= goodnesses[k][i];
                current_goodness += goodnesses[k][i+q];
                if (current_goodness > max_goodness) {
                    max_goodness = current_goodness;
                    window_pos = i+1;
                }
            }
            window_positions[k] = window_pos;
        }
    }

    // Calculate lookahead scores for all matrices
    score_t **T = NULL;
    T = (score_t **) malloc(msize*sizeof(score_t*));
    if (T == NULL) {
        goto mlf_fail;
    }

    score_t *C = NULL;
    for (k = 0; k < msize; ++k) {
        C = (score_t *) calloc(m[k], sizeof(score_t));
        if (C == NULL) {
            goto mlf_fail;
        }
        for (j = m[k] - 1; j > 0; --j) {
            score_t max = SCORE_MIN;
            for (i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            C[j - 1] = C[j] + max;
        }
        T[k] = C;
    }

    // Pre-window scores
    score_t *P = NULL;
    P = (score_t *) malloc(msize*sizeof(score_t));
    if (P == NULL) {
        goto mlf_fail;
    }
    for (k = 0; k < msize; ++k) {
        score_t B = 0;
        for (j = 0; j < window_positions[k]; ++j) {
            score_t max = SCORE_MIN;
            for (i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j]) {
                    max = matrices[k][i][j];
                }
            }
            B += max;
        }
        P[k] = B;
    }

    // Arrange matrix indeces not in window by entropy, for use in scanning
    int** orders = NULL;
    orders = (int **) malloc(msize*sizeof(int *));
    if (orders == NULL) {
        goto mlf_fail;
    }
    score_t **L = NULL;
    L = (score_t **) malloc(msize*sizeof(score_t *));
    if (L == NULL) {
        goto mlf_fail;
    }
    int *order;
    for (k = 0; k < msize; ++k) {
        if (q >= m[k]) {
            orders[k] = NULL;
            L[k] = NULL;
        } else {
            // intArray order(m[k]-q, 0);
            *order = NULL;
            order = (int *) calloc(m[k]-q, sizeof(int));
            if (order == NULL) {
                goto mlf_fail;
            }

            for (i = 0, il = window_positions[k]; i < il; ++i) {
                order[i] = i;
            }
            for (i = window_positions[k]+q, il=m[k]; i < il; ++i) {
                order[i - q] = i;
            }

            // compareRows comp;
            // comp.goodness = &(goodnesses[k]);

            // sort(order.begin(), order.end(), comp);
            ks_mergesort_double(m[k]-q, order, 0);

            orders[k] = order;

            score_t * K = NULL;
            K = (score_t *) calloc(m[k] - q, sizeof(score_t));
            if (K == NULL) {
                goto mlf_fail;
            }

            for (j = m[k]-q-1; j > 0; --j) {
                score_t max = INT_MIN;
                for (i = 0; i < numA; ++i) {
                    if (max < matrices[k][i][order[j]]) {
                        max = matrices[k][i][order[j]];
                    }
                }
                K[j - 1] = K[j] + max;
            }
            L[k] = K;
        }
    }

    const bits_t size = 1 << (BITSHIFT * q); // numA^q
    // const bits_t BITAND = size - 1;
    // vector<vector< OutputListElementMulti> > output(size);

    OutputListElementMulti_vec_t* output = NULL;
    output = (OutputListElementMulti_vec_t*) calloc(size, sizeof(OutputListElementMulti_vec_t));
    if (output == NULL) {
        goto mlf_fail;
    }
    bit_t * sA = (bit_t *) calloc(q, sizeof(bit_t));
    if (sA == NULL) {
        goto mlf_fail;
    }
    OutputListElementMulti_t *temp_ptr;
    OutputListElementMulti_t temp;
    while (1) {
        bits_t code = 0;
        for (int j = 0; j < q; ++j) {
            code = (code << BITSHIFT) | sA[j];
        }

        for (k = 0; k < msize; ++k ) {
            if (m[k] <= q) {
                score_t score = 0;
                const int il = m[k];
                for (i = 0; i < il; ++i) {
                    score += matrices[k][sA[i]][i];
                }
                if (score >= tol[k]) {
                    temp.full = 1;
                    temp.matrix = k;
                    temp.score = score;
                    // push plus takes care of memory allocation
                    kv_push_safe(OutputListElementMulti_vec_t, output[code], temp, temp_ptr, mlf_fail);
                }
            } else {
                score_t score = 0;
                for (i = 0; i < q; ++i) {
                    score += matrices[k][sA[i]][i + window_positions[k]];
                }
                if (score + P[k] + T[k][q + window_positions[k]-1] >= thresholds[k]) {
                    temp.full = 0;
                    temp.matrix = k;
                    temp.score = score;
                    kv_push_safe(OutputListElementMulti_vec_t, output[code], temp, temp_ptr, mlf_fail);
                }
            }
        }

        int pos = 0;
        while (pos < q) {
            if (sA[pos] < numA - 1) {
                ++sA[pos];
                break;
            } else {
                sA[pos] = 0;
                ++pos;
            }
        }

        if (pos == q) {
            break;
        }
    }

    mlf_fail:
        // do stuff

}

match_data_t * doScan(const unsigned char *s, 
    const int q, const score_matrix_vec_t *matrices, 
    OutputListElementMulti **output, 
    const int_vec_t *window_positions, const int_vec_t *m, int_matrix_t *orders, 
    const score_matrix_t *L,
    const score_vec_t *thresholds)
{
    const int BITSHIFT = 2;
    const bits_t size = 1 << (BITSHIFT * q); // numA^q
    const bits_t BITAND = size - 1;
    
    const position_t n = s.size();
    // Scanning

    vector<matchArray> ret;
    for (unsigned int i = 0; i < matrices.size(); ++i) {
        matchArray temp;
        ret.push_back(temp);
    }

    matchData hit;

    score_t score;
    position_t k;
    position_t limit;
    position_t ii;
    score_t tolerance;
    intArray::iterator z;

    bits_t code = 0;
    for (position_t ii = 0; ii < q - 1; ++ii) {
        code = (code << BITSHIFT) + s[ii];
    }

    for (position_t i = 0; i < n - q + 1; ++i) {
        code = ((code << BITSHIFT) + s[i + q - 1]) & BITAND;

        if (!output[code].empty()) {
            for (vector<OutputListElementMulti>::iterator y = output[code].begin(); y != output[code].end(); ++y) {
                if (y->full) { // A Hit for a matrix of length <= q
                    hit.position = i;
                    hit.score = y->score;
                    ret[y->matrix].push_back(hit);
                    continue;
                }
                // A possible hit for a longer matrix. Don't check if matrix can't be positioned entirely on the sequence here
                if (i - window_positions[y->matrix] >= 0 && i + m[y->matrix] - window_positions[y->matrix] <= n) { 
                    score = y->score;
                    k = y->matrix;
                    limit = m[k] - q;
                    ii = i - window_positions[k];
                    tolerance = thresholds[k];
                    z = orders[k].begin();
                    for (int j = 0; j < limit  ;++j) {
                        score += matrices[k][s[ii+(*z)]][*z];
                        if (score + L[k][j] < tolerance) {
                            break;
                        }
                        ++z;
                    }
                    if (score >= tolerance) {
                        hit.position = i - window_positions[k];
                        hit.score = score;
                        ret[k].push_back(hit);
                    }
                }
            }
        }
    }

    for (position_t i = n - q + 1; i < n; ++i) { // possible hits for matrices shorter than q near the end of sequence
        code = (code << BITSHIFT) & BITAND; // dummy character to the end of code

        if (!output[code].empty()) {

            for (vector<OutputListElementMulti>::iterator y = output[code].begin(); y != output[code].end(); ++y) {
                if (y->full && m[y->matrix] < n - i + 1) { // only sufficiently short hits are considered
                    hit.position = i;
                    hit.score = y->score;
                    ret[y->matrix].push_back(hit);
                }
            }
        }
    }

    return ret;
}