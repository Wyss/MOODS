
#include "pssm_algoritm.h"
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
    int i, k, rc;
    const int msize = (int) kv_size(matrices);
    for (i = 0; i < msize; ++i) {
        m[i] = (int) kv_size(kv_A(kv_A(matrices, i), 0));
    }

    // Calculate entropies for all matrices
    double **goodnesses;
    goodnesses = (double **) malloc(msize*sizeof(double *));
    for (i = 0; i < msize; ++i) {
        goodnesses[i] = expectedDifferences(kv_A(matrices, i), bg, &rc);
        if (rc < 0) {
            goto mlf_fail;
        }
    }

    // intArray window_positions;
    // window_positions.reserve(matrices.size());
    for (k = 0; k < msize; ++k) {
        if (q >= kv_A(m, k)) {
            window_positions[0] = 0;
        } else {
            double current_goodness = 0;
            for (i = 0; i < q; ++i) {
                current_goodness += goodnesses[k][i];
            }

            double max_goodness = current_goodness;
            int window_pos = 0;

            for (int i = 0; i < m[k] - q; ++i) {
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
    score_t *C = NULL;
    for (int k = 0; k < msize; ++k) {
        C = (score_t *) calloc(m[k], sizeof(score_t));
        if (C == NULL) {
            goto mlf_fail;
        }
        for (int j = m[k] - 1; j > 0; --j) {
            score_t max = SCORE_MIN;
            for (unsigned int i = 0; i < numA; ++i) {
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

    for (int k = 0; k < msize; ++k) {
        score_t B = 0;
        for (int j = 0; j < window_positions[k]; ++j) {
            score_t max = SCORE_MIN;
            for (unsigned int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j]) {
                    max = matrices[k][i][j];
                }
            }
            B += max;
        }
        P[k] = B;
    }

    // Arrange matrix indeces not in window by entropy, for use in scanning
    // intMatrix orders;
    // orders.reserve(matrices.size());
    // scoreMatrix L;
    // L.reserve(matrices.size());

    for (unsigned short k = 0; k < msize; ++k) {
        if (q >= m[k]) {
            // orders.push_back(temp1);
            orders[k] = NULL;
            // L.push_back(temp2);
            L[k] = NULL;
        } else {
            // intArray order(m[k]-q, 0);
            int *order = NULL;
            order = (int *) calloc(m[k]-q, sizeof(int));
            if (order == NULL) {
                goto mlf_fail;
            }

            for (int i = 0; i < window_positions[k]; ++i) {
                order[i] = i;
            }
            for (int i = window_positions[k]+q; i < m[k]; ++i) {
                order[i-q] = i;
            }

            compareRows comp;
            comp.goodness = &(goodnesses[k]);

            sort(order.begin(), order.end(), comp);

            // orders.push_back(order);
            orders[k] = order;


            score_t * K = NULL;
            K = (score_t *) calloc(m[k] - q, sizeof(score_t));
            if (K == NULL) {
                goto mlf_fail;
            }

            for (int j = m[k]-q-1; j > 0; --j) {
                score_t max = INT_MIN;
                for (unsigned int i = 0; i < numA; ++i) {
                    if (max < matrices[k][i][order[j]]) {
                        max = matrices[k][i][order[j]];
                    }
                }
                K[j - 1] = K[j] + max;
            }
            // L.push_back(K);
            L[k] = K;
        }
    }

    // const bits_t size = 1 << (BITSHIFT * q); // numA^q
    // const bits_t BITAND = size - 1;
    // vector<vector< OutputListElementMulti> > output(size);

    bit_t * sA = (bit_t *) calloc(q, sizeof(bit_t));
    OutputListElementMulti *temp_ptr;
    while (true) {
        bits_t code = 0;
        for (int j = 0; j < q; ++j) {
            code = (code << BITSHIFT) | sA[j];
        }

        for (unsigned int k = 0; k < msize; ++k ) {
            if (m[k] <= q) {
                score_t score = 0;
                for (int i = 0; i < m[k]; ++i) {
                    score += matrices[k][sA[i]][i];
                }
                if (score >= tol[k]) {
                    OutputListElementMulti temp;
                    temp.full = true;
                    temp.matrix = k;
                    temp.score = score;
                    kv_push_plus(OutputListElementMulti, output[code], temp, temp_ptr, mlf_fail);
                }
            } else {
                score_t score = 0;
                for (int i = 0; i < q; ++i) {
                    score += matrices[k][sA[i]][i + window_positions[k]];
                }
                if (score + P[k] + T[k][q + window_positions[k]-1] >= thresholds[k]) {
                    OutputListElementMulti temp;
                    temp.full = false;
                    temp.matrix = k;
                    temp.score = score;
                    kv_push_plus(OutputListElementMulti, output[code], temp, temp_ptr, mlf_fail);
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

}

vector<matchArray> doScan( const charArray &s, 
    const int q, const vector<scoreMatrix> &matrices,
    vector<vector< OutputListElementMulti> > &output, 
    const intArray &window_positions, const intArray &m, intMatrix &orders,
    const scoreMatrix &L,
    const scoreArray &thresholds)
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