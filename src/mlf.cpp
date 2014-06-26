
int MOODS_MLF::multipleMatrixLookaheadFiltrationDNASetup() {

    const int BITSHIFT = 2;
    const unsigned int numA = 4; // 2**BITSIFT

    m.resize(matrices.size(), 0);
    for (int i = 0; i < (int) matrices.size(); ++i) {
        m[i] = matrices[i][0].size();
    }

    // Calculate entropies for all matrices
    vector<doubleArray> goodnesses;
    goodnesses.reserve(matrices.size());

    for (int i = 0; i < (int)matrices.size(); ++i) {
        goodnesses.push_back(expectedDifferences(matrices[i], bg));
    }

    window_positions.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k) {
        if (q >= m[k]) {
            window_positions.push_back(0);
        } else {
            double current_goodness = 0;
            for (int i = 0; i < q; ++i) {
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
            window_positions.push_back(window_pos);
        }
    }

    // Calculate lookahead scores for all matrices
    scoreMatrix T;
    T.reserve(matrices.size());

    for (int k = 0; k < (int)matrices.size(); ++k) {
        scoreArray C(m[k],0);
        for (int j = m[k] - 1; j > 0; --j) {
            score_t max = SCORE_MIN;
            for (unsigned int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            C[j - 1] = C[j] + max;
        }
        T.push_back(C);
    }

    // Pre-window scores
    scoreArray P;
    P.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k) {
        score_t B = 0;
        for (int j = 0; j < window_positions[k]; ++j) {
            score_t max = SCORE_MIN;
            for (unsigned int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            B += max;
        }
        P.push_back(B);
    }

    // Arrange matrix indeces not in window by entropy, for use in scanning
    orders.reserve(matrices.size());
    L.reserve(matrices.size());

    for (unsigned short k = 0; k < (int) matrices.size(); ++k) {
        if (q >= m[k]) {
            intArray temp1;
            orders.push_back(temp1);
            scoreArray temp2;
            L.push_back(temp2);
        }
        else {
            intArray order(m[k]-q, 0);
            for (int i = 0; i < window_positions[k]; ++i) {
                order[i] = i;
            }
            for (int i = window_positions[k]+q; i < m[k]; ++i) {
                order[i-q] = i;
            }

            compareRows comp;
            comp.goodness = &(goodnesses[k]);

            sort(order.begin(), order.end(), comp);

            orders.push_back(order);

            scoreArray K(m[k]-q, 0);
            for (int j = m[k]-q-1; j > 0; --j) {
                score_t max = INT_MIN;
                for (unsigned int i = 0; i < numA; ++i) {
                    if (max < matrices[k][i][order[j]]) {
                        max = matrices[k][i][order[j]];
                    }
                }
                K[j - 1] = K[j] + max;
            }
            L.push_back(K);
        }
    }

    // const bits_t size = 1 << (BITSHIFT * q); // numA^q
    // const bits_t BITAND = size - 1;
    // vector<vector< OutputListElementMulti> > output(size);

    {
        bitArray sA(q,0);
        while (true) {
            bits_t code = 0;
            for (int j = 0; j < q; ++j) {
                code = (code << BITSHIFT) | sA[j];
            }

            for (unsigned int k = 0; k < matrices.size(); ++k ) {
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
                        output[code].push_back(temp);
                    }
                } else {
                    score_t score = 0;
                    for (int i = 0; i < q; ++i) {
                        score += matrices[k][sA[i]][i + window_positions[k]];
                    }
                    if (score + P[k] + T[k][q + window_positions[k]-1] >= tol[k]) {
                        OutputListElementMulti temp;
                        temp.full = false;
                        temp.matrix = k;
                        temp.score = score;
                        output[code].push_back(temp);
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
    }

    // return doScan(q, matrices, output, window_positions, m, orders, L, tol);
}

vector<matchArray> MOODS_MLF::doScan(const charArray &s) {
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
                    tolerance = tol[k];
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