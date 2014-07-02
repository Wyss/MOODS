// pssm_algorithms - collection of algorithms for finding PSSM matches from sequences
// Copyright (C) 2007-2009 Pasi Rastas, Janne Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.

#ifndef _MLF_H
#define _MLF_H

#include "pssm_algorithms.hpp"

struct MOODS_MLF {
    int q; 
    std::vector<scoreMatrix> matrices;
    intArray window_positions;
    intArray m; 
    intMatrix orders; 
    scoreMatrix L;
    scoreArray thresholds;
    scoreArray bg;
    
    static const int BITSHIFT = 2;
    static const unsigned int numA = 4; // 2**BITSIFT

    bits_t size; // numA^q
    std::vector<std::vector< OutputListElementMulti> > output; 
    bits_t BITAND;

    // matrices could be required by the constructor, or I can assign it 
    // and look it up later.
    MOODS_MLF(int q_in) : q(q_in),  size(1 << (BITSHIFT * q)), output(size), BITAND(size -1) {
    }

    // must assign member matrices before calling the below methods
    int multipleMatrixLookaheadFiltrationDNASetup();
    std::vector<matchArray> doScan(const charArray &s, int *rc);
};

// for C ease of use
typedef struct MOODS_MLF MOODS_MLF;

#endif
