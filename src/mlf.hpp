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
    std::vector<std::vector< OutputListElementMulti> > output; 
    intArray window_positions;
    intArray m; 
    intMatrix orders; 
    scoreMatrix L;
    scoreArray thresholds;
    scoreArray bg;

    MOODS_MLF(void) {
        // q = 0;
        // matrices = new std::vector<scoreMatrix>();
        // output = new std::vector<std::vector< OutputListElementMulti> >();
        // window_positions = new intArray();
        // m = new intArray();
        // orders = new intMatrix();
        // L = new scoreMatrix();
        // thresholds = new scoreArray();
        // bg = new scoreArray();
    } 

    int multipleMatrixLookaheadFiltrationDNASetup();
    match_data_t doScan(MOODS_MLF &in, int *rc);
};

typedef struct MOODS_MLF MOODS_MLF;

#endif
