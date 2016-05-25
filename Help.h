//
// Created by nikita on 13.05.16.
//

#ifndef SPHSM6_READWRITE_H
#define SPHSM6_READWRITE_H
#include "Particile.h"

namespace ReadWrite{
    void parseBound( string& current, vector<vector<double>>& geometry);
    //
    Point parseInit(const string& current);
}

namespace Help {
    void fillBoundLine(vector<Particle>& current, vector<double>& line);
}
#endif //SPHSM6_READWRITE_H
