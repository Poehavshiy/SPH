//
// Created by nikita on 13.05.16.
//

#ifndef SPHSM6_READWRITE_H
#define SPHSM6_READWRITE_H

#include "Particile.h"
class SpaceParsing;

namespace ReadWrite {
    void parse_bound(string &current, vector<vector<double>> &geometry);

    //
    Point parse_init(const string &current);

    //
    void data_write(ostream &os, vector<vector<Particle>> &boundaries);

    //
    void data_write(ostream &os, vector<Particle>& data);

}

namespace Help {
    void fill_bound_line(vector<Particle> &current, vector<double> &line);
}
#endif //SPHSM6_READWRITE_H
