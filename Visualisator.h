//
// Created by nikita on 12.09.16.
//

#ifndef SPHSM6_VISUALISATOR_H
#define SPHSM6_VISUALISATOR_H

#include "Flow.h"

class Visualisator {
    Flow *target;
    string path;
    string name = "/test.csv.";
    int counter;
    std::ofstream outfile;
public:
    Visualisator(Flow *to_vis, string dir);

    void write_step();
};


#endif //SPHSM6_VISUALISATOR_H
