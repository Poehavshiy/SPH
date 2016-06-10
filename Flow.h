//
// Created by nikita on 13.05.16.
//

#ifndef SPHSM6_FLOW_H
#define SPHSM6_FLOW_H

#include <sstream>
#include <fstream>
#include "Calculator.h"


class Flow {
    SpaceParsing *s_distribution;
    //
    Calculator* calculator;
    //
    vector<Particle> data;
    //
    vector<vector<Particle>> boundaries;//each vector represents line, consisted from Particles
    //1- x1, 2-y1, 3-x2, 4-y2, 5-number of 1st type particles
    vector<vector<double>> geometry;
    //
    void set_bound(const string &boundaryFile);
    //
    void set_init(const string &initFile);
    //
    void build_bound_part();
    //
    void build_init_part(vector<Point> &branch);
    //
    void calculate_step();
public:
<<<<<<< HEAD
    friend class Flow_Drawer;
    //
    Flow(const string &boundaryFile, const string &initFile);
    //
    virtual void calculate();
    //
    const void* get_part(int type){
        if( type == 0) return &boundaries;
        else if(type == 1) return &data;
        else if(type == 2) return s_distribution->get_parsing();
    }

=======
    //
    Flow(const string &boundaryFile, const string &initFile);
    //
    void calculate();
    //
>>>>>>> 883cc5708c435c9115b864d6fdf26606dffc1703
};


#endif //SPHSM6_FLOW_H
