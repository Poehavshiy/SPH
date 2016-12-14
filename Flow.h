//
// Created by nikita on 13.05.16.
//

#ifndef SPHSM6_FLOW_H
#define SPHSM6_FLOW_H

#include <sstream>
#include <fstream>
#include "Calculator.h"

extern Particle* for_debugin;
extern QGraphicsScene *scene_debug;
#define PI 3.14159265

class Flow {
    bool sphere ;

    double maxP, maxp, maxe;

    double smooth_length;

    int N;
    //
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
    void set_bound(const string &boundaryFile, int& x, int& y);
    //
    void set_init(const string &initFile);
    //
    void build_bound_part();
    //
    void build_init_part(vector<Point> &branch, vector<double> &cond, int step,bool shape);
public:
    friend class Flow_Drawer;
    //
    friend class Visualisator;
    //
    Flow(const string &boundaryFile, const string &initFile);
    //
    const void* get_part(int type){
        if( type == 0) return &boundaries;
        else if(type == 1) return &data;
        else if(type == 2) return s_distribution->get_parsing();
    }

};


#endif //SPHSM6_FLOW_H
