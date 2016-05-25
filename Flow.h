//
// Created by nikita on 13.05.16.
//

#ifndef SPHSM6_FLOW_H
#define SPHSM6_FLOW_H

#include "Particile.h"
#include <sstream>
#include <fstream>

struct BoundaryPart {
    vector<Particle> particiles;
    double r;
    BoundaryPart(vector<Particle>& income, double R){
        particiles=income;
        r=R;
    }
};



class Flow {
    vector <Particle> data;
    //
    vector<BoundaryPart> boundaries;//index -
    //1- x1, 2-y1, 3-x2, 4-y2, 5-number of 1st type particles
    vector<vector<double>> geometry;
    //
    void setBound(const string& boundaryFile);
    //
    void setInit(const string& initFile);
    //
    void buildBoundPart();
    //
    void buildInitPart(vector<Point>& branch);
    //
    void dataWrite(ostream& os);
    //
    void calculateStep();
public:
    //
    Flow(const string& boundaryFile, const string& initFile);
    //
    void calculate();
    //


};


#endif //SPHSM6_FLOW_H
