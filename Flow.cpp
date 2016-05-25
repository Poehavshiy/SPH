//
// Created by nikita on 13.05.16.
//

#include "Flow.h"
#include "Help.h"


Flow::Flow(const string &boundaryFile, const string &initFile) {
    setBound(boundaryFile);

    setInit(initFile);
}
////parsing boundary and initial conditions files
void Flow::setBound(const string& boundaryFile) {
    fstream bound;
    bound.open(boundaryFile);
    string current;
    while (std::getline(bound, current)) {
        if (current.find_first_not_of('*') != '/') {
            ReadWrite::parseBound(current, geometry);
        }
    }

    buildBoundPart();
    bound.close();
}
//
void Flow::buildBoundPart() {

    for (int i = 0; i < geometry.size(); ++i) {
        int size = geometry[i][4];
        vector<Particle> current(size, Particle(1));//размер и тип ячейки
        Help::fillBoundLine(current, geometry[i]);
        double x = abs(geometry[i][0] - geometry[i][2]);
        double y = abs(geometry[i][1] - geometry[i][3]);
        double r = sqrt(x * x + y * y);
        boundaries.push_back(BoundaryPart(current, r));
    }
}
//
void Flow::setInit(const string& initFile) {
    fstream initCond;
    initCond.open(initFile);
    vector<Point> branch;
    double x, y;
    while (initCond>>x>>y) {
        // if(current.find_first_not_of('*')!='/') {
        Point cur(x, y);
        branch.push_back(cur);
        //}
    }
    buildInitPart(branch);
    initCond.close();
}
//fill branc with particles with initial patameters
void Flow::buildInitPart(vector<Point> &branch) {
    int number = 20;
    double X = abs(branch[0].x - branch[3].x);
    double Y = abs(branch[0].y - branch[1].y);
    double XYratio=X/Y;
    int per_y=2*sqrt(number)/(XYratio+1);
    int per_x=number/per_y;
    double xstep=X/per_x;
    double ystep=Y/per_y;
    //
    double P = 1;
    double p = 0.1;
    double mass = 0.5;
    double e = P / (0.4 * p);
    for (int i = 0; i < number; ++i) {
        data.push_back(Particle(0, p, P, e, 10, 1, mass));
    }
    for(int i=0; i<per_y; ++i) {
        for (int j = 0; j < per_x; ++j) {
            double curX = branch[0].x+j * xstep;
            double curY = branch[0].y+i * ystep;
            data[per_x * i + j].set_pos(curX, curY);
        }
    }

}
//
void Flow::calculateStep() {
    for(int i=0; i<data.size(); ++i) {
        double P=data[i].P();
        P=rand() % 150 +1;
        data[i].set_pressure(P);
    }
}
//
void Flow::calculate() {
    fstream result;
    result.open("/home/nikita/SPHSm6/result.txt");
    for(int i=0; i<100; ++i) {
        calculateStep();
        dataWrite(result);
    }
    result.close();
}
//
void Flow::dataWrite(ostream& os) {
    os<<"next"<<endl;
    for(int i=0; i<data.size(); ++i) {
        os<<data[i];
    }
}