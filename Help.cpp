//
// Created by nikita on 13.05.16.
//

#include "Help.h"
#include "SpaceParsing.h"

//1- x1, 2-y1, 3-x2, 4-y2, 5-number of 1st type particles //seems ok
void ReadWrite::parse_bound(string &current, vector<vector<double>> &geometry) {
    int j = 0, signs = 0;
    int n = 0;
    vector<double> curG;
    curG.resize(5);
    for (int i = 0; i <= current.size(); ++i) {
        // char cur=current[i];
        //help+=cur;
        if (current[i] == ' ' || current[i] == '\0') {
            j = i - signs;
            int cur = stoi(current.substr(j, i));
            curG[n] = cur;
            signs = 0;
            ++n;
        }
        ++signs;
    }
    geometry.push_back(curG);
}

//
Point ReadWrite::parse_init(const string &current) {
    double x;
    double y;
    return Point(0, 0);
}

////seems ok
void Help::fill_bound_line(vector<Particle> &current, vector<double> &line) {
    double xstep = (line[2] - line[0]) / line[4];
    double ystep = (line[3] - line[1]) / line[4];
    //
    for (int i = 0; i < current.size(); ++i) {
        current[i].set_pos(line[0] + xstep * i, line[1] + ystep * i);
    }
    current.push_back(Particle(1));
    current[current.size() - 1].set_pos(line[2], line[3]);
}
//