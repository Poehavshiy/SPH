//
// Created by nikita on 13.05.16.
//

#include "Help.h"
#include "SpaceParsing.h"
#include <sstream>
#include <stdio.h>
#include <stdlib.h>



//1- x1, 2-y1, 3-x2, 4-y2, 5-number of 1st type particles //seems ok
void ReadWrite::parse_bound(string &current, vector<vector<double>> &geometry) {
    std::string segment;
    std::vector<std::string> seglist;
    vector<double> curG(5);
    std::stringstream test(current);
    while(std::getline(test, segment, ' '))
    {
        seglist.push_back(segment);
    }
    std::string::size_type sz;
    for(int i = 0; i < 5; ++i){
        const char* check_c = seglist[i].c_str();
        double check = atof(check_c);//stod(seglist[i], &sz);
        curG[i] = check;
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