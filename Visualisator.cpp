//
// Created by nikita on 12.09.16.
//

#include "Visualisator.h"
/*
x coord, y coord, z coord, scalar
0, 0, 0, 0
1, 0, 0, 1
0, 1, 0, 2
1, 1, 0, 3
-0.5, -0.5, 1, 4
0.5, -0.5, 1, 5
-0.5, 0.5, 1, 6
0.5, 0.5, 1, 7
*/


//
Visualisator::Visualisator(Flow *to_vis, string dir):counter(0), path(dir),target(to_vis) {
}

void Visualisator::write_step() {
    stringstream ss;
    ss << counter;
    string str_count = ss.str();
    std::ofstream outfile( path + name + str_count);

    for(int i=0; i<target->data.size(); ++i) {
        outfile<<target->data[i].X()<<", ";
        outfile<<target->data[i].Y()<<", ";
        outfile<<"0, ";
        outfile<<target->data[i].P();
        outfile<<endl;
    }
    counter++;
    outfile.close();
}