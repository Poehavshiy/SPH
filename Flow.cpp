//
// Created by nikita on 13.05.16.
//

#include "Flow.h"
#include "Help.h"

Particle *for_debugin;

Flow::Flow(const string &boundaryFile, const string &initFile) {
    int per_x, per_y;

    N = 2000;

    sphere = false;

    set_bound(boundaryFile, per_x, per_y);

    set_init(initFile);

    s_distribution = SpaceParsing::init(geometry, boundaries, data, per_x, per_y);
}

////parsing boundary and initial conditions files
void Flow::set_bound(const string &boundaryFile, int &x, int &y) {
    fstream bound;
    bound.open(boundaryFile);
    string current;
    int per_x, per_y;
    bound >> per_x >> per_y;
    getline(bound, current);
    while (std::getline(bound, current)) {
        if (current.find_first_not_of('*') != '/') {
            ReadWrite::parse_bound(current, geometry);
        }
    }
    build_bound_part();
    bound.close();
    x = per_x;
    y = per_y;
}

//
void Flow::build_bound_part() {

    for (int i = 0; i < geometry.size(); ++i) {
        int size = geometry[i][4];
        vector<Particle> current(size, Particle(1));//размер и тип ячейки
        Help::fill_bound_line(current, geometry[i]);
        boundaries.push_back(current);
    }
}

//
void Flow::set_init(const string &initFile) {
    fstream initCond;
    initCond.open(initFile);
    vector<vector<Point>> branches;
    double x, y;
    string line;
    string point_str;
    int count = 0;
    while (std::getline(initCond, line)) {
        if (line == "Branch") {
            std::getline(initCond, point_str);
            vector<Point> branch_cur;
            do {
                initCond >> x >> y;
                Point cur(x, y);
                if (x != 999) branch_cur.push_back(cur);
            } while (x != 999);
            count++;
            branches.push_back(branch_cur);
        }

    }
    vector<vector<double>> conditions;
    maxP = 10;
    maxp = 1;
    maxe = maxP / (0.4 * maxp);
    double r = 20;
    double height = 144;
    double M;
    if(sphere == false) M = maxp* r *r* height;
    else M = maxp*PI * pow(r, 3)*4/3 ;
    double Mi = M / N;
    conditions = {
            {maxP,      maxp, Mi, 0, 0},
            {maxP / 10, maxp, Mi, 0, 0}
    };
    if(sphere == false) {
        for (int i = 0; i < branches.size(); ++i) {
            build_init_part(branches[i], conditions[i], i);
        }
    }
    else {
        for (int i = 0; i < branches.size(); ++i) {
            build_init_part(branches[i], conditions[i], i, r);
        }
    }

    initCond.close();
}

//fill branc with particles with initial patameters
void Flow::build_init_part(vector<Point> &branch, vector<double> &cond, int step, double r) {
    int number = N;
    //стороны прямоугольника, который я заполняю частицами
    double X = abs(branch[0].x - branch[3].x);
    double Y = abs(branch[0].y - branch[1].y);
    // double XYratio = X / Y;
    int per_y = sqrt(number);
    int per_x = per_y;//number / per_y;
    double xstep = X / per_x;
    double ystep = Y / per_y;
    //
    double P = cond[0];
    double p = cond[1];
    double mass = cond[2];
    double Vx = cond[3];
    double Vy = cond[4];
    double e = P / (0.4 * p);
    vector<Particle> buffer;
    for (int i = 0; i < number; ++i) {
        buffer.push_back(Particle(0, p, P, e, Vx, Vy, mass));
    }

    for (int i = 0; i <= per_y; ++i) {
        for (int j = 0; j <= per_x; ++j) {
            double curX = branch[0].x + j * xstep;
            double curY = branch[0].y + i * ystep;
            buffer[per_x * i + j].set_pos(curX, curY);
        }
    }
    double center_x = 300;
    double center_y = 300;
    for (int i = 0; i < number; ++i) {
        if(r > 0) {
            if (hypot(center_x - buffer[i].X(), center_y - buffer[i].Y()) < r)
                data.push_back(buffer[i]);
        }
        else
            data.push_back(buffer[i]);
    }

    for_debugin = &data[71];

}
//
void Flow::build_init_part(vector<Point> &branch, vector<double> &cond, int step) {
    int number = N;
    //стороны прямоугольника, который я заполняю частицами
    double X = abs(branch[0].x - branch[3].x);
    double Y = abs(branch[0].y - branch[1].y);
    // double XYratio = X / Y;
    int per_y = sqrt(number);
    int per_x = per_y;//number / per_y;
    double xstep = X / (per_x+1);
    double ystep = Y / (per_y+1) ;
    //
    double P = cond[0];
    double p = cond[1];
    double mass = cond[2];
    double Vx = cond[3];
    double Vy = cond[4];
    double e = P / (0.4 * p);
    for (int i = 0; i < number; ++i) {
        data.push_back(Particle(0, p, P, e, Vx, Vy, mass));
    }

    for (int i = 0; i <= per_y+1; ++i) {
        for (int j = 0; j <= per_x+1; ++j) {
            double curX = branch[0].x + j * xstep;
            double curY = branch[0].y + i * ystep;
            data[step*number+per_x * i + j].set_pos(curX, curY);
        }
    }
    for_debugin = &data[71];

}