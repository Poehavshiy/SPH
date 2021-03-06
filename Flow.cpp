//
// Created by nikita on 13.05.16.
//

#include "Flow.h"
#include "Help.h"

Particle *for_debugin;

void Flow::test_W() {
    std::ofstream outfile("/home/nikita/ClionProjects/dummy_sph/python/check_W.txt");
    outfile << "begin.\n";
    double x = -2;
    double sum = 0;
    double dx = 0.01;
    while (x < 2) {
        double value = W_functions::dW_disser(x, 1);
        outfile << x << " " << value << '\n';
        x += dx;
        sum += dx * value;
    }
// assert(fabs(sum - 1) < 0.001);
    outfile << "end.\n";
}



Flow::Flow(const string &boundaryFile, const string &initFile) {
    test_W();
    int per_x, per_y;

    N = 1500;

    sphere = true;

    set_bound(boundaryFile, per_x, per_y);

    set_init(initFile);

    s_distribution = SpaceParsing::init(geometry, boundaries, data, per_x, per_y, h);
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
    vector<bool> shapes;
    double x, y;
    string line;
    string point_str;
    int count = 0;
    while (std::getline(initCond, line))
    {
        if(line == "Rectangle" || line == "Circle") {
            if(line == "Rectangle"){
                shapes.push_back(false);
            }
            else shapes.push_back(true);

            std::getline(initCond, point_str);
            vector<Point> branch_cur;
            do {
                initCond >> x >> y;
                Point cur(x, y);
                if(x!=999) branch_cur.push_back(cur);
            }
            while(x!=999);
            count ++;
            branches.push_back(branch_cur);
        }
    }
    vector<vector<double>> conditions;

    // шар мужиков r= 0.1
//    maxp=1.63;
//    maxe = 4.29 * 10E6;
//    maxP = 0.4 * maxp * maxe;
  //  maxp = 1;
  //  maxP = 1;

    maxp = 1;
    maxP = 1;

    conditions = {
            {maxP, maxp, 1, 0, 0},
            {maxP/10, maxp/10, 1, 0, 0}
    };
    for(int i = 0; i < branches.size(); ++i) {
        build_init_part(branches[i], conditions[i], i, shapes[i]);
    }
    initCond.close();
}

//fill branc with particles with initial patameters
//0- rectangle 1-circle
void Flow::build_init_part(vector<Point> &branch,vector<double>& cond, int step, bool shape) {
    int number = 200;
    //стороны прямоугольника, который я заполняю частицами
    double X = abs(branch[0].x - branch[3].x);
    double Y = abs(branch[0].y - branch[1].y);
    double XYratio = X / Y;
    int per_y = sqrt(number * Y / X);
    int per_x = per_y * XYratio;//number / per_y;
    double xstep = X / per_x;
    double ystep = Y / per_y ;
    //
    h = 4 * xstep;
    //
    double P = cond[0];
    double p = cond[1];
    double mass = cond[2];
    if(shape){
        double r = (branch[2].x - branch[0].x) / 2;
        mass = (p * PI * pow(r, 2) / 2) / number;
    }
    else{
        mass = (p * X * Y) / number;
    }


    double Vx = cond[3];
    double Vy = cond[4];
    double e = P / (0.4 * p);
    for (int i = 0; i < number; ++i) {
        data.push_back(Particle(0, p, P, e, Vx, Vy, mass));
    }

    for (int i = 0; i <= per_y; ++i) {
        for (int j = 0; j <= per_x; ++j) {
            double curX = branch[0].x + j * xstep;
            double curY = branch[0].y + i * ystep;
            data[step*number+per_x * i + j].set_pos(curX, curY);
        }
    }
    if (shape){
        vector<Particle> data_sphere;
        double center_x = (branch[2].x + branch[0].x) / 2;
        double center_y = (branch[2].y + branch[0].y) / 2;
        double r = (branch[2].x - branch[0].x) / 2;
        for (int i = 0; i < number; ++i)
        {
            double curR = hypot((center_x - data[i].X()) , (center_y - data[i].Y()));
            if (curR <= r)
                data_sphere.push_back(data[i]);
        }
        data = data_sphere;
    }
    for_debugin = &data[71];

}
