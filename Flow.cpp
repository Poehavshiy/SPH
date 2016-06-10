//
// Created by nikita on 13.05.16.
//

#include "Flow.h"
#include "Help.h"

Flow::Flow(const string &boundaryFile, const string &initFile) {
    set_bound(boundaryFile);

    set_init(initFile);

    s_distribution = SpaceParsing::init(geometry, boundaries, data);
}

////parsing boundary and initial conditions files
void Flow::set_bound(const string &boundaryFile) {
    fstream bound;
    bound.open(boundaryFile);
    string current;
    while (std::getline(bound, current)) {
        if (current.find_first_not_of('*') != '/') {
            ReadWrite::parse_bound(current, geometry);
        }
    }
    build_bound_part();
    bound.close();
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
    vector<Point> branch;
    double x, y;
    while (initCond >> x >> y) {
        // if(current.find_first_not_of('*')!='/') {
        Point cur(x, y);
        branch.push_back(cur);
        //}
    }
    build_init_part(branch);
    initCond.close();
}

//fill branc with particles with initial patameters
void Flow::build_init_part(vector<Point> &branch) {
    int number = 100;
    //стороны прямоугольника, который я заполняю частицами
    double X = abs(branch[0].x - branch[3].x);
    double Y = abs(branch[0].y - branch[1].y);
   // double XYratio = X / Y;
    int per_y = sqrt(number);
    int per_x = per_y;//number / per_y;
    double xstep = X / per_x;
    double ystep = Y / (per_y-1);
    //
    double P = 100;
    double p = 2;
    double mass = 1;
    double e = P / (0.4 * p);
    for (int i = 0; i < number; ++i) {
        data.push_back(Particle(0, p, P, e, 10, 0, mass));
    }
    for (int i = 0; i <= per_y; ++i) {
        for (int j = 0; j <= per_x; ++j) {
            double curX = branch[0].x + j * xstep;
            double curY = branch[0].y + i * ystep;
            data[per_x * i + j].set_pos(curX, curY);
        }
    }

}
//
void Flow::calculate_step() {
    for (int i = 0; i < data.size(); ++i) {
        double P = data[i].P();
        P = rand() % 150 + 1;
        data[i].set_pressure(P);
    }
}

//
void Flow::calculate() {
<<<<<<< HEAD

    calculator = new Calculator(s_distribution);
    while (calculations::current_time < 1) {
        calculator->calculate();
        calculations::current_time+=calculations::deltaT;
    }
=======
    ofstream check("/home/nikita/SPHSm6/result.txt");
    ofstream result("/home/nikita/SPHSm6/result1.txt");
    calculator = new Calculator(s_distribution);
    for (int i = 0; i < 1; ++i) {
        //calculate_step();
//        ReadWrite::data_write(result, data);
        check<<"next"<<endl;
        check<<*s_distribution;
   //     ReadWrite::data_write(result, boundaries);
    }

    while (calculations::current_time < 1) {
       // cout<<calculations::current_time<<endl;
        //cout<<data.size()<<endl;
        //calculator->calculate();
        calculations::current_time+=calculations::deltaT;
        for(int i = 0; i<data.size() ; ++i) {
            double newx=data[i].X()+0.1;
            double newy=data[i].Y()+0.1;
            data[i].set_pos(newx, newy);
        }

        ReadWrite::data_write(result, data);
        ReadWrite::data_write(result, boundaries);
    }
    result.close();
>>>>>>> 883cc5708c435c9115b864d6fdf26606dffc1703
}
