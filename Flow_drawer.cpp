//
// Created by nikita on 15.06.16.
//

#include "Flow_drawer.h"
#include "myqgraphicsview.h"

Flow_Drawer::COLOUR Flow_Drawer::GetColour(double v,double vmin,double vmax)
{
    COLOUR c = {1.0,1.0,1.0}; // white
    double dv;

    if (v < vmin)
        v = vmin;
    if (v > vmax)
        v = vmax;
    dv = vmax - vmin;

    if (v < (vmin + 0.25 * dv)) {
        c.r = 0;
        c.g = 4 * (v - vmin) / dv;
    } else if (v < (vmin + 0.5 * dv)) {
        c.r = 0;
        c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
    } else if (v < (vmin + 0.75 * dv)) {
        c.r = 4 * (v - vmin - 0.5 * dv) / dv;
        c.b = 0;
    } else {
        c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
        c.b = 0;
    }

    return(c);
}


Flow_Drawer::Flow_Drawer(const string &boundaryFile, const string &initFile) :
        Flow(boundaryFile, initFile) {
    calculator = new Calculator(s_distribution);
}

void Flow_Drawer::draw_boundary(QGraphicsScene *scene) {
    const vector<vector<Particle>> *data = reinterpret_cast<const vector<vector<Particle>> *>(this->get_part(0));
    double rad = 2;
    for (int i = 0; i < data->size(); ++i) {
        for (int j = 0; j < data->operator[](i).size(); ++j) {
            Particle *curent = const_cast<Particle *>(&data->operator[](i)[j]);
            scene->addEllipse(curent->X() - rad, curent->Y() - rad, rad * 2.0, rad * 2.0,
                              QPen(Qt::black), QBrush(Qt::SolidPattern));
        }
    }
}

void Flow_Drawer::draw_data(QGraphicsScene *scene) {
    const vector<Particle> *data = reinterpret_cast<const vector<Particle> *>(this->get_part(1));
    //   const vector<Particle> *data = (vector<Particle> *)this->get_part(1);
    double rad = 2;
    for (int i = 0; i < data->size(); ++i) {
        const Particle *curent = &data->operator[](i);
        scene->addEllipse(curent->X() - rad, curent->Y() - rad, rad * 2.0, rad * 2.0,
                          QPen(Qt::green), QBrush(Qt::SolidPattern));
    }
}

void Flow_Drawer::draw_grid(QGraphicsScene *scene) {
    const vector<vector<Cell>> *data = reinterpret_cast<const vector<vector<Cell>> *>(this->get_part(2));
    double rad = 2;
    for (int i = 0; i < data->size(); ++i) {
        for (int j = 0; j < data->operator[](i).size(); ++j) {
            const Cell *cur = &data->operator[](i)[j];
            const Point l_b = cur->get_left_bottom();
            const Point r_t = cur->get_right_top();
            //
            int colour[3] = {0, 0, 0};
            colour[1] = 128;
            QColor myC(colour[0], colour[1], colour[2]);
            QPointF pt1(l_b.x, l_b.y);
            QPointF pt2(l_b.x, r_t.y);
            QPointF pt3(r_t.x, r_t.y);
            QPointF pt4(r_t.x, l_b.y);
            //
            QLineF lt1(pt1, pt2);
            QLineF lt2(pt2, pt3);
            QLineF lt3(pt3, pt4);
            QLineF lt4(pt4, pt1);
            //
            scene->addLine(lt1, QPen(myC));
            scene->addLine(lt2, QPen(myC));
            scene->addLine(lt3, QPen(myC));
            scene->addLine(lt4, QPen(myC));

            //
        }
    }
}

void Flow_Drawer::draw_shadow(QGraphicsScene *scene) {
    const vector<vector<Cell>> *data = reinterpret_cast<const vector<vector<Cell>> *>(this->get_part(2));
    double rad = 2;
    int colour[3] = {0, 0, 0};
    for (int i = 0; i < data->size(); ++i) {
        for (int j = 0; j < data->operator[](i).size(); ++j) {
            const Cell *current = &data->operator[](i)[j];
            const vector<Particle> *curent_shadow = current->get_shadow();
            for (int i = 0; i < curent_shadow->size(); ++i) {
                scene->addEllipse(curent_shadow->operator[](i).X() - rad,
                                  curent_shadow->operator[](i).Y() - rad,
                                  rad * 2.0, rad * 2.0, QPen(QColor(colour[0], colour[1], colour[2])),
                                  QBrush(Qt::SolidPattern));
            }
        }
    }
}

void Flow_Drawer::draw_data_bycells(QGraphicsScene *scene) {
    const vector<vector<Cell>> *data = reinterpret_cast<const vector<vector<Cell>> *>(this->get_part(2));
    double rad = 2;
    int c[6][3] = {
            {255, 0,   0},
            {0,   255, 0},
            {0,   0,   255},
            {255,   255,   0},
            {0,   255,   255},
            {255,   0,   255}
    };
    QColor colors[6];
    for (int i = 0; i < 6; ++i) {
        colors[i].setRgb(c[i][0], c[i][1], c[i][2]);
    }
    int colour_counter = 0;
    for (int i = 0; i < data->size(); ++i) {
        for (int j = 0; j < data->operator[](i).size(); ++j) {
            if (colour_counter == 6) colour_counter = 0;
            const Cell *current = &data->operator[](i)[j];
            const vector<Particle *> *curent_real = current->get_real();
            for (int i = 0; i < curent_real->size(); ++i) {
                scene->addEllipse(curent_real->operator[](i)->X() - rad,
                                  curent_real->operator[](i)->Y() - rad,
                                  rad * 2.0, rad * 2.0, colors[colour_counter], QBrush(Qt::SolidPattern));
            }
            if(curent_real->size()!=0) ++colour_counter;
        }
    }
}

void Flow_Drawer::draw_dataP(QGraphicsScene *scene){
    const vector<vector<Cell>> *data = reinterpret_cast<const vector<vector<Cell>> *>(this->get_part(2));
    double rad = 2;
    for (int i = 0; i < data->size(); ++i) {
        for (int j = 0; j < data->operator[](i).size(); ++j) {
            const Cell *current = &data->operator[](i)[j];
            const vector<Particle *> *curent_real = current->get_real();
            for (int i = 0; i < curent_real->size(); ++i) {
                double value = curent_real->operator[](i)->P();
                COLOUR currentC = GetColour(value, 0 ,maxP);
                QColor QTcurrentC;
                QTcurrentC.setRgb(255*currentC.r,255*currentC.g, 255*currentC.b);
                scene->addEllipse(curent_real->operator[](i)->X() - rad,
                                  curent_real->operator[](i)->Y() - rad,
                                  rad * 2.0, rad * 2.0, QTcurrentC , QBrush(Qt::SolidPattern));
            }
        }
    }
}

void Flow_Drawer::max_minP(QGraphicsScene *scene ){
    const vector<vector<Cell>> *data = reinterpret_cast<const vector<vector<Cell>> *>(this->get_part(2));
    double rad = 2;
    double Pmax_cur = -10050000;
    double Pmin_cur = 10050000;

    for (int i = 0; i < data->size(); ++i) {
        for (int j = 0; j < data->operator[](i).size(); ++j) {
            const Cell *current = &data->operator[](i)[j];
            const vector<Particle *> *curent_real = current->get_real();
            for (int i = 0; i < curent_real->size(); ++i) {
                double P = curent_real->operator[](i)->P();
                if(P > Pmax_cur) Pmax_cur = P;
                if(P < Pmin_cur) Pmin_cur = P;
                if(Pmax_cur > maxP && critical_iter==0) {
                    crit_i = i;
                    crit_j = j;
                    critical_iter = iteration;

                }
            }
        }
    }
    QString PmaxStr = QString::number(Pmax_cur);
    QString PminStr = QString::number(Pmin_cur);
    QString res = "Pmax=" + PmaxStr + "\n" + "Pmin=" + PminStr;
    QGraphicsTextItem *text = scene->addText(res);
    text->setPos(-50, -100);
};

void Flow_Drawer::show_iter(QGraphicsScene *scene) {
    QString res ="iteration=" + QString::number(iteration);
    QGraphicsTextItem *text = scene->addText(res);
    text->setPos(-50, -120);

    res ="crit_iteration=" + QString::number(critical_iter) + " i=" + QString::number(crit_i) + " j=" + QString::number(crit_j);
    QGraphicsTextItem *text1 = scene->addText(res);
    text1->setPos(-50, -140);
}

void Flow_Drawer::calculate_step(QGraphicsScene *scene) {
    draw_boundary(scene);
    draw_grid(scene);
    draw_shadow(scene);
    //

    max_minP(scene);

    show_iter(scene);
    //
    calculator->calculate();
    draw_dataP(scene);
    // cout<<data.size();
    iteration++;
}