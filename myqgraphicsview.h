#ifndef MyQGraphicsView_H
#define MyQGraphicsView_H

#include "Flow.h"


class Flow_Drawer;

class MyQGraphicsView : public QGraphicsView {
Q_OBJECT
public:
    explicit MyQGraphicsView(QWidget *parent = 0);

    ~MyQGraphicsView();

signals:

public slots:

    //drawing
    void mousePressEvent(QMouseEvent *e);

    void ClearWindow();

    void mouseDoubleClickEvent(QMouseEvent *e);

private:
    string file1 = "/home/nikita/sph_interface/txt/testB.txt";
    string file2 = "/home/nikita/sph_interface/txt/testI.txt";
    QGraphicsScene *scene;
    Flow_Drawer *flow;
    //
    int counter_help=0;
public:
    void draw();
    void clear();


};

class Flow_Drawer : public Flow {
    void draw_boundary( QGraphicsScene *scene) {
        const vector<vector<Particle>> *data = reinterpret_cast<const vector<vector<Particle>> *>(this->get_part(0));
        double rad = 2;
        for (int i = 0; i < data->size(); ++i) {
            for (int j = 0; j < data->operator[](i).size(); ++j) {
                Particle *curent = const_cast<Particle *>(&data->operator[](i)[j]);
                scene->addEllipse(curent->X() * 3 - rad, curent->Y() * 3 - rad, rad * 2.0, rad * 2.0,
                                  QPen(Qt::red), QBrush(Qt::SolidPattern));
            }
        }
    }

    void draw_data(QGraphicsScene *scene) {
        const vector<Particle> *data = reinterpret_cast<const vector<Particle> *>(this->get_part(1));
     //   const vector<Particle> *data = (vector<Particle> *)this->get_part(1);
        double rad = 2;
        for (int i = 0; i < data->size(); ++i) {
            const Particle *curent = &data->operator[](i);
            scene->addEllipse(curent->X() * 3 - rad, curent->Y() * 3 - rad, rad * 2.0, rad * 2.0,
                              QPen(Qt::green), QBrush(Qt::SolidPattern));
        }
    }

    void draw_grid(QGraphicsScene *scene) {
        const vector<vector<Cell>> *data = reinterpret_cast<const vector<vector<Cell>>*>(this->get_part(2));
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
                QPointF pt1(l_b.x * 3, l_b.y * 3);
                QPointF pt2(l_b.x * 3, r_t.y * 3);
                QPointF pt3(r_t.x * 3, r_t.y * 3);
                QPointF pt4(r_t.x * 3, l_b.y * 3);
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

    void draw_shadow(QGraphicsScene *scene) {
        const vector<vector<Cell>> *data = reinterpret_cast<const vector<vector<Cell>>*>(this->get_part(2));
        double rad = 2;
        for (int i = 0; i < data->size(); ++i) {
            for (int j = 0; j < data->operator[](i).size(); ++j) {
                const Cell* current = &data->operator[](i)[j];
                const vector<Particle>* curent_shadow = current->get_shadow();
                for(int i = 0; i<curent_shadow->size(); ++i) {
                    scene->addEllipse(curent_shadow->operator[](i).X() * 3 - rad,
                                      curent_shadow->operator[](i).Y() * 3 - rad,
                                      rad * 2.0, rad * 2.0, QPen(Qt::blue), QBrush(Qt::SolidPattern));
                }
            }
        }
    }

public:
    Flow_Drawer(const string &boundaryFile, const string &initFile) :
            Flow(boundaryFile, initFile) {
        calculator = new Calculator(s_distribution);
    }

    virtual void calculate_step(QGraphicsScene *scene) {
        draw_boundary(scene);
        draw_grid(scene);
        draw_shadow(scene);
        //
        calculator->calculate();
        calculations::current_time += calculations::deltaT;
        draw_data(scene);
    }


};


#endif
