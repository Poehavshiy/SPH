//
// Created by nikita on 15.06.16.
//

#ifndef SPHSM6_FLOW_DRAWER_H
#define SPHSM6_FLOW_DRAWER_H

#include "Flow.h"
class QGraphicsScene;

class Flow_Drawer : public Flow {
    void draw_boundary(QGraphicsScene *scene);

    void draw_data(QGraphicsScene *scene) ;

    void draw_grid(QGraphicsScene *scene) ;

    void draw_shadow(QGraphicsScene *scene);

    void draw_data_bycells(QGraphicsScene *scene) ;

public:
    Flow_Drawer(const string &boundaryFile, const string &initFile);

    virtual void calculate_step(QGraphicsScene *scene) ;

};



#endif //SPHSM6_FLOW_DRAWER_H
