//
// Created by nikita on 15.06.16.
//

#ifndef SPHSM6_FLOW_DRAWER_H
#define SPHSM6_FLOW_DRAWER_H

#include "Flow.h"
class QGraphicsScene;

class Flow_Drawer : public Flow {
    /*
   Return a RGB colour value given a scalar v in the range [vmin,vmax]
   In this case each colour component ranges from 0 (no contribution) to
   1 (fully saturated), modifications for other ranges is trivial.
   The colour is clipped at the end of the scales if v is outside
   the range [vmin,vmax]
*/

    int iteration = 0;
    int critical_iter = 0;
    int crit_i, crit_j;
    typedef struct {
        double r,g,b;
    } COLOUR;

    COLOUR GetColour(double v,double vmin,double vmax);

    void show_iter(QGraphicsScene *scene);

    void draw_dataP(QGraphicsScene *scene);

    void draw_boundary(QGraphicsScene *scene);

    void draw_data(QGraphicsScene *scene) ;

    void draw_grid(QGraphicsScene *scene) ;

    void draw_shadow(QGraphicsScene *scene);

    void draw_data_bycells(QGraphicsScene *scene) ;

    void max_minP(QGraphicsScene *scene );

public:
    Flow_Drawer(const string &boundaryFile, const string &initFile);

    virtual void calculate_step(QGraphicsScene *scene) ;

};



#endif //SPHSM6_FLOW_DRAWER_H
