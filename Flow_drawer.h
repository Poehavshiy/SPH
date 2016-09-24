//
// Created by nikita on 15.06.16.
//

#ifndef SPHSM6_FLOW_DRAWER_H
#define SPHSM6_FLOW_DRAWER_H

#include "Flow.h"
class QGraphicsScene;

class Flow_Drawer : public Flow {
    int cof; int cof1;
    /*
   Return a RGB colour value given a scalar v in the range [vmin,vmax]
   In this case each colour component ranges from 0 (no contribution) to
   1 (fully saturated), modifications for other ranges is trivial.
   The colour is clipped at the end of the scales if v is outside
   the range [vmin,vmax]
*/
    typedef struct {
        double r,g,b;
    } COLOUR;

    string mathplot_path = "/home/nikita/SPH/python/";

    COLOUR GetColour(double v,double vmin,double vmax);

    void show_information(QGraphicsScene *scene);

    void draw_dataP(QGraphicsScene *scene);

    void draw_boundary(QGraphicsScene *scene);

    void draw_data(QGraphicsScene *scene) ;

    void draw_grid(QGraphicsScene *scene) ;

    void draw_shadow(QGraphicsScene *scene);

    void draw_data_bycells(QGraphicsScene *scene) ;

    void set_text(QString&& text, QGraphicsScene *scene,pair<int, int>&& position, double&& value = 0);

    void write_py_data(int target_itertion);
public:
    Flow_Drawer(const string &boundaryFile, const string &initFile);

    virtual void calculate_step(QGraphicsScene *scene) ;

};



#endif //SPHSM6_FLOW_DRAWER_H
