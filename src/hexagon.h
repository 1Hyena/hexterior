#ifndef HEXAGON_H_20_09_2016
#define HEXAGON_H_20_09_2016

#include "vector_3d.h"

class HEXAGON {
    public:
    static double size;
    static double gap;

    HEXAGON() {}
    HEXAGON(size_t peel, size_t index) : peel(peel), index(index) {
        polar_to_axial(peel, index, &axial_x, &axial_y);
    }
   ~HEXAGON() {}

    size_t get_peel()  const {return peel;}
    size_t get_index() const {return index;}

    void to_vertices(std::vector<ALLEGRO_VERTEX> *) const;
    bool is_master(size_t) const;

    static void   polar_to_axial (size_t, size_t, int *, int *);
/*    static size_t polar_to_vortex(size_t, size_t);
    static void   vortex_to_polar(size_t, size_t *, size_t *);
*/
    private:
    size_t peel;
    size_t index;
    int axial_x;
    int axial_y;
};

double HEXAGON::size = 1.0;
double HEXAGON::gap  = 0.0;

inline bool HEXAGON::is_master(size_t level) const {
    // Master hexagons have influence over their neighbouring hexagons. This
    // can be used to generate Perlin noise. Each level has a different set of
    // master hexagons. The higher the level the greater the vicinity.
    if (peel == 0 || level == 0) return true;

    if (level % 2 == 0) {
        size_t magic = std::pow(3, level/2);
        if (peel % magic == 0 && (index % magic) == 0) return true;
        return false;
    }

    size_t interval = std::pow(3, (level+1)/2);
    if (peel < 2*(interval/3) || peel % (interval/3) != 0) return false;
    size_t index_on_edge = index % peel;
    if (((index_on_edge + (peel % interval)) % interval) == 0) return true;
    return false;
}

inline void HEXAGON::to_vertices(std::vector<ALLEGRO_VERTEX> * to) const {
    const double distance = 2.0*sqrt(size*size - (size/2.0)*(size/2.0));
    double center_x = 0.0;
    double center_z = 0.0;
    size_t sector   = 0;
    if (peel > 0) {
        sector = index / peel;
        size_t first = sector * peel;
        size_t diff  = index - first;
        double dir_1 = ALLEGRO_PI/2.0 + ((ALLEGRO_PI/3.0) * (sector+1));
        double dir_2 = ALLEGRO_PI/2.0 + ((ALLEGRO_PI/3.0) * (sector+2));
        double p = diff / double(peel);
        center_x = (1.0-p)*sin(dir_1) + p*sin(dir_2);
        center_z = (1.0-p)*cos(dir_1) + p*cos(dir_2);
        center_x*= peel*distance;
        center_z*= peel*distance;
    }

    bool highlight1 = false || is_master(1);
    bool highlight2 = false || is_master(2);
    bool highlight3 = false || is_master(3);
    bool highlight4 = false || is_master(4);
    bool highlight5 = false || is_master(5);
    bool highlight6 = false || is_master(6);
    bool highlight7 = false || is_master(7);

    double dir = ALLEGRO_PI/2.0 + (ALLEGRO_PI/3.0)/2.0;
    for (size_t t=0; t<6; ++t) {
        for (size_t n=0; n<3; ++n) {
            double v_x = center_x;
            double v_z = center_z;
            if (n != 0) {
                double dx = sin(dir);
                double dz = cos(dir);
                v_x += dx*size*(1.0 - gap);
                v_z += dz*size*(1.0 - gap);
            }
            ALLEGRO_VERTEX vtx;
            vtx.x = v_x;
            vtx.y = 0.2;
            vtx.z = v_z;
            vtx.color = highlight7 ? al_map_rgb_f(0.0, 0.0, 0.0) :
                        highlight6 ? al_map_rgb_f(1.0, 1.0, 1.0) :
                        highlight5 ? al_map_rgb_f(0.0, 1.0, 1.0) :
                        highlight4 ? al_map_rgb_f(1.0, 0.0, 1.0) :
                        highlight3 ? al_map_rgb_f(1.0, 1.0, 0.0) :
                        highlight2 ? al_map_rgb_f(0.0, 0.0, 1.0) :
                        highlight1 ? al_map_rgb_f(1.0, 0.0, 0.0) :
                                     al_map_rgb_f(0.0, (t+1)/6.0, 0.0);
            to->push_back(vtx);
            if (n == 1) dir += ALLEGRO_PI/3.0;
        }
    }
}

inline void HEXAGON::polar_to_axial(size_t peel, size_t index, int *x, int *y) {
    if (peel == 0) {
        *x = 0;
        *y = 0;
        return;
    }
    size_t sector = index / peel;
    size_t index_on_edge = index % peel;
    switch (sector) {
        case   0: *x = -index_on_edge;      *y = -peel + index_on_edge; break;
        case   1: *x = -peel;               *y = index_on_edge;         break;
        case   2: *x = -peel+index_on_edge; *y = peel;                  break;
        case   3: *x = index_on_edge;       *y = peel - index_on_edge;  break;
        case   4: *x = peel;                *y = -index_on_edge;        break;
        case   5: *x = peel-index_on_edge;  *y = -peel;                 break; 
        default : break;
    }
}
/*
inline size_t HEXAGON::polar_to_vortex(size_t peel, size_t index) {
    // f(x)=6*(((x-1)*(x+2)+2)/2)+1
    // where x is the number of peels, it gives the number of hexagons in total
    if (peel == 0) return 0;
    size_t index_on_edge = index % peel;

    

    1            =  1                1   
    1+6          =  7                6  6                                        1
    1+6+12       = 19               18  6 6 6                                    3 
    1+6+12+18    = 37               36  6 6 6 6 6 6                              6
    1+6+12+18+24 = 61               60  6 6 6 6 6 6 6 6 6 6                     10
    1+6+12+18+24+30 = 91            90  6 6 6 6 6 6 6 6 6 6 6 6 6 6 6           15
                   +36 = 127    
}

inline void HEXAGON::vortex_to_polar(size_t vortex, size_t * peel, size_t * index) {

}*/

#endif

