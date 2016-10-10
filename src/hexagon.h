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
        axial_to_cube(axial_x, axial_y, &cube_x, &cube_y, &cube_z);
        vorthex = polar_to_vorthex(peel, index);
    }
   ~HEXAGON() {}

    size_t get_peel()   const {return peel;}
    size_t get_index()  const {return index;}

    void to_vertices(std::vector<ALLEGRO_VERTEX> *) const;
    void to_masters (size_t level, std::map<size_t, size_t> *) const;
    bool is_master(size_t) const;

    static void   polar_to_axial  (size_t, size_t, int *, int *);
    static void   cube_to_polar   (int, int, int, size_t *, size_t *);
    static void   axial_to_cube   (int, int, int *, int *, int *);
    static size_t polar_to_vorthex(size_t, size_t);
    static void   vorthex_to_polar(size_t, size_t *, size_t *);

    bool highlight;

    private:
    size_t peel;
    size_t index;
    size_t vorthex;
    int axial_x;
    int axial_y;
    int cube_x;
    int cube_y;
    int cube_z;
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

inline void HEXAGON::to_masters(size_t level, std::map<size_t, size_t> *to) const {
    if (is_master(level)) {
        (*to)[vorthex]++;
        return;
    }
    if (level == 1) {
        size_t p, i;
        cube_to_polar(cube_x,   cube_y+1, cube_z-1, &p, &i); (*to)[polar_to_vorthex(p, i)]++;
        cube_to_polar(cube_x-1, cube_y,   cube_z+1, &p, &i); (*to)[polar_to_vorthex(p, i)]++;
        cube_to_polar(cube_x+1, cube_y-1, cube_z,   &p, &i); (*to)[polar_to_vorthex(p, i)]++;
    }
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
    bool highlight4 = false && is_master(4);
    bool highlight5 = false && is_master(5);
    bool highlight6 = false && is_master(6);
    bool highlight7 = false && is_master(7);

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
            vtx.color = highlight  ? al_map_rgb_f(1.0,0.75, 0.8) :
                        highlight7 ? al_map_rgb_f(0.0, 0.0, 0.0) :
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

inline void HEXAGON::cube_to_polar(int x, int y, int z, size_t *peel, size_t *index) {
    int r = std::max(std::abs(x), std::max(std::abs(y), std::abs(z)));
    *peel = r;

         if (r == y) *index = 0*r-x; // sector 0
    else if (r ==-x) *index = 1*r+z; // sector 1
    else if (r == z) *index = 2*r-y; // sector 2
    else if (r ==-y) *index = 3*r+x; // sector 3
    else if (r == x) *index = 4*r-z; // sector 4
    else if (r ==-z) *index = 5*r+y; // sector 5
    else             *index = 0;
}

inline void HEXAGON::axial_to_cube(int axial_x, int axial_y, int *cube_x, int *cube_y, int *cube_z) {
    *cube_x = axial_x;
    *cube_z = axial_y;
    *cube_y = -axial_x-axial_y;
}

inline size_t HEXAGON::polar_to_vorthex(size_t peel, size_t index) {
    if (peel == 0) return 0;
    size_t prev_peel = peel-1;
    size_t vorthex_size = 3*prev_peel*prev_peel+3*prev_peel+1;
    size_t index_on_rim = index % (peel*6);
    return vorthex_size + index_on_rim;
}

inline void HEXAGON::vorthex_to_polar(size_t vorthex, size_t * peel, size_t * index) {
    if (vorthex == 0) {
        *peel  = 0;
        *index = 0;
        return;
    }
    const double a = 3.0;
    const double b = 3.0;
    double c = 1.0 - vorthex;
    // We solve the quadratic equation f(x) = 3x^2 + 3x + 1 here.
    double temp = -0.5 * (b + (std::signbit(b) ? -1.0 : 1.0) * sqrt(b*b - 4.0*a*c));
    double x1 = temp / a;
    double x2 = c / temp;

    size_t prev_peel = std::max(x1, x2);
    size_t first_index = polar_to_vorthex(prev_peel+1, 0);
    *peel = prev_peel + 1;
    *index = vorthex - first_index;
}

#endif

