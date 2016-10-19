#ifndef HEXAGON_H_20_09_2016
#define HEXAGON_H_20_09_2016

#include <random>
#include <set>

#include "vector_3d.h"

class HEXAGON {
    public:
    static double size;
    static double gap;

    HEXAGON() {}
    HEXAGON(size_t peel, size_t index) : peel(peel), index(index) {
        polar_to_axial(peel, index, &axial_x, &axial_y);
        axial_to_cube(axial_x, axial_y, &cube_x, &cube_y, &cube_z);
        vortex = polar_to_vortex(peel, index);
        highlight = false;
        height    = 0.0;
    }
   ~HEXAGON() {}

    size_t get_peel()        const {return peel;}
    size_t get_index()       const {return index;}
    size_t get_vortex()      const {return vortex;}
    int    get_x()           const {return axial_x;}
    int    get_y()           const {return axial_y;}
    double get_noise()       const;
    double get_noise(size_t) const;

    void to_vertices(std::vector<ALLEGRO_VERTEX> *) const;
    void to_masters (size_t level, std::map<size_t, size_t> *) const;
    bool is_master(size_t) const;

    double get_height() const   {return height;}
    void   set_height(double h) {height = h;}

    static void   polar_to_axial (size_t, size_t, int *, int *);
    static void   cube_to_polar  (int, int, int, size_t *, size_t *);
    static void   axial_to_cube  (int, int, int *, int *, int *);
    static size_t polar_to_vortex(size_t, size_t);
    static void   vortex_to_polar(size_t, size_t *, size_t *);
    static bool   check_level    (size_t, size_t, size_t);
    static size_t cube_distance  (int, int, int, int, int, int);

    bool highlight;

    private:
    size_t peel;
    size_t index;
    size_t vortex;
    int axial_x;
    int axial_y;
    int cube_x;
    int cube_y;
    int cube_z;

    mutable std::mt19937 gen{0};
    mutable std::uniform_real_distribution<double> dis{0, 1};

    double height;
};

double HEXAGON::size = 1.0;
double HEXAGON::gap  = 0.0;

inline double HEXAGON::get_noise() const {
    return dis(gen);
}

inline double HEXAGON::get_noise(size_t level) const {
    gen.seed(vortex);
    gen.discard(level);
    return dis(gen);
}

inline bool HEXAGON::check_level(size_t peel, size_t index, size_t level) {
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

inline bool HEXAGON::is_master(size_t level) const {
    return check_level(peel, index, level);
}

inline void HEXAGON::to_masters(size_t level, std::map<size_t, size_t> *to) const {
    size_t p, i;
    if (level == 0) return;

    std::map<size_t, size_t> masters;

    if (is_master(level)) {
      /*  size_t power;
        (*to)[vortex]+=3;
        if (level % 2) {
            power = std::pow(3, level/2);
            cube_to_polar(cube_x+2*power, cube_y-1*power, cube_z-1*power, &p, &i); (*to)[polar_to_vortex(p, i)]++;
            cube_to_polar(cube_x-1*power, cube_y+2*power, cube_z-1*power, &p, &i); (*to)[polar_to_vortex(p, i)]++;
            cube_to_polar(cube_x-1*power, cube_y-1*power, cube_z+2*power, &p, &i); (*to)[polar_to_vortex(p, i)]++;
            cube_to_polar(cube_x+1*power, cube_y+1*power, cube_z-2*power, &p, &i); (*to)[polar_to_vortex(p, i)]++;
            cube_to_polar(cube_x-2*power, cube_y+1*power, cube_z+1*power, &p, &i); (*to)[polar_to_vortex(p, i)]++;
            cube_to_polar(cube_x+1*power, cube_y-2*power, cube_z+1*power, &p, &i); (*to)[polar_to_vortex(p, i)]++;
        }
        else {
            power = std::pow(3, (level+1)/2);
            cube_to_polar(cube_x+1*power, cube_y,         cube_z-1*power, &p, &i); (*to)[polar_to_vortex(p, i)]++;
            cube_to_polar(cube_x-1*power, cube_y+1*power, cube_z,         &p, &i); (*to)[polar_to_vortex(p, i)]++;
            cube_to_polar(cube_x,         cube_y-1*power, cube_z+1*power, &p, &i); (*to)[polar_to_vortex(p, i)]++;
            cube_to_polar(cube_x,         cube_y+1*power, cube_z-1*power, &p, &i); (*to)[polar_to_vortex(p, i)]++;
            cube_to_polar(cube_x-1*power, cube_y,         cube_z+1*power, &p, &i); (*to)[polar_to_vortex(p, i)]++;
            cube_to_polar(cube_x+1*power, cube_y-1*power, cube_z,         &p, &i); (*to)[polar_to_vortex(p, i)]++;
        }
        return;
    */
        /*
        cube_to_polar(cube_x+1, cube_y,             cube_z-1, &p, &i); masters[polar_to_vortex(p, i)]++;
        cube_to_polar(cube_x-1, cube_y+1,           cube_z,   &p, &i); masters[polar_to_vortex(p, i)]++;
        cube_to_polar(cube_x,             cube_y-1, cube_z+1, &p, &i); masters[polar_to_vortex(p, i)]++;
        cube_to_polar(cube_x,             cube_y+1, cube_z-1, &p, &i); masters[polar_to_vortex(p, i)]++;
        cube_to_polar(cube_x-1, cube_y,             cube_z+1, &p, &i); masters[polar_to_vortex(p, i)]++;
        cube_to_polar(cube_x+1, cube_y-1,           cube_z,   &p, &i); masters[polar_to_vortex(p, i)]++;
        */
        (*to)[vortex]++;
        return;
    }
    else masters[vortex] = 1;

    for (size_t n=1; n<=level; ++n) {
        std::map<size_t, size_t> buf;
        size_t power;

        if (n % 2 != 0) {
            power = std::pow(3, n/2);
            for (const auto & a : masters) {

                vortex_to_polar(a.first, &p, &i);
                if (check_level(p, i, n)) {
                    buf[a.first]++;
                    continue;
                }

                int ax, ay;
                int cx, cy, cz;
                polar_to_axial (p, i, &ax, &ay);
                axial_to_cube  (ax, ay, &cx, &cy, &cz);

                ax/=int(power);
                ay/=int(power);
                size_t modx = std::signbit(ax) ? (3 + ax % 3) % 3 : ax % 3;
                size_t mody = std::signbit(ay) ? (3 + ay % 3) % 3 : ay % 3;
                size_t type = (modx+2*mody) % 3;

                switch (type) {
                    case  1: { 
                                cube_to_polar(cx+1*power, cy,   cz-1*power, &p, &i); buf[polar_to_vortex(p, i)]++;
                                cube_to_polar(cx-1*power, cy+1*power, cz,   &p, &i); buf[polar_to_vortex(p, i)]++;
                                cube_to_polar(cx,   cy-1*power, cz+1*power, &p, &i); buf[polar_to_vortex(p, i)]++;
                                break; // one master at top right
                             }
                    case  2: { 
                                cube_to_polar(cx,   cy+1*power, cz-1*power, &p, &i); buf[polar_to_vortex(p, i)]++;
                                cube_to_polar(cx-1*power, cy,   cz+1*power, &p, &i); buf[polar_to_vortex(p, i)]++;
                                cube_to_polar(cx+1*power, cy-1*power, cz,   &p, &i); buf[polar_to_vortex(p, i)]++;
                                break; // one master at bottom left
                             }
                    default: buf[a.first]++; break; // we are master
                }            
            }
            masters.swap(buf);
            continue;
        }

        power = std::pow(3, (n-1)/2);
        for (const auto & a : masters) {
            vortex_to_polar(a.first, &p, &i);
            if (check_level(p, i, n)) buf[a.first]++;
            else {
                int ax, ay;
                int cx, cy, cz;
                polar_to_axial (p, i, &ax, &ay);
                axial_to_cube  (ax, ay, &cx, &cy, &cz);

                ax/=int(power);
                ay/=int(power);
                if ( (std::signbit(ay) && (std::abs(ay)) % 3 == 2) || (!std::signbit(ay) && ay % 3 == 1) ) {
                    // point down
                    cube_to_polar(cx+2*power, cy-1*power, cz-1*power, &p, &i); buf[polar_to_vortex(p, i)]++;
                    cube_to_polar(cx-1*power, cy+2*power, cz-1*power, &p, &i); buf[polar_to_vortex(p, i)]++;
                    cube_to_polar(cx-1*power, cy-1*power, cz+2*power, &p, &i); buf[polar_to_vortex(p, i)]++;
                }
                else {
                    // point up
                    cube_to_polar(cx+1*power, cy+1*power, cz-2*power, &p, &i); buf[polar_to_vortex(p, i)]++;
                    cube_to_polar(cx-2*power, cy+1*power, cz+1*power, &p, &i); buf[polar_to_vortex(p, i)]++;
                    cube_to_polar(cx+1*power, cy-2*power, cz+1*power, &p, &i); buf[polar_to_vortex(p, i)]++;
                }
            }
        }
        masters.swap(buf);
    }

    std::set<size_t> weights;
    std::map<size_t, size_t> distances;
    size_t max_d = 0;
    for (const auto & a : masters) {
        int ax, ay;
        int cx, cy, cz;
        vortex_to_polar(a.first, &p, &i);
        polar_to_axial (p, i, &ax, &ay);
        axial_to_cube  (ax, ay, &cx, &cy, &cz);
        size_t d = cube_distance(cx, cy, cz, cube_x, cube_y, cube_z);
        distances[a.first] = d;
        if (d > max_d) max_d = d;
    }

    for (const auto & a : distances) {
        weights.insert(max_d - a.second + 1);
    }

    while (weights.size() > 1) {
        size_t w = *weights.begin();
        weights.erase(w);
    }

    size_t added=0;
    for (const auto & a : distances) {
        size_t w = max_d - a.second + 1;
        if (weights.count(w) > 0) {
            (*to)[a.first] += w;
            if (++added >= 1) break;
        }
    }

    /*for (const auto & a : masters) {
        (*to)[a.first]+=a.second;
    }*/
/*
    std::set<size_t> weights;
    for (const auto & a : masters) weights.insert(a.second);

    while (weights.size() > 4) {
        size_t w = *weights.begin();
        weights.erase(w);
    }

    size_t added=0;
    for (const auto & a : masters) {
        if (weights.count(a.second) > 0) {
            (*to)[a.first] += a.second;
            if (++added == 3) break;
        }
    }
*/
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

    gen.seed(vortex);
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
            if (highlight) {
                vtx.color = highlight7 ? al_map_rgb_f(0.0, 0.8, 0.5) :
                            highlight6 ? al_map_rgb_f(1.0, 1.0, 1.0) :
                            highlight5 ? al_map_rgb_f(0.0, 1.0, 1.0) :
                            highlight4 ? al_map_rgb_f(1.0, 0.0, 1.0) :
                            highlight3 ? al_map_rgb_f(1.0, 1.0, 0.0) :
                            highlight2 ? al_map_rgb_f(0.0, 0.0, 1.0) :
                            highlight1 ? al_map_rgb_f(1.0, 0.0, 0.0) : al_map_rgb_f(0.0, (t+1)/6.0, 0.0);
            }
            else vtx.color = al_map_rgb_f(height, height, height);
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
        default : *x = 0;                   *y = 0;                     break;
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

inline size_t HEXAGON::polar_to_vortex(size_t peel, size_t index) {
    if (peel == 0) return 0;
    size_t prev_peel = peel-1;
    size_t vortex_size = 3*prev_peel*prev_peel+3*prev_peel+1;
    size_t index_on_rim = index % (peel*6);
    return vortex_size + index_on_rim;
}

inline void HEXAGON::vortex_to_polar(size_t vortex, size_t * peel, size_t * index) {
    if (vortex == 0) {
        *peel  = 0;
        *index = 0;
        return;
    }
    const double a = 3.0;
    const double b = 3.0;
    double c = 1.0 - vortex;
    // We solve the quadratic equation f(x) = 3x^2 + 3x + 1 here.
    double temp = -0.5 * (b + (std::signbit(b) ? -1.0 : 1.0) * sqrt(b*b - 4.0*a*c));
    double x1 = temp / a;
    double x2 = c / temp;

    size_t prev_peel = std::max(x1, x2);
    size_t first_index = polar_to_vortex(prev_peel+1, 0);
    *peel = prev_peel + 1;
    *index = vortex - first_index;
}

inline size_t HEXAGON::cube_distance(int x1, int y1, int z1, int x2, int y2, int z2) {
    return std::max(std::abs(x1 - x2), std::max(std::abs(y1 - y2), std::abs(z1 - z2)));
}

#endif

