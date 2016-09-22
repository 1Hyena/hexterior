#ifndef PEEL_H_20_09_2016
#define PEEL_H_20_09_2016

#include "hexagon.h"

class PEEL {
    public:
    PEEL(const PEEL &p) : nr(p.nr) {max_sz = (nr == 0 ? 1 : nr*6);}
    PEEL(size_t nr = 0) : nr(nr)   {max_sz = (nr == 0 ? 1 : nr*6);}
   ~PEEL() {}

    HEXAGON &operator[] (size_t i) {
        size_t index  = i % max_sz;
        if (hexagons.count(index) == 0) hexagons[index] = HEXAGON(nr, index);
        return hexagons[index];
    }

    size_t size() const {return hexagons.size();}

    void fill();
    void to_vertices(std::vector<ALLEGRO_VERTEX> *) const;
    size_t get_nr() const {return nr;}

    private:
    size_t nr;
    size_t max_sz;
    std::map<size_t, HEXAGON> hexagons;   
};

inline void PEEL::fill() {
    for (size_t i=0; i<max_sz; ++i) {
        (*this)[i];
    }
}

inline void PEEL::to_vertices(std::vector<ALLEGRO_VERTEX> * to) const {
    for (const auto & h : hexagons) h.second.to_vertices(to);
}

#endif

