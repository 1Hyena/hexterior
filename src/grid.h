#include <map>
#include <vector>

#include "peel.h"

class GRID {
    public:
    GRID() {}
   ~GRID() {}

    PEEL &operator[] (size_t i) {
        if (peels.count(i) == 0) peels[i] = PEEL(i);
        return peels[i];
    }

    void fill(size_t peels);
    void to_vertices(std::vector<ALLEGRO_VERTEX> *);

    private:
    std::map<size_t, PEEL> peels;

/*
    void populate(unsigned int layer) {
        size_t count = (layer == 0 ? 1 : layer*6);
        tiles[layer].clear();
        tiles[layer].reserve(count);

        for (size_t i=0; i<count; ++i) tiles[layer].push_back(HEXAGON());
    }

    std::map<unsigned int, std::vector<HEXAGON>> tiles;
*/
};

inline void GRID::fill(size_t peel_count) {
    for (size_t i=0; i<peel_count; ++i) (*this)[i].fill();
}

inline void GRID::to_vertices(std::vector<ALLEGRO_VERTEX> * to) {
    for (const auto & p : peels) p.second.to_vertices(to);
}

