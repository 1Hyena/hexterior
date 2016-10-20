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
    void generate_heightmap(size_t);

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

inline void GRID::generate_heightmap(size_t /*levels*/) {
    size_t peel_count = peels.size();
    for (size_t i=0; i<peel_count; ++i) {
        size_t sz = (*this)[i].size();
        for (size_t j=0; j<sz; ++j) {
            std::map<size_t, size_t> last_masters;
            (*this)[i][j].to_masters(4, &last_masters);

            size_t weight=0;
            double noise = 0.0;
            for (const auto & a : last_masters) weight += a.second;
            for (const auto & a : last_masters) {
                size_t peel, index;
                HEXAGON::vortex_to_polar(a.first, &peel, &index);
                if (peel >= peel_count) {
                    weight -= a.second;
                    continue;
                }
                noise += (*this)[peel][index].get_noise(1) * a.second/double(weight);
            }

            (*this)[i][j].set_height(noise);
/*
            for (size_t k=2; k<=levels; ++k) {
                std::map<size_t, size_t> masters;
                masters.swap(last_masters);
                

                for (const auto & a : masters) {
                    size_t peel, index;
                    vortex_to_polar(a.first, &peel, &index);
                    
                }
            }
*/
        }
        
    }
}

