
#ifndef CMAKE_RNG_H
#define CMAKE_RNG_H

#include <random>

//Make this better later
class RNG{
    public:
        RNG(): dist{0.0,1.0}{
            gen = new std::mt19937_64(time(0));
        };
        ~RNG(){
            delete(gen);
        };

        double rand(){
            return dist(*gen);
        }

    private:
        std::uniform_real_distribution<double> dist;
        std::mt19937_64* gen;
};




#endif