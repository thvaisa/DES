#include <iostream>
#include "desengine.hpp"

int main(int argc, const char *argv[]) {
    int steps = 2;
    double startTime = 123.0;             
    double endTime = 1234.0;
    double mu = 1.0;
    double maxDistFromComet = 2000.0f;
    int nPartsPerStep = 
    DESEngine engine = DESEngine();
    engine.initialize();
    engine.generate_snapshots(steps,startTime,endTime);
    engine.evaluate_particles(sTime, eTime,mu,nPartsPerStep, maxDistFromComet;
    return 0;
}
