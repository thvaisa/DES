#include <iostream>
#include <iomanip>
#include "desengine.hpp"

int main(int argc, const char *argv[]) {


    std::cout << std::setprecision(16);

    int steps = 1;
    //Starting time Rosetta closest approach 
    double startTime =  2456875.954;
    //Ends +1
    double endTime = startTime+0.0;
    //solar radiation
    double mu = 1.0;
    //maximum distance from the comet that the particles are not tracked (in meters)
    double maxDistFromComet = 20000000000000;
    //how many particles per time step
    int nPartsPerStep = 1;
    //number of time steps betwen start and end time
    int nSteps = 10;
    //Create the dust environment system
    DESEngine engine = DESEngine();
    //Initialize the engine

    engine.initialize();

    //for(int i=0;i<1;++i){
    //  engine.print_comet_position((startTime+i)*1.0f);
    //}


    //Generate snapshots of the comet at the certain intervals
    engine.generate_snapshots(steps,startTime,endTime);
    //Create and print particle positions at time endTime
    engine.evaluate_particles(startTime, endTime,mu,nPartsPerStep, maxDistFromComet);
    return 0;
}
