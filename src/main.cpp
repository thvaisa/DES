#include <iostream>
#include <iomanip>
#include "desengine.hpp"
#include "rng.hpp"
int main(int argc, const char *argv[]) {


    std::cout << std::setprecision(16);

    //Starting time Rosetta closest approach 
    double startTime =  2456875.954;
    //Ends +1
    double endTime = startTime+0.5;
    //solar radiation
    double mu = 0.6;
    //maximum distance from the comet that the particles 
    //are not tracked (in meters)
    double maxDistFromComet = 1.0e10;
    //how many particles per time step
    int nPartsPerStep = 100;
    //number of time steps betwen start and end time
    int nSteps = 24*10;
    //Create the dust environment system
    DESEngine engine = DESEngine();
    //Initialize the engine

    engine.initialize();

  

    

    //Eq. 1 (Fulle) 1/k
    double coeff = 1/6.0;

    //Generate snapshots of the comet at the certain intervals
    engine.generate_snapshots(nSteps,startTime,endTime);
    //Create and print particle positions at time endTime
    engine.evaluate_particles(startTime, endTime,mu,nPartsPerStep, 
                                maxDistFromComet,coeff);


    //for(int i=0;i<50;++i){
    //  engine.print_comet_position(startTime+i*0.01);
    //}

    return 0;
}
