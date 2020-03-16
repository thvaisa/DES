#include <iostream>
#include <iomanip>
#include "desengine.hpp"
#include "rng.hpp"
#include <fstream>


double read_value(char* arg){
    return std::stod(arg); 
}


/*
Arguments: 
start of the simulation in Julian days
end simulation (Julian days)

*/
int main(int argc, const char *argv[]) {

    //set precision
    std::cout << std::setprecision(16);

    //The simulation starts
    double startTime =  2456875.954;
    //The simulation ends
    double endTime = startTime+0.5;
    //solar radiation factor
    double mu = 0.6;
    //maximum distance from the comet that the particles 
    //are not tracked (in meters)
    double maxDistFromComet = 5.0e5;
    //how many particles per time step
    int nPartsPerStep = 100;
    //number of time steps betwen start and end time
    int nSteps = 24*10;
    //Eq. 1 (Fulle) 1/k
    double coeff = 1/6.0;




    //Create comet 67P. Parameters from wikipedia  
    double a = 518060000000.0;          //semi-major axis
    double e = 0.64102;                 //eccentricity
    double omega = in_rads(12.780);     //argument of perapsis
    double ohm = in_rads(50.147);       //longitude of ascending node
    double i = in_rads(7.0405);         //inclination
    double M0 = in_rads(313.71);        //mean anomaly
    double t0 = 2456879.5;              //at t0



    std::vector<double> timeStamps = std::vector<double>();
    //Require that all the input is given
    if(argc>1){
        /*
        a = read_value(argv[1]);    
        e = (int)read_value(argv[2]); 
        omega = (int)read_value(argv[3]); 
        ohm = read_value(argv[4]); 
        i = read_value(argv[5]);    
        M0 = (int)read_value(argv[6]); 
        t0 = (int)read_value(argv[7]);     
    
    
        mu = read_value(argv[8]);    
        maxDistFromComet = read_value(argv[9]);    
        nPartsPerStep = (int)read_value(argv[10]); 
        nSteps = (int)read_value(argv[11]); 
        coeff = read_value(argv[12]); 
        */
        std::ifstream myfile (argv[1]);
        std::string str;
        if(myfile.is_open())
        {
            myfile >> a;
            myfile >> e;
            myfile >> omega;
            myfile >> ohm;
            myfile >> i;
            myfile >> t0;
            myfile >> M0;
            myfile >> mu;
            myfile >> maxDistFromComet;
            myfile >> nPartsPerStep;
            myfile >> nSteps;
            myfile >> coeff;


            double a;
            while (myfile >> a){
                timeStamps.push_back(a);
            }
            myfile.close();
        }
    }else{
        std::cout << "The program will use default arguments" << std::endl;
    }


    std::cout << "a: " << a << std::endl;
    std::cout << "e: " << e << std::endl;
    std::cout << "omega: " << omega << std::endl;
    std::cout << "ohm: " << ohm << std::endl;
    std::cout << "i: " << i << std::endl;
    std::cout << "t0: " << t0 << std::endl;  
    std::cout << "M0: " << M0 << std::endl;   
    

    //print arguments just in case
    //std::cout << "startTime: " << startTime << std::endl; 
    //std::cout << "endTime: " << endTime << std::endl;
    std::cout << "solar radiation factor (mu): " << mu << std::endl;
    std::cout << "maxDistFromComet: " << maxDistFromComet << std::endl;
    std::cout << "nPartsPerStep: " << nPartsPerStep << std::endl;
    std::cout << "nSteps: " << nSteps << std::endl;
    std::cout << "some Fulle coefficient: " << coeff << std::endl;


 
    startTime = timeStamps[0]-1;
    endTime = timeStamps[timeStamps.size()-1];
    std::cout << startTime << " " << endTime << std::endl;

    //Create the dust environment system
    DESEngine engine = DESEngine();
    //Initialize the engine

    engine.initialize(a,e,omega,ohm,i,M0,t0);


    //Generate snapshots of the comet at the certain intervals
    engine.generate_snapshots(nSteps,startTime,endTime);


    for(int i=0;i<timeStamps.size();++i){
        //Create and print particle positions at time endTime
        engine.evaluate_particles(startTime, timeStamps[i], 
                                mu, nPartsPerStep, 
                                maxDistFromComet,coeff);
        std::string fname = std::to_string(i);
        engine.write_to_file(timeStamps[i], 1.0, fname);

    }

    return 0;
}
