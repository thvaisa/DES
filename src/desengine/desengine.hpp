
#ifndef CMAKE_DESENGINE_H
#define CMAKE_DESENGINE_H

#include "orbit.hpp"
#include "particlesystem.hpp"
#include "common.hpp"
#include <iostream>

class DESEngine{
    public:

        //Constructor
        DESEngine(){
            //Create orbit engine (Handle comet)
            orbitEngine = OrbitEngine();
            //Create particle engine (Handle particles)
            particleSystem = ParticleSystem();
        };

        //Destructor
        ~DESEngine(){};

        //Create comet 67P. Parameters from wikipedia  
        void initialize(){
            double a = 518060000000.0;          //semi-major axis
            double e = 0.64102;                 //eccentricity
            double omega = in_rads(12.780);     //argument of perapsis
            double ohm = in_rads(50.147);       //longitude of ascending node
            double i = in_rads(7.0405);         //inclination
            double M0 = in_rads(313.71);        //mean anomaly
            double t0 = 2456879.5;              //at t0

            //We add the object to the orbit engine
            orbitEngine.add_object(a,e,omega,ohm,i,M0,t0);
        }

        //Generate nSteps snapshots
        void generate_snapshots(int nSteps,
                                    double startTime,
                                    double endTime){
            //Store the endTime to the DESengine
            maxTime = endTime;
            //Get the step size
            double stepSize = (startTime-endTime)/nSteps;

            //Evalue all steps
            for(int i=0;i<nSteps;++i){

                //Time which should be evaluated
                double time = startTime+stepSize*i;
                if(i==nSteps-1) time=endTime;
                
                //Evaluate the objects (index=0) position and store the snapshot 
                orbitEngine.get_velocity_and_transform(
                                particleSystem.create_and_get_snapshot(),0,time);

            }
        }


        //Initialize particles
        void initialize_particles(  double endTime,
                                    int nPartsPerStep,
                                    double mu,
                                    double GM){
            //Don't try to create particles from snapshot
            //that does not exist
            if(maxTime<endTime){
                std::cout << "Error: snapshot does not exist." << std::endl;
                return;
            }


            //std::cout << particleSystem.empty() << std::endl;
            //In case of multiple initializations, use the older state of 
            //the engine if possible. The particles are evaluated from the
            //snapshots and if the end time is shorter than the last evaluated
            //snapshot, the data could be reused
            //if(!particleSystem.empty()){
            //    if(particleSystem.valid_time(particleSystem.snapshot_count()-1,endTime)) return;
            //}
            
            //Evaluate all snapshots
            for(int i=0; i<particleSystem.snapshot_count();++i){
                //if(!particleSystem.valid_time(i,endTime)) break; //snapshot was created after the given time
                ///evaluate snapshot
                particleSystem.evaluate_snapshot(i,nPartsPerStep,GM,mu);
            }

        }
            
        void print_comet_position(double endTime){
            double cPos[3];
            orbitEngine.get_cartesian_position(cPos, 0, endTime);
            std::cout << cPos[0] << " " << cPos[1] << " " << cPos[2] << std::endl;
        }


        //Evaluate all the particles from the snapshots
        void evaluate_particles(double endTime,
                                    double mu,
                                    double GM,
                                    double maxDistFromComet){
                                        
            //maxium distance from the comet squared
            double maxDist2 = std::pow(maxDistFromComet,2);
            //Get comet position
            double cPos[3];
            orbitEngine.get_cartesian_position(cPos, 0, endTime);
            //Check all the pariticles
            particleSystem.evaluate_positions(endTime,maxDist2,cPos,GM,mu);
        }

        //Print all particle positions at time eTime
        void evaluate_particles(double sTime, double eTime, 
                                double mu, int nPartsPerStep, double maxDistFromComet
                                ){
            initialize_particles(eTime,nPartsPerStep,mu,SOLAR_MASS_PARAMETER);
            evaluate_particles(eTime,mu,SOLAR_MASS_PARAMETER,maxDistFromComet);
            print_particles(eTime);
        }


        //print all particle position to the console
        void print_particles(double endTime){
            double scale = 1.0e-10;
            double cPos[3];
            orbitEngine.get_cartesian_position(cPos, 0, endTime);
            std::cout << cPos[0]*scale << " " << cPos[1]*scale << " " << cPos[2]*scale << std::endl;
            for(auto it=particleSystem.getIterator();
                    it!=particleSystem.endIteration();
                    ++it){
                std::cout << (*it).x*scale<< " " << (*it).y*scale << " " << (*it).z*scale << std::endl;
            }
        }


    private:
        OrbitEngine orbitEngine;
        ParticleSystem particleSystem;
        double maxTime;

};
#endif