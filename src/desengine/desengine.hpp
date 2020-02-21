
#ifndef CMAKE_DESENGINE_H
#define CMAKE_DESENGINE_H

#include "orbit.hpp"
#include "particlesystem.hpp"
#include "common.hpp"

class DESEngine{

    DESEngine(){
        orbitEngine = OrbitEngine();
        particleSystem = ParticleSystem();
    };
    ~DESEngine(){};

    void initialize(){
        //Create comet 67P
        double a = 518060000000.0;      //semi-major axis
        double e = 0.64102;             //eccentricity
        double omega = in_rads(45);     //argument of perapsis
        double ohm = in_rads(45);       //longiute of asceinding node
        double i = in_rads(47.0405);    //inclination
        double M0 = in_rads(0);         //mean anomaly
        double t0 = 2457247.65907;      //at t0

        orbitEngine.add_object(a,e,omega,ohm,i,M0,t0);
    }

    void generate_snapshots(int steps,
                                double startTime,
                                double endTime){
        maxTime = endTime;
        double stepSize = (startTime-endTime)/steps;

        for(int i=0,i<steps;++i){
            double time = startTime+stepSize*i;
            if(i==steps-1): time=endTime;

            orbitEngine.get_velocity_and_transform(
                            particleSystem.create_and_get_snapshot(),0,time);

        }
    }

    void initialize_particles(  double endTime,
                                int nPartsPerStep,
                                double mu,
                                double GM){
        if(maxTime<endTime){
            std::cout << "Error" << std::endl;
            return;
        }

        if(!particleSystem.empty()){
            if(particleSystem.get_last_snapshot().t<endTime) return;
        }

        for(int i=0; i<particleSystem.snapshot_count();++i){
            if(!particleSystem.valid_time(i,endTime)) continue; //snapshot was created after the given time
            particleSystem.evaluate_snapshot(i,nPartsPerStep,GM,mu);
        }

    }
        
    void evaluate_particles(double endTime,
                                double mu,
                                double GM,
                                double maxDistFromComet){
        if(maxTime<endTime){
            std::cout << "Error" << std::endl;
            return;
        }
        double maxDist = std::pow(maxDistFromComet,2)
        double cPos[3];
        orbitEngine.get_cartesian_position(cPos, 0, endTime);
        evaluate_positions(endTime,maxDist2,cPos,GM,mu);
    }


    void evaluate_particles(double sTime, double eTime, 
                            double mu, int nPartsPerStep, double maxDistFromComet,
                            ){
        initialize_particles(eTime,nPartsPerStep,mu,GM);
        evaluate_particles(eTime,mu,SOLAR_MASS_PARAMETER,maxDistFromComet);
        print_particles();
    }

    void print_particles(){
        for(auto it=particleSystem.getPositionIterator();
                it!=particleSystem.endIteration();
                ++it){
            std::cout << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << std::endl;
        }
    }


    private:
        OrbitEngine orbitEngine;
        ParticleSystem particleSystem;
        double maxTime;






#endif