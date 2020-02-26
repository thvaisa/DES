#ifndef CMAKE_PARTICLESYSTEM_H
#define CMAKE_PARTICLESYSTEM_H

#include "common.hpp"
#include <vector>
#include <iostream>

struct Vector{
    double x;
    double y;
    double z;
    Vector(double* xyz){
        x = xyz[0];
        y = xyz[1];
        z = xyz[2];
    }

    Vector(const Vector& obj){
        x = obj.x;
        y = obj.y;
        z = obj.z;
    }
};

struct ParticleData{
    Snapshot* origin;
    double R[6];
    double t0;
    double M0;
    double a;
    double e;

    public:
        ParticleData(){};
        ~ParticleData(){};
        ParticleData(const ParticleData &obj){
            origin = obj.origin;
            std::memcpy(R,obj.R,sizeof(R));
            t0 = obj.t0;
            M0 = obj.M0;
            a = obj.a;
            e = obj.e;
        }

        ParticleData& operator=(const ParticleData& obj) {
            origin = obj.origin;
            std::memcpy(R,obj.R,sizeof(R));
            t0 = obj.t0;
            M0 = obj.M0;
            a = obj.a;
            e = obj.e;
            return *this;
        }

};



class ParticleSystem{
    public:
        ParticleSystem() : snapshots(std::vector<Snapshot>()),
                            positions(std::vector<Vector>()){}
        ~ParticleSystem(){}
        
        Snapshot* create_and_get_snapshot();



        //Elliptical velocity vector and distance to the central body
        void fill_particle_parameters(ParticleData* particle,
                                        Snapshot* snapshot,double* v,
                                        double GM, double mu);

        void generate_dust_velocity(double* arr,double* v_c, 
                                    double* v, double mu, double coeff);



        void evaluate_snapshot(int indx,int NParts,
                                double GM, double mu, 
                                double coeff);


        void evaluate_positions(double t,double maxDist2,
                                double* cPos, double GM, 
                                double mu);



        

        bool empty(){
            return (snapshots.size()==0);
        }

        Snapshot* get_last_snapshot(){
            return &(snapshots[snapshots.size()-1]);
        }


        int snapshot_count(){
            return snapshots.size();
        }


        bool valid_time(int i,double endTime){
            return snapshots[i].t < endTime;
        }

        std::vector<Vector>::iterator getIterator(){
            return positions.begin();
        }

        std::vector<Vector>::iterator endIteration(){
            return positions.end();
        }


    private:
        std::vector<ParticleData> particles;
        std::vector<Snapshot> snapshots;
        std::vector<Vector> positions;
};




#endif