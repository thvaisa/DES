
#ifndef CMAKE_ORBIT_H
#define CMAKE_ORBIT_H

#include <vector>
#include <cmath>
#include "common.hpp"
#include <iostream>



struct OrbitData{
    double E;       //eccentric anomaly
    double t;       //time
    double r_c;     //distance from the central body
    double v;       //true anomaly
};

//Be careful with the multiple revolutions, fmod(). Performance hit


//Keplerian orbit
class OrbitEngine{
    public:
        OrbitEngine() : objects(std::vector<KeplerianObject>()) {};
        ~OrbitEngine(){};

        void solve_elliptical_velocities(double* arr,double r_c,
                                        double mu, double a, 
                                        double e,double E);


        KeplerianObject* get_object(int indx);
        //Rene Schwarz, keplerian orbits => cartesian state vectors
        double get_E_at_t(KeplerianObject* obj, double t);

        void add_object(double a,double e,double omega,
                        double ohm,double i,double M0,double t0);


        void get_elliptical_parameters(double* o_xyz, double* o_vel, double* r_c,
                                        KeplerianObject* obj, double t);


        void get_cartesian_position(double* xyz, int indx, double t);

        void get_cartesian_position_and_velocity(double* xyz,double* vel, int indx, double t);

        void convert_to_cartesian_coords(double* xyz, double* v_xyz, 
                                    double* o_xyz, double* o_vel, KeplerianObject* obj);

        //solve position at time t for particle[indx]
        void get_velocity_and_transform(Snapshot* snapshot, int indx,double t);

    private:
        std::vector<KeplerianObject> objects;


        
};

#endif