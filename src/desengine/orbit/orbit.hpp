
#ifndef CMAKE_ORBIT_H
#define CMAKE_ORBIT_H

#include <vector>
#include <cmath>
#include "common.hpp"




struct KeplerianObject{
    double a;       //semi-major axis
    double e;       //eccentricity
    double omega;   //argument of periapsis
    double ohm;     //longitude of ascending node
    double i;       //inclination
    double M0;      //Mean anomaly at time t0
    double t0;      //t0

    KeplerianObject(double a,double e,double omega,double ohm,double i,double M0,double t0){
        this->a = a;
        this->e = e;
        this->omega = omega;
        this->ohm = ohm;
        this->i = i;
        this->M0 = M0;
        this->t0 = t0;
    }

};

struct OrbitData{
    double E;       //eccentric anomaly
    double t;       //time
    double r_c;     //distance from the central body
    double v;       //true anomaly
};

//Be careful with the multiple revolutions, fmod(). Performance hit


//Keplerian orbit
class OrbitEngine{
    OrbitEngine() : objects(std::vector<KeplerianObject>()) {};
    ~OrbitEngine(){};

    void add_object(double a,double e,double omega,
                    double ohm,double i,double M0,double t0);

    //Rene Schwarz, keplerian orbits => cartesian state vectors
    double get_E_at_t(KeplerianObject* obj, double t);

    void get_elliptical_parameters(double* o_xyz, double* o_vel, double* r_c,
                                    KeplerianObject* obj, double t);

    void get_cartesian_position(double* xyz, int indx, double t);
    void convert_to_cartesian_coords(double* xyz, double* v_xyz, 
                                    double* o_xyz, double* o_vel);

    //solve position at time t for particle[indx]
    void get_velocity_and_transform(Snapshot* snapshot, int indx,double t);



    void rotation_matrix(double* R,KeplerianObject* obj);
    void convert_elliptical_to_cartesian(double* xyz,double* o_xyz,
                                        double* R);
    void solve_elliptical_coords(double* arr,double r_c,double v);
    void solve_elliptical_velocity(double* arr,double r_c,
                                    double mu, double a, 
                                    double e,double E);

    KeplerianObject* get_object(int indx);


    private:
        std::vector<KeplerianObject> objects;


        
};

/*

//true anomaly
        
        //Store orbit data
        OrbitData oData = OrbitData();
        oData.E = E;
        oData.t = t;
        oData.r_c = r_c;
        oData.o_xyz = {r_c*(cos),r_c*(),0.0};

*/

#endif