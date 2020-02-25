#ifndef CMAKE_COMMON_H
#define CMAKE_COMMON_H

#include <cmath>
#include <iostream>
#include <assert.h> 

//scale AU per meter
const double AU_SCALE = 1.49597870700e+11;  
const double SOLAR_MASS_PARAMETER = 1.32712442099e+20; 


inline double in_rads(double degree){
    return degree/180.0*M_PI;
}


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


struct Snapshot{
    KeplerianObject* obj;   //to which object
    double r_c;             //distance to the central body
    double R_t[3*3];
    double vel[3];          //velocity vector (v_r (comet=>sun),v_theta,v_z (elliptica))
    double t;               //time taken
};


//solve elliptical coords from r_c, v
inline void solve_elliptical_coords(double* arr,double r_c,double v){
    arr[0] = r_c*std::cos(v);
    arr[1] = r_c*std::sin(v);
    arr[2] = 0.0;
} 

//coeff can be used to flip the vector (-1)
inline void normalize(double* arr, double* vec1,double coeff){
    double N = std::sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2]);
    arr[0] = coeff*vec1[0]/N;
    arr[1] = coeff*vec1[1]/N;
    arr[2] = coeff*vec1[2]/N;
}

inline void cross_product(double* arr,double* A, double* B) 
{ 
    arr[0] = A[1]*B[2]-A[2]*B[1]; 
    arr[1] = A[2]*B[0]-A[0]*B[2]; 
    arr[2] = A[0]*B[1]-A[1]*B[0]; 
} 

//Same as the outer product
inline void get_transform(double* R, 
                        double* x0, double* x1, 
                        double* y0, double* y1,
                        double* z0, double* z1){
    
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            int indx = j+i*3;
            R[indx] = x0[i]*x1[j]+y0[i]*y1[j]+z0[i]*z1[j];
        }    
    }
}


//matrix[3*3] multiplication with vector[3]
inline void matvec_multiply(double* arr, double* M, double* vec){
    assert(arr!=vec);
    for(int i=0;i<3;++i){
        arr[i] = M[i*3+0]*vec[0]+M[i*3+1]*vec[1]+M[i*3+2]*vec[2];
    }
}


//Wrap value between 0 to 2Pi
inline double wrap_value_0_to_2PI(double value){
    return fmod(2.0*M_PI + fmod(value, 2.0*M_PI), 2.0*M_PI);
}


//Rotation matrix from elliptical coordinates to the cartesian
static void rotation_matrix(double* R,double w, double i, double o){

    double cw=std::cos(w);
    double sw=std::sin(w);

    double co=std::cos(o);
    double so=std::sin(o);

    double ci=std::cos(i);
    double si=std::sin(i);


    R[0] = (cw*co-sw*ci*so);
    R[1] = -(sw*co+cw*ci*so);
    R[2] = (cw*so+sw*ci*co);
    R[3] = (cw*ci*co-sw*so);
    R[4] = (sw*si);
    R[5] = (cw*si);
}

static double M_from_t(double t, double t0, double M0, double a,  double GM){
    //time difference (JD=>s)
    double delta_t = 86400*(t-t0);
    //Mean anomaly at time t
    double M = M0+delta_t*std::sqrt(GM/std::pow(a,3));
    return wrap_value_0_to_2PI(M);
}

static double E_at_t(double t0,double M0,double a,
                                    double e, double t, const double GM){

    double M = M_from_t(t,t0,M0,a,GM);
    //Solve Kepler
    double Ej = M;
    double Ej1 = M+10;
    double diff_E = Ej1-Ej;


    while(std::abs(diff_E)>1.0e-16){
        Ej1=Ej-(Ej-e*std::sin(Ej)-M)/(1.0-e*std::cos(Ej));
        diff_E = Ej1-Ej;
        Ej = Ej1;
    }
    return Ej1;        
}

inline double eccentric_to_true_anomaly(double E, double e){
    return 2.0*std::atan2(std::sqrt(1+e)*std::sin(E*0.5),std::sqrt(1-e)*std::cos(E*0.5));
}

inline void convert_elliptical_to_cartesian(double* xyz,double* o_xyz,double* R){
    assert(xyz!=o_xyz);
    xyz[0] = o_xyz[0]*R[0]+o_xyz[1]*R[1];
    xyz[1] = o_xyz[0]*R[2]+o_xyz[1]*R[3];
    xyz[2] = o_xyz[0]*R[4]+o_xyz[1]*R[5];
}



#endif

