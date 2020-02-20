
#ifndef CMAKE_ORBIT_H
#define CMAKE_ORBIT_H

#include <vector>
#include <cmath>


inline void normalized(double* arr, double* vec1,double coeff){
    double N = coeff*std::sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2]);
    arr[0] = vec1[0]/N;
    arr[1] = vec1[1]/N;
    arr[2] = vec1[2]/N;
}

inline void cross(double* arr,double* A, double* B) 
{ 
    arr[0] = A[1]*B[2]-A[2]*B[1]; 
    arr[1] = A[0]*B[2]-A[2]*B[0]; 
    arr[2] = A[0]*B[1]-A[1]*B[0]; 
} 

inline void get_transform(double* R, 
                        double* x0, double x1, 
                        double* y0, double y1,
                        double* z0, double z0){
    
    for(i=0;i<3;++i){
        for(j=0;j<3;++j){
            int indx = j+i*3;
            R[indx] = x0[i]*x1[j]+y0[i]*y1[j]+z0[i]*z1[j];
        }    
    }
}

inline void matmul(double* arr, double* M, double* vec){
    for(int i=0;i<3;++i){
        arr[i] = M[i*3+0]*vec[0]+M[i*3+1]*vec[1]+M[i*3+2]*vec[2];
    }
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

    void add_object(double a,double e,double omega,double ohm,double i,double M0,double t0){
        objects.push_back(KeplerianObject(a,e,omega,ohm,i,M0,t0));
    }


    //Rene Schwarz, keplerian orbits => cartesian state vectors
    double get_E_at_t(KeplerianObject* obj, double t){
        //time difference (JD=>s)
        double delta_t = 86400*(t-obj->t0);
        //Mean anomaly at time t
        double M = obj->M0+delta_t*std::sqrt(SOLAR_MASS_PARAMETER/std::pow(obj->a,3));
        M = fmod(2.0*M_PI + fmod(M, 2.0*M_PI), 2.0*M_PI);

        //Solve Kepler
        double Ej = M;
        double Ej1 = M+10;
        while(std::abs(Ej1-Ej)<1.0e-14){
            Ej1=Ej-(Ej-obj->e*std::sin(Ej)-M)/(1.0-obj->e*std::cos(Ej));
        }
        return Ej1;        
    }


    //solve position at time t for particle[indx]
    double velocity_vector(int indx,double t){
        KeplerianObject* obj = &(objects[indx]);
        double E = get_E_at_t(obj, t);
        //True anomaly
        double v = 2.0*std::atan2(std::sqrt(1+obj->e)*std::sin(E*0.5),std::sqrt(1-obj->e)*std::cos(E*0.5));
        //distance to the central body
        double r_c = obj->a*(1.0-obj->e*std::cos(E));

        double o_xyz[3];
        solve_elliptical_coords(o_xyz,r_c,v);

        double o_vel[3];
        solve_elliptical_velocities(o_vel,r_c,SOLAR_MASS_PARAMETER,obj->a,obj->e, E)

        double R[6];
        rotation_matrix(R,obj);

        double xyz[3];
        convert_elliptical_to_cartesian(xyz,o_xyz,R);


        //radial and theta vector in a elliptical plane
        double k_r[3], k_theta[3];
        double k_x[3] = {1,0,0};
        double k_y[3] = {0,1,0};
        double k_z[3] = {0,0,1};
        normalized(k_r,o_xyz,-1);
        cross(k_theta,k_r,k_z)      


        double R_t[3*3], new_vel;
        //get transform matrix 
        get_transform(R_t,k_r,k_x,k_theta,k_y,k_z,k_z);
        matmul(new_vel,R_t,o_vel);

        //double[3] o_vel = solve_elliptical_velocity(obj,E,r_c,SOLAR_MASS_PARAMETER);
    }




    void rotation_matrix(double* R,KeplerianObject* obj){
        double w = obj->omega;
        double i = obj->i;
        double o = obj->ohm;

        double cw=std::cos(w);
        double sw=std::sin(w);

        double co=std::cos(o);
        double so=std::sin(o);

        double ci=std::cos(i);
        double si=std::sin(i);


        R[0] = (cw*co-sw*ci*so);
        R[1] = (sw*co+cw*ci*so);
        R[2] = (cw*so+sw*ci*co);
        R[3] = (cw*ci*co-sw*so);
        R[4] = (sw*si);
        R[5] = (cw*si);



    }



 


    void convert_elliptical_to_cartesian(double* xyz,double* o_xyz,double* R){
        xyz[0] = o_xyz[0]*R[0]-o_xyz[1]*R[1];
        xyz[1] = o_xyz[0]*R[2]+o_xyz[1]*R[3];
        xyz[2] = o_xyz[0]*R[4]+o_xyz[1]*R[5];
    }

    void solve_elliptical_coords(double* arr,double r_c,double v){
        arr[0] = r_c*std::cos(v);
        arr[1] = r_c*std::sin(v);
        arr[2] = 0.0;
    } 


    void solve_elliptical_velocity(double* arr,double r_c,double mu, double a, double e,double E){
        double coeff = std::sqrt(mu*a)/rc;
        
        arr[0] = -coeff*std::sin(E);
        arr[1] = coeff*std::cos(std::sqrt(1.0-std::pow(e,2))*E);
        arr[2] = 0.0;
    } 


    std::vector<KeplerianObject> objects;

    private:
        //scale AU per meter
        const double AU_SCALE = 1.49597870700e+11;  
        
        const double SOLAR_MASS_PARAMETER = 1.32712442099e+20;
        




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