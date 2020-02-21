#ifndef CMAKE_COMMON_H
#define CMAKE_COMMON_H


//scale AU per meter
const double AU_SCALE = 1.49597870700e+11;  
const double SOLAR_MASS_PARAMETER = 1.32712442099e+20; 


inline double in_rads(double degree){
    return degree*180.0/M_PI;
}



struct Snapshot{
    KeplerianObject* obj;   //to which object
    double r_c;             //distance to the central body
    double R_t[3*3];
    double vel[3];          //velocity vector (v_r (comet=>sun),v_theta,v_z (elliptica))
    double t;               //time taken
};


inline void normalize(double* arr, double* vec1,double coeff){
    double N = coeff*std::sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2]);
    arr[0] = vec1[0]/N;
    arr[1] = vec1[1]/N;
    arr[2] = vec1[2]/N;
}

inline void cross_product(double* arr,double* A, double* B) 
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

inline void matvec_multiply(double* arr, double* M, double* vec){
    for(int i=0;i<3;++i){
        arr[i] = M[i*3+0]*vec[0]+M[i*3+1]*vec[1]+M[i*3+2]*vec[2];
    }
}

inline double wrap_value_0_to_2PI(double value){
    return fmod(2.0*M_PI + fmod(value, 2.0*M_PI), 2.0*M_PI);
}


//Rotation matrix from elliptical coordinates to the cartesian
void rotation_matrix(double* R,double w, double i, double o){

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

double get_E_at_t(double t0,double M0,double a,
                                    double e, double t, double GM){
    //time difference (JD=>s)
    double delta_t = 86400*(t-t0);
    //Mean anomaly at time t
    double M = M0+delta_t*std::sqrt(GM/std::pow(a,3));
    M = wrap_value_from_0_to_2PI(M);

    //Solve Kepler
    double Ej = M;
    double Ej1 = M+10;
    while(std::abs(Ej1-Ej)<1.0e-14){
        Ej1=Ej-(Ej-e*std::sin(Ej)-M)/(1.0-e*std::cos(Ej));
    }
    return Ej1;        
}

double eccentricity_to_true_anomaly(double E, double e){
    return 2.0*std::atan2(std::sqrt(1+e)*std::sin(E*0.5),std::sqrt(1-e)*std::cos(E*0.5));
}

void convert_elliptical_to_cartesian(double* xyz,double* o_xyz,double* R){
    xyz[0] = o_xyz[0]*R[0]-o_xyz[1]*R[1];
    xyz[1] = o_xyz[0]*R[2]+o_xyz[1]*R[3];
    xyz[2] = o_xyz[0]*R[4]+o_xyz[1]*R[5];
}



#endif

