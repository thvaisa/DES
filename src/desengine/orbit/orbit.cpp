
#include "orbit.hpp"




//solve elliptical velocities from r_c,mu,a,e_E
void OrbitEngine::solve_elliptical_velocities(double* arr,double r_c,
                                                double mu, double a, 
                                                double e,double E){
    double coeff = std::sqrt(mu*a)/r_c;
    
    arr[0] = -coeff*std::sin(E);
    arr[1] = coeff*std::sqrt(1.0-std::pow(e,2))*std::cos(E);
    arr[2] = 0.0;
} 

//Get Keplerian object from the list of stored objects
KeplerianObject* OrbitEngine::get_object(int indx){
    return &(objects[indx]);
}

//Return eccentric anomaly from t for the object
double OrbitEngine::get_E_at_t(KeplerianObject* obj, double t){
    return E_at_t(obj->t0,obj->M0,obj->a,
                        obj->e,t, SOLAR_MASS_PARAMETER);    
}

//Store Keplerian object
void OrbitEngine::add_object(double a,double e,
                            double omega,double ohm,double i,
                            double M0,double t0){
    objects.push_back(KeplerianObject(a,e,omega,ohm,i,M0,t0));
}



//Solve elliptical parameters
void OrbitEngine::get_elliptical_parameters(double* o_xyz, double* o_vel, double* r_c,
                                KeplerianObject* obj, double t){
    double E = get_E_at_t(obj, t);
    
    double v = eccentric_to_true_anomaly(E,obj->e);
    //distance to the central body
    *r_c = obj->a*(1.0-obj->e*std::cos(E));

    solve_elliptical_coords(o_xyz,*r_c,v);
    solve_elliptical_velocities(o_vel,*r_c,SOLAR_MASS_PARAMETER,obj->a,obj->e, E);
}


void OrbitEngine::get_cartesian_position(double* xyz, int indx, double t){
    KeplerianObject* obj = get_object(indx);
    double r_c,o_xyz[3],o_vel[3];
    get_elliptical_parameters(o_xyz,o_vel,&r_c,obj,t);

    double R[6];
    rotation_matrix(R,obj->omega,obj->i,obj->ohm);
    convert_elliptical_to_cartesian(xyz,o_xyz,R);
}

void OrbitEngine::get_cartesian_position_and_velocity(double* xyz,double* vel, int indx, double t){
    KeplerianObject* obj = get_object(indx);
    double r_c,o_xyz[3],o_vel[3];
    get_elliptical_parameters(o_xyz,o_vel,&r_c,obj,t);

    double R[6];
    rotation_matrix(R,obj->omega,obj->i,obj->ohm);
    convert_elliptical_to_cartesian(xyz,o_xyz,R);
    convert_elliptical_to_cartesian(vel,o_vel,R);
}




void OrbitEngine::convert_to_cartesian_coords(double* xyz, double* v_xyz, 
                                double* o_xyz, double* o_vel, KeplerianObject* obj){
    double R[6];
    rotation_matrix(R,obj->omega,obj->i,obj->ohm);
    convert_elliptical_to_cartesian(xyz,o_xyz,R);
    convert_elliptical_to_cartesian(v_xyz,o_vel,R);
}


//solve position at time t for particle[indx]
void OrbitEngine::get_velocity_and_transform(Snapshot* snapshot, int indx,double t){


  

    KeplerianObject* obj = get_object(indx);
    double o_xyz[3],o_vel[3],v_xyz[3],xyz[3], r_c;
    get_elliptical_parameters(o_xyz, o_vel,&r_c, obj,t);

    convert_to_cartesian_coords(xyz, v_xyz, o_xyz, o_vel,obj);

    snapshot->obj = obj;
    snapshot->r_c = r_c;
    snapshot->t = t;
    //radial and theta vector in a elliptical plane
    double k_r[3], k_theta[3];
    double k_x[3] = {1,0,0};
    double k_y[3] = {0,1,0};
    double k_z[3] = {0,0,1};
    normalize(k_r,o_xyz,-1);
    cross_product(k_theta,k_r,k_z); 
   

    double R_t1[3*3], new_vel[3];

    //Particle ellipse to ejector ellipse
    get_transform(R_t1,k_r,k_x,k_theta,k_y,k_z,k_z);



    //Velocity vector(radial,theta,z) (in elliptical coordinates)
    matvec_multiply(snapshot->vel,R_t1,o_vel);
    //double[3] o_vel = solve_elliptical_velocity(obj,E,r_c,SOLAR_MASS_PARAMETER);


    //ejector ellipse to inertial frame
    double k_ellipse_z[3],k_theta2[3];
    for(int i=0;i<3;++i) xyz[i] = -xyz[i];
    cross_product(k_ellipse_z,xyz,v_xyz);
    normalize(k_ellipse_z,k_ellipse_z,1);
    cross_product(k_theta2,k_ellipse_z,xyz);
    normalize(k_theta2,k_theta2,1.0);
    normalize(k_r,xyz,-1);


    /*
    std::cout << k_r[0] << " " << k_r[1] << " " << k_r[2] << std::endl;
    std::cout << k_x[0] << " " << k_x[1] << " " << k_x[2] << std::endl;

    std::cout << k_theta2[0] << " " << k_theta2[1] << " " << k_theta2[2] << std::endl;
    std::cout << k_y[0] << " " << k_y[1] << " " << k_y[2] << std::endl;

    std::cout << k_ellipse_z[0] << " " << k_ellipse_z[1] << " " << k_ellipse_z[2] << std::endl;
    std::cout << k_z[0] << " " << k_z[1] << " " << k_z[2] << std::endl;
    */


    get_transform(snapshot->R_t,k_r,k_x,k_theta2,k_y,k_ellipse_z,k_z);
    
}


