#include "particlesystem.hpp"

Snapshot* ParticleSystem::create_and_get_snapshot(){
    snapshots.push_back(Snapshot());
    return &(snapshots.back());
}


/*
void ParticleSystem::set_radiation_factor(double mu){
    this->mu = mu;
}

void ParticleSystem::set_central_body_gravity(double GM){
    this->GM = GM;
}
*/

//Elliptical velocity vector and distance to the central body
void ParticleSystem::fill_particle_parameters(ParticleData* particle,
                                Snapshot* snapshot,double* v,
                                double GM, double mu){
    double r_c = snapshot->r_c;

    double si = v[2]/std::sqrt(std::pow(v[1],2)+std::pow(v[2],2));
    double ci = v[1]/std::sqrt(std::pow(v[1],2)+std::pow(v[2],2));

    double a = (std::pow(v[0],2)+std::pow(v[1],2)+std::pow(v[2],2))/(mu*GM)-2/r_c;
    double b = std::pow(r_c,2)*(std::pow(v[1],2)+std::pow(v[2],2))/(mu*GM);
    
    double e_d = mu/std::abs(mu)*std::sqrt(1+a*b);

    double q_d = b/(1+e_d);

    double c = std::sqrt(q_d*(1+e_d)/(std::pow(e_d,2)*mu*GM));

    double sinw_d = v[0]*c;
    double cosw_d = c*std::sqrt(std::pow(v[1],2)+std::pow(v[2],2))-1/e_d;

    double aa = q_d/(1-e_d);
    double ohm = 0.0;
    double i = std::atan2(si,ci);
    double omega = std::atan2(sinw_d,cosw_d);
    double vv = 2*M_PI-omega;
    //This is bad estimation from the wikipedia that is originally 
    //from Celestial mechanics books from 50s :D
    //double M = vv-2*e_d*std::sin(vv)+
    //            (3/4.0*std::pow(e_d,2)+1/8.0*std::pow(e_d,4))*std::sin(2*vv)-
    //            1/3.0*std::pow(e_d,3)*std::sin(3*vv)+
    //            5/32.0*std::pow(e_d,4)*std::sin(4*vv);


    //Better approximation
    //Clément Gazzino. Dynamics of a Geostationary Satellite. [Research Report] Rapport LAAS n° 17432,LAAS-CNRS. 2017. ￿hal-01644934v2￿
    double M = 2*std::atan(std::sqrt((1-e_d)/(1+e_d))*std::tan(vv*0.5))
                -e_d*std::sin(2*std::atan(std::sqrt((1-e_d)/(1+e_d))*std::tan(vv*0.5)));
    //std::cout << snapshot->t << " " << wrap_value_0_to_2PI(M) << r_c << " " << mu*GM << std::endl;
    particle->M0 = wrap_value_0_to_2PI(M); 
    particle->t0 = snapshot->t;
    particle->a = aa;
    particle->e = e_d;
    particle->origin = snapshot;
    rotation_matrix(particle->R,omega,i,ohm);
}



//Fulle paper
void ParticleSystem::generate_dust_velocity(double* arr,double* v_c, double* v, double mu, double coeff){
    arr[0]=v_c[0]+v[0]*std::pow(1.0-mu,coeff);
    arr[1]=v_c[1]+v[1]*std::pow(1.0-mu,coeff);
    arr[2]=v_c[2]+v[2]*std::pow(1.0-mu,coeff);
}

void ParticleSystem::evaluate_snapshot(int indx,int NParts,double GM, double mu, double coeff){
    Snapshot* snapshot = &(snapshots[indx]);
    double velocity[3];

    for(int i=0;i<NParts;++i){
        double scale  = 1;
        double v[3] = {(-1+2*rng.rand())*scale,(-1+2*rng.rand())*scale,(-1+2*rng.rand())*scale};
        //std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
        generate_dust_velocity(velocity,snapshot->vel, v, mu,coeff);
        particles.push_back(ParticleData());
        fill_particle_parameters(&(particles.back()),snapshot,
                                velocity,GM,mu); 
    }
}

void ParticleSystem::evaluate_positions(double t,double maxDist2,
                        double* cPos, double GM, 
                        double mu){
    int omitted_particles = 0;
    for(int i=0;i<particles.size();++i){
        //std::cout << i << " " << particles.size() << std::endl;
        ParticleData* particle = &(particles[i]);
        if(t<particles[i].t0) return;
        double E = E_at_t(particle->t0,particle->M0,particle->a,
                                particle->e, t, GM*mu);
        //std::cout << E << " " << particle->t0 << " " << particle->M0 << " " 
        //            << particle->a << " " << particle->e << " " << GM*mu << std::endl;
        
        double v = eccentric_to_true_anomaly(E,particle->e);
        double r_c = particle->a*(1.0-particle->e*std::cos(E));

        double o_xyz[3],xyz_r[3],xyz[3];
        solve_elliptical_coords(o_xyz,r_c,v);
        convert_elliptical_to_cartesian(xyz_r,o_xyz,particle->R);
        matvec_multiply(xyz, particle->origin->R_t, xyz_r);
        
        double dist = 0;

        for(int j=1;j<3;++j) dist+=std::pow(xyz[j]-cPos[j],2);
       
        if(dist<maxDist2){
            positions.push_back(Vector(xyz));
        }else{
            ++omitted_particles;
        }
    }
    //std::cout << "omitted particles " << omitted_particles << std::endl;
}

