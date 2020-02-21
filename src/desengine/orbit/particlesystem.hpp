#ifndef CMAKE_PARTICLESYSTEM_H
#define CMAKE_PARTICLESYSTEM_H

#include "common.hpp"

struct Vector{
    double x;
    double y;
    double z;
}

struct ParticleData{
    Snapshot* origin;
    double[6] R;
    double t0;
    double M0;
    double a;
    double e;
}



class ParticleSystem{

    ParticleSystem() : snapshot(std::vector<Snapshot>(),
                        positions(std::vector<positions>())){}
    ~ParticleSystem(){}
    
    Snapshot* create_and_get_snapshot(){
        return snapshot.push_back(Snapshot()).back();
    }

    set_radiation_factor(double mu){
        this->mu = mu;
    }

    set_central_body_gravity(double GM){
        this->GM = GM;
    }


    //Elliptical velocity vector and distance to the central body
    void fill_particle_parameters(ParticleData* particle,
                                    Snapshot* snapshot,double* v,
                                    double GM, double mu){
        double r_c = snapshot->r_c;

        double si = v[2]/std::sqrt(v[1]**2+v[2]**2)
        double ci = v[1]/std::sqrt(v[1]**2+v[2]**2)

        double a = (v[0]**2+v[1]**2+v[2]**2)/(mu*GM)-2/r_c
        double b = r_c**2*(v[1]**2+v[2]**2)/(mu*GM)
        
        double e_d = mu/std::abs(mu)*std::sqrt(1+a*b)

        double q_d = b/(1+e_d)

        double c = std::sqrt(q_d*(1+e_d)/(e_d**2*mu*GM))

        double sinw_d = v[0]*c
        double cosw_d = c*std::sqrt(v[1]**2+v[2]**2)-1/e_d

        double a = q_d/(1-e_d);
        double ohm = 0.0;
        double i = std::atan2(si,ci);
        double omega = std:.atan2(sinw_d,cosw_d);
        double v = 2*M_PI-omega;
        double M = v-2*e_d*std::sin(v)+
                    (3/4.0*e_d**2+1/8.0*e_d**4)*std::sin(2*v)-
                    1/3.0*e_d**3*std::sin(3*v)+
                    5/32.0*e_d**4*std::sin(4*v);
        particle->M0 = wrap_value_0_to_2PI(M); 
        particle->t0 = snapshot->t0;
        particle->a = a;
        particle->e = e_d;
        particle->origin = snapshot;
        rotation_matrix(particle->R,omega,i,ohm);
    }

    void generate_dust_velocity(double* arr,double v){
        arr[0]=v[0];
        arr[1]=v[1];
        arr[2]=v[2];
        return arr;
    }

    void evaluate_snapshot(int indx,int NParts,double GM, double mu){
        Snapshot* snapshot = snapshots[indx];
        double[3] velocity;
        for(int i=0;i<NParts;++i){
            generate_dust_velocity(velocity,snapshot->vel);
            particles.push_back(ParticleData());
            fill_particle_parameters(&(particles.back()),snapshot,
                                    velocity,GM,mu); 
        }
    }

    void evaluate_positions(double t,double maxDist2,
                            double* cPos, double GM, 
                            double mu){
        int omitted_particles = 0;
        for(int i=0;i<particles.count();++i){
            ParticleData* particle = particles[i];
            if(t>particles[i].t0) return;
            Vector3 xyz = Vector3();
            get_E_at_t(particle->t0,particle->M0,particle->a,
                                    particle->e, t, GM*mu);
            v = eccentricity_to_true_anomaly(E,e);
            double r_c = particle->a*(1.0-particle->e*std::cos(E));

            double o_xyz[3],xyz[3];
            solve_elliptical_coords(o_xyz,r_c,v);
            convert_elliptical_to_cartesian(xyz,o_xyz,particleR);
            matvec_multiply(xyz, particle->snapshot->R_t, xyz);
            
            double dist;
            for(int j=1;j<3;++j) dist+=std::pow(xyz[i]-cPos[i],2)
            if(dist<maxDist2){
                position.push_back(xyz);
            }else{
                ++omitted_particles;
            }
        }
    }


    std::vector<ParticleData> particles;
    std::vector<Vector3> position;


    private:
        std::vector<Snapshot>() snapshots;
        std::vector<Vector>() positions;
        double mu = 1.0;
        double GM = SOLAR_MASS_PARAMETER;
}




#endif