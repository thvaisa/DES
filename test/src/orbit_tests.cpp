//
// Created by Konstantin Gredeskoul on 5/16/17.
//
#include "common.hpp"
#include "orbit.hpp"
#include "gtest/gtest.h"
#include <cmath>

using namespace std;


#define VI vector<long long>
const double abs_error = 0.00000001;

class OrbitTest : public ::testing::Test {


protected:
  //VI numerators   = {5, 9, 17, 933345453464353416L};
  //VI denominators = {2, 3, 19, 978737423423423499L};
  //VI divisions    = {2, 3, 0, 0};
  //VI remainders   = {1, 0, 17, 933345453464353416};

  virtual void SetUp() {
  };

  virtual void TearDown() {
  };

  virtual void verify(int index) {
    //Fraction       f        = Fraction{numerators.at(index), denominators.at(index)};
    //DivisionResult expected = DivisionResult{divisions.at(index), remainders.at(index)};
    //DivisionResult result   = Division(f).divide();
    //EXPECT_EQ(result.division, expected.division);
   // EXPECT_EQ(result.remainder, expected.remainder);
    EXPECT_EQ(1, 1);
  }
};




TEST(OrbitTest, solve_elliptical_coords) {
    OrbitEngine orbitE = OrbitEngine();
    double arr[3];
    double r_c = 2;
    double mu = 4;
    double a = 3;
    double e = 0.6;
    double E = 0.5;
    orbitE.solve_elliptical_velocities(arr,r_c,mu,a,e,E);
    EXPECT_NEAR(arr[0],-0.830389391308554 , abs_error);
    EXPECT_NEAR(arr[1], 1.21601406802447 , abs_error);
    EXPECT_NEAR(arr[2],0.0 , abs_error); 
}



TEST(OrbitTest, add_object) {
    OrbitEngine orbitE = OrbitEngine();
    double a = 1;
    double e = 2;
    double omega = 3;
    double ohm = 4;
    double i = 5;
    double M0 = 6;
    double t0 = 7;
    orbitE.add_object(a,e,omega,ohm,i,M0,t0);
    KeplerianObject* obj = orbitE.get_object(0);

    EXPECT_EQ(obj->a,a);
    EXPECT_EQ(obj->e,e);
    EXPECT_EQ(obj->omega,omega);
    EXPECT_EQ(obj->ohm,ohm);
    EXPECT_EQ(obj->i,i);
    EXPECT_EQ(obj->M0,M0);
    EXPECT_EQ(obj->t0,t0);
}



TEST(OrbitTest, get_E_at_t) {
    OrbitEngine orbitE = OrbitEngine();
    double t = 2449788.986+390;
    double t0 = 2449788.986+350;
    double M0 = in_rads(55.68836);
    double a = 504594094000;
    double e = 0.5;
    double GM = SOLAR_MASS_PARAMETER;

    double omega = 3;
    double ohm = 4;
    double i = 5;

    KeplerianObject* obj = new KeplerianObject(a,e,omega,ohm,i,M0,t0);



    double E = orbitE.get_E_at_t(obj, t);

    EXPECT_NEAR(E, in_rads(90.69853), 0.00001); 
    delete(obj);
}

TEST(OrbitTest,get_elliptical_parameters){
    OrbitEngine orbitE = OrbitEngine();
    double o_xyz[3], o_vel[3],r_c;

    double a = 518060000000.0;          //semi-major axis
    double e = 0.64102;                 //eccentricity
    double omega = in_rads(12.780);     //argument of perapsis
    double ohm = in_rads(50.147);       //longitude of ascending node
    double i = in_rads(7.0405);         //inclination
    double M0 = in_rads(313.71);        //mean anomaly
    double t0 = 2456879.5;              //at t0

    //We add the object to the orbit engine
    orbitE.add_object(a,e,omega,ohm,i,M0,t0);
    double t = 2456879.5+10;
    orbitE.get_elliptical_parameters(o_xyz,o_vel,&r_c,orbitE.get_object(0),t);


    EXPECT_NEAR(o_xyz[0],-2.51403209e+11 , 100000);
    EXPECT_NEAR(o_xyz[1], -3.92771690e+11 , 100000);
    EXPECT_NEAR(o_xyz[2],0.0 , abs_error); 


    EXPECT_NEAR(o_vel[0],17563.49953143 , 0.1);
    EXPECT_NEAR(o_vel[1], 2125.39833385 , 0.1);
    EXPECT_NEAR(o_vel[2],0.0 , abs_error); 

    EXPECT_NEAR(r_c,466340191062.3299, 1000); 

}

TEST(OrbitTest, get_cartesian_position){
    OrbitEngine orbitE = OrbitEngine();

    double a = 518060000000.0;          //semi-major axis
    double e = 0.64102;                 //eccentricity
    double omega = in_rads(12.780);     //argument of perapsis
    double ohm = in_rads(50.147);       //longitude of ascending node
    double i = in_rads(7.0405);         //inclination
    double M0 = in_rads(313.71);        //mean anomaly
    double t0 = 2456879.5;              //at t0

    //We add the object to the orbit engine
    orbitE.add_object(a,e,omega,ohm,i,M0,t0);
    double t = 2456879.5+10;

    double xyz[3];
    double o_xyz[3] = {1,2,3};
    double v_xyz[3];
    double o_vel[3] = {4,5,6};
    
    orbitE.convert_to_cartesian_coords(xyz, v_xyz,o_xyz,o_vel,orbitE.get_object(0));

    EXPECT_NEAR(xyz[0],-1.31315857 , abs_error);
    EXPECT_NEAR(xyz[1], 1.79018475, abs_error);
    EXPECT_NEAR(xyz[2],0.26618249 , abs_error); 

    EXPECT_NEAR(v_xyz[0],-2.59828768 , abs_error);
    EXPECT_NEAR(v_xyz[1], 5.80949965, abs_error);
    EXPECT_NEAR(v_xyz[2],0.70612674 , abs_error);

}


TEST(OrbitTest,get_velocity_and_transform){
    OrbitEngine orbitE = OrbitEngine();
    std::cout << std::setprecision(16);
    double a = 518060000000.0;          //semi-major axis
    double e = 0.64102;                 //eccentricity
    double omega = in_rads(12.780);     //argument of perapsis
    double ohm = in_rads(50.147);       //longitude of ascending node
    double i = in_rads(7.0405);         //inclination
    double M0 = in_rads(313.71);        //mean anomaly
    double t0 = 2456879.5;              //at t0

    //We add the object to the orbit engine
    orbitE.add_object(a,e,omega,ohm,i,M0,t0);
    double t = 2456879.5+10;

    Snapshot* snapshot = new Snapshot();

    orbitE.get_velocity_and_transform(snapshot, 0,t);


    EXPECT_EQ(snapshot->obj, orbitE.get_object(0));

    double* vel = snapshot->vel;
    double time = snapshot->t;
    double* R_t = snapshot->R_t;

    EXPECT_NEAR(snapshot->r_c, 466340190986.0647, 0.001);
    EXPECT_NEAR(time,t , abs_error); 



    EXPECT_NEAR(vel[0],11258.55455934 , abs_error);
    EXPECT_NEAR(vel[1],13646.93328567, abs_error);
    EXPECT_NEAR(vel[2],0.0 , abs_error);


    EXPECT_NEAR(R_t[0],0.49915387 , abs_error);
    EXPECT_NEAR(R_t[1],0.86138914, abs_error);
    EXPECT_NEAR(R_t[2],-0.09409659 , abs_error);
    EXPECT_NEAR(R_t[3], -0.8588089  , abs_error);
    EXPECT_NEAR(R_t[4],0.50623889, abs_error);
    EXPECT_NEAR(R_t[5],0.0785459 , abs_error);
    EXPECT_NEAR(R_t[6],-0.11529394 , abs_error);
    EXPECT_NEAR(R_t[7],-0.0416045, abs_error);
    EXPECT_NEAR(R_t[8],-0.99245976, abs_error);
    delete(snapshot);
}