//
// Created by Konstantin Gredeskoul on 5/16/17.
//
#include "common.hpp"
#include "gtest/gtest.h"
#include <cmath>

using namespace std;


#define VI vector<long long>
const double abs_error = 0.00000001;

class CommonTest : public ::testing::Test {


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




TEST(CommonTest, solve_elliptical_coords) {
    double r_c = 1.0;
    double v = M_PI/3;
    double arr[3];

    solve_elliptical_coords(arr,r_c,v);
    EXPECT_NEAR(arr[0], r_c*std::cos(v), abs_error); 
    EXPECT_NEAR(arr[1], r_c*std::sin(v), abs_error); 
    EXPECT_EQ(arr[2], 0.0);
}

TEST(CommonTest, normalize) {
    double scale = 1.0;
    double arr[3] = {1,1,1};
    double arr2[3] = {1,1,1};
    normalize(arr,arr,scale);

    double value = std::sqrt(arr[0]*arr[0]+arr[1]*arr[1]+arr[2]*arr[2]);
    EXPECT_NEAR(value,1.0, abs_error); 

    normalize(arr2,arr,-scale);

    EXPECT_NEAR(arr[0], -arr2[0], abs_error); 
    EXPECT_NEAR(arr[1], -arr2[1], abs_error); 
    EXPECT_NEAR(arr[2], -arr2[2], abs_error); 

}

TEST(CommonTest, cross_product) {
    double arr[3] = {3,-3,1};
    double arr2[3] = {4,9,2};
    double arr3[3];
    cross_product(arr3,arr,arr2);
    EXPECT_NEAR(arr3[0], -15.0, abs_error); 
    EXPECT_NEAR(arr3[1], -2.0, abs_error); 
    EXPECT_NEAR(arr3[2], 39.0, abs_error); 
}




TEST(CommonTest, test_transform) {
    double R[3*3];
    double x0[3] = {0.1, 0.2, 0.3};
    double x1[3] = {1,2,3};
    double y0[3] = {0.4 , 0.5 ,0.6};
    double y1[3] = {4,5,6};
    double z0[3] = {0.7 , 0.8 , 0.9};
    double z1[3] = {7,8,9};

    get_transform(R,x0,x1,y0,y1,z0,z1);

    double abs_error = 0.0000001;
    EXPECT_NEAR(R[0], 6.6, abs_error); 
    EXPECT_NEAR(R[1], 7.8, abs_error); 
    EXPECT_NEAR(R[2], 9, abs_error); 
    EXPECT_NEAR(R[3], 7.8, abs_error); 
    EXPECT_NEAR(R[4], 9.3, abs_error); 
    EXPECT_NEAR(R[5], 10.8 , abs_error); 
    EXPECT_NEAR(R[6], 9.0, abs_error); 
    EXPECT_NEAR(R[7], 10.8, abs_error); 
    EXPECT_NEAR(R[8],  12.6, abs_error); 
}

TEST(CommonTest, matvec_multiply) {
    double R[3*3] = {1,2,3,4,5,6,7,8,9};
    double vec[3] = {2,1,3};
    double arr2[3];


    matvec_multiply(arr2,R,vec);
    EXPECT_NEAR(arr2[0], 13, abs_error); 
    EXPECT_NEAR(arr2[1], 31, abs_error);    
    EXPECT_NEAR(arr2[2], 49, abs_error); 
}


TEST(CommonTest, wrap_value_0_to_2PI) {
    EXPECT_NEAR(wrap_value_0_to_2PI(2.0f*M_PI), 0.0, abs_error); 
    EXPECT_NEAR(wrap_value_0_to_2PI(3.0f*M_PI), M_PI, abs_error); 
    EXPECT_NEAR(wrap_value_0_to_2PI(0.0f*M_PI), 0.0, abs_error);
    EXPECT_NEAR(wrap_value_0_to_2PI(M_PI), M_PI, abs_error);
    EXPECT_NEAR(wrap_value_0_to_2PI(-M_PI*0.5f), 3.0/2.0*M_PI, abs_error);
}



TEST(CommonTest, rotation_matrix) {
    double R[6];
    double w = 0.22305307840487532;
    double o = 0.8752302599975964;
    double i = 0.12287990598666076;

    rotation_matrix(R,w,i,o);

    EXPECT_NEAR(R[0], 0.4564058335119641, abs_error); 
    EXPECT_NEAR(R[1], -0.8847822032713326, abs_error); 
    EXPECT_NEAR(R[2],0.8893585124874801, abs_error); 
    EXPECT_NEAR(R[3], 0.4504131193176011, abs_error); 
    EXPECT_NEAR(R[4], 0.02711367555876587, abs_error); 
    EXPECT_NEAR(R[5],  0.11953440810458958, abs_error); 
}


TEST(CommonTest, eccentric_to_true_anomaly){
    double E = in_rads(48.43);
    double e = 0.5;
    EXPECT_NEAR(eccentric_to_true_anomaly(E, e), in_rads(75.84), 0.01); 
}



TEST(CommonTest, M_from_t){
    double t = 2449788.986+390;
    double t0 = 2449788.986+350;
    double M0 = in_rads(55.68836);
    double a = 504594094000;
    double GM = SOLAR_MASS_PARAMETER;

    EXPECT_NEAR(M_from_t(t,t0,M0,a,GM), in_rads(62.052744), 0.00001); 
}



TEST(CommonTest, E_at_t){
    double t = 2449788.986+390;
    double t0 = 2449788.986+350;
    double M0 = in_rads(55.68836);
    double a = 504594094000;
    double e = 0.5;
    double GM = SOLAR_MASS_PARAMETER;

    double E = wrap_value_0_to_2PI(E_at_t(t0,M0,a,e,t,GM));
    EXPECT_NEAR(E,  in_rads(90.69853), 0.00001); 
}

TEST(CommonTest, convert_elliptical_to_cartesian){
  double xyz[3];
  double o_xyz[3] = {1,2,0};
  double R[6] = {1,2,3,4,5,6};
  convert_elliptical_to_cartesian(xyz,o_xyz,R);
  EXPECT_NEAR(xyz[0], 5.0 ,abs_error ); 
  EXPECT_NEAR(xyz[1], 11.0 , abs_error); 
  EXPECT_NEAR(xyz[2], 17.0 , abs_error); 
}
