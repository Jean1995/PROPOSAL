#include "gtest/gtest.h"
#include "PROPOSAL/Integral.h"
#include "math.h"

double Testfkt(double r){
  return r*r;
}

double Testexp(double r){
    return exp(r);
}

bool relErr(double Is, double HasToBe, double RelError){
    return fabs((Is - HasToBe)/HasToBe)<RelError;
}

TEST(IntegralValue , Zero_to_Three_of_xx ) {
    Integral* Int = new Integral();
    EXPECT_EQ(Int->IntegrateClosed(0,3,Testfkt),9);
    delete Int;
}

TEST(IntegralValue, EqualBorders) {
    Integral* Int = new Integral();

    EXPECT_EQ(Int->IntegrateClosed(3,3,Testfkt),0);

    delete Int;
}

TEST(IntegralValue, SmallError) {
    Integral* Int = new Integral();

    EXPECT_TRUE(  relErr(   Int->IntegrateClosed(0,3,Testexp),exp(3)-1
                            ,1.e-6) //Later GetPrecision();
                );

    delete Int;
}

TEST(IntegralValue, FloatEqual) {
    Integral* Int = new Integral();

    //Last 4 digits can differ. relError<1E-4
    ASSERT_FLOAT_EQ(Int->IntegrateClosed(0,3,Testexp),exp(3)-1);

    delete Int;
}

TEST(IntegralValue, MultiplePrecisions) {
    double xmin=0,xmax=3;
    double  ExactIntegral=exp(3)-1;
    double  CalcIntegral=0;

    double precision = 1E-5;
    for(double precision = 1E-5; precision>1E-16;precision/=10){

        Integral* Int = new Integral(5,20,precision);
        for(int i=1000;i>0;i--)CalcIntegral = i* Int->IntegrateClosed(xmin,xmax,Testexp);

        ASSERT_NEAR(CalcIntegral,ExactIntegral, ExactIntegral*precision);

        delete Int;
    }
}

TEST(IntegralValue, HasToFail) {
    EXPECT_TRUE(false);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}