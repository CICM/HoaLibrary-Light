/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <Hoa.hpp>
#include "catch.hpp"
#include <cfloat>
#include <cmath>
//#include <iostream>

using namespace hoa;

TEST_CASE("Encoder 2D", "[Encoder] [2D]")
{
    Encoder<Hoa2d, float> encoder(7);
    const float epsilon = FLT_EPSILON * 10.f;
    float              input(1.);
    std::vector<float> ouputs(15);
    
    SECTION("Azimuth 0")
    {
        encoder.setAzimuth(0.);
        encoder.process(&input, ouputs.data());
        
        CHECK(ouputs[0] == 1.);
        CHECK(ouputs[1] == 0.);
        CHECK(ouputs[2] == 1.);
        CHECK(ouputs[3] == 0.);
        CHECK(ouputs[4] == 1.);
        CHECK(ouputs[5] == 0.);
        CHECK(ouputs[6] == 1.);
        CHECK(ouputs[7] == 0.);
        CHECK(ouputs[8] == 1.);
        CHECK(ouputs[9] == 0.);
        CHECK(ouputs[10] == 1.);
        CHECK(ouputs[11] == 0.);
        CHECK(ouputs[12] == 1.);
        CHECK(ouputs[13] == 0.);
        CHECK(ouputs[14] == 1.);
    }
    
    SECTION("Azimuth 0.25pi")
    {
        encoder.setAzimuth(0.25*HOA_PI);
        encoder.process(&input, ouputs.data());
        
        CHECK(std::abs(ouputs[0] - float(1)) < epsilon);
        CHECK(std::abs(ouputs[1] - float(0.707106812)) < epsilon);
        CHECK(std::abs(ouputs[2] - float(0.707106812)) < epsilon);
        CHECK(std::abs(ouputs[3] - float(1)) < epsilon);
        CHECK(std::abs(ouputs[4] - float(0.)) < epsilon);
        CHECK(std::abs(ouputs[5] - float(0.707106812)) < epsilon);
        CHECK(std::abs(ouputs[6] + float(0.707106812)) < epsilon);
        CHECK(std::abs(ouputs[7] - float(0.)) < epsilon);
        CHECK(std::abs(ouputs[8] + float(1)) < epsilon);
        CHECK(std::abs(ouputs[9] + float(0.707106812)) < epsilon);
        CHECK(std::abs(ouputs[10] + float(0.707106812)) < epsilon);
        CHECK(std::abs(ouputs[11] + float(1)) < epsilon);
        CHECK(std::abs(ouputs[12] - float(0)) < epsilon);
        CHECK(std::abs(ouputs[13] + float(0.707106812)) < epsilon);
        CHECK(std::abs(ouputs[14] - float(0.7071067812)) < epsilon);
    }
    
    SECTION("Azimuth 0.5pi")
    {
        encoder.setAzimuth(0.5*HOA_PI);
        encoder.process(&input, ouputs.data());
        
        CHECK(std::abs(ouputs[0] - 1.) < epsilon);
        CHECK(std::abs(ouputs[1] - 1.) < epsilon);
        CHECK(std::abs(ouputs[2] - 0.) < epsilon);
        CHECK(std::abs(ouputs[3] - 0.) < epsilon);
        CHECK(std::abs(ouputs[4] + 1.) < epsilon);
        CHECK(std::abs(ouputs[5] + 1.) < epsilon);
        CHECK(std::abs(ouputs[6] + 0.) < epsilon);
        CHECK(std::abs(ouputs[7] - 0.) < epsilon);
        CHECK(std::abs(ouputs[8] - 1.) < epsilon);
        CHECK(std::abs(ouputs[9] - 1.) < epsilon);
        CHECK(std::abs(ouputs[10] + 0.) < epsilon);
        CHECK(std::abs(ouputs[11] + 0.) < epsilon);
        CHECK(std::abs(ouputs[12] + 1.) < epsilon);
        CHECK(std::abs(ouputs[13] + 1.) < epsilon);
        CHECK(std::abs(ouputs[14] - 0.) < epsilon);
    }
    
    SECTION("Azimuth 0.75pi")
    {
        encoder.setAzimuth(0.75*HOA_PI);
        encoder.process(&input, ouputs.data());
        
        CHECK(std::abs(ouputs[0] - 1) < epsilon);
        CHECK(std::abs(ouputs[1] - 0.707106812) < epsilon);
        CHECK(std::abs(ouputs[2] + 0.707106812) < epsilon);
        CHECK(std::abs(ouputs[3] + 1) < epsilon);
        CHECK(std::abs(ouputs[4] - 0.) < epsilon);
        CHECK(std::abs(ouputs[5] - 0.707106812) < epsilon);
        CHECK(std::abs(ouputs[6] - 0.707106812) < epsilon);
        CHECK(std::abs(ouputs[7] - 0.) < epsilon);
        CHECK(std::abs(ouputs[8] + 1) < epsilon);
        CHECK(std::abs(ouputs[9] + 0.707106812) < epsilon);
        CHECK(std::abs(ouputs[10] - 0.707106812) < epsilon);
        CHECK(std::abs(ouputs[11] - 1) < epsilon);
        CHECK(std::abs(ouputs[12] - 0) < epsilon);
        CHECK(std::abs(ouputs[13] + 0.707106812) < epsilon);
        CHECK(std::abs(ouputs[14] + 0.7071067812) < epsilon);
    }
    
    SECTION("Radius 1")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(1.);
        encoder.process(&input, ouputs.data());
        
        CHECK(ouputs[0] == 1.);
        CHECK(ouputs[1] == 0.);
        CHECK(ouputs[2] == 1.);
        CHECK(ouputs[3] == 0.);
        CHECK(ouputs[4] == 1.);
        CHECK(ouputs[5] == 0.);
        CHECK(ouputs[6] == 1.);
        CHECK(ouputs[7] == 0.);
        CHECK(ouputs[8] == 1.);
        CHECK(ouputs[9] == 0.);
        CHECK(ouputs[10] == 1.);
        CHECK(ouputs[11] == 0.);
        CHECK(ouputs[12] == 1.);
        CHECK(ouputs[13] == 0.);
        CHECK(ouputs[14] == 1.);
    }
    
    SECTION("Radius 0")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(0.);
        encoder.process(&input, ouputs.data());
        
        CHECK(ouputs[0] == 8);
        CHECK(ouputs[1] == 0.);
        CHECK(ouputs[2] == 0.);
        CHECK(ouputs[3] == 0.);
        CHECK(ouputs[4] == 0);
        CHECK(ouputs[5] == 0.);
        CHECK(ouputs[6] == 0.);
        CHECK(ouputs[7] == 0.);
        CHECK(ouputs[8] == 0.);
        CHECK(ouputs[9] == 0.);
        CHECK(ouputs[10] == 0.);
        CHECK(ouputs[11] == 0.);
        CHECK(ouputs[12] == 0.);
        CHECK(ouputs[13] == 0.);
        CHECK(ouputs[14] == 0.);
    }
    
    SECTION("Radius 0.25")
    {
        encoder.setAzimuth(0.25);
        encoder.setRadius(0.);
        encoder.process(&input, ouputs.data());
        
        float sum_of_elems = 0.;
        for(std::vector<float>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CHECK(sum_of_elems == 8);
    }
    
    SECTION("Radius 0.5")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(0.5);
        encoder.process(&input, ouputs.data());
        
        float sum_of_elems = 0.;
        for(std::vector<float>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CHECK(sum_of_elems == 8);
    }
    
    SECTION("Radius 0.75")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(0.75);
        encoder.process(&input, ouputs.data());
        
        float sum_of_elems = 0.;
        for(std::vector<float>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CHECK(sum_of_elems == 8);
    }
    
    SECTION("Radius 2.")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(2.);
        encoder.process(&input, ouputs.data());
        
        float sum_of_elems = 0.;
        for(std::vector<float>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CHECK(sum_of_elems == 4);
    }
}


TEST_CASE("Encoder 3D", "[Encoder] [3D]")
{
    Encoder<Hoa3d, float> encoder(7);
    const float epsilon = FLT_EPSILON * float(10.f);
    float              input(1.);
    std::vector<float> ouputs(64);

    SECTION("Azimuth 0 Elevation 0")
    {
        encoder.setAzimuth(0.);
        encoder.setElevation(0.);
        encoder.process(&input, ouputs.data());
        
        CHECK(std::abs(ouputs[0] - 1.) < epsilon);
        
        CHECK(std::abs(ouputs[1] - 0.) < epsilon);
        CHECK(std::abs(ouputs[2] - 0.) < epsilon);
        CHECK(std::abs(ouputs[3] - -1.) < epsilon);
        
        CHECK(std::abs(ouputs[4] - 0.) < epsilon);
        CHECK(std::abs(ouputs[5] - 0.) < epsilon);
        CHECK(std::abs(ouputs[6] - -0.5) < epsilon);
        CHECK(std::abs(ouputs[7] - 0.) < epsilon);
        CHECK(std::abs(ouputs[8] - 0.8660254037844386) < epsilon);
        
        CHECK(std::abs(ouputs[9] - 0.) < epsilon);
        CHECK(std::abs(ouputs[10] - 0.) < epsilon);
        CHECK(std::abs(ouputs[11] - 0.) < epsilon);
        CHECK(std::abs(ouputs[12] - 0.) < epsilon);
        CHECK(std::abs(ouputs[13] - 0.6123724356957945) < epsilon);
        CHECK(std::abs(ouputs[14] - 0.) < epsilon);
        CHECK(std::abs(ouputs[15] - -0.7905694150420948) < epsilon);
        
        CHECK(std::abs(ouputs[25] - 0.) < epsilon);
        CHECK(std::abs(ouputs[26] - 0.) < epsilon);
        CHECK(std::abs(ouputs[27] - 0.) < epsilon);
        CHECK(std::abs(ouputs[28] - 0.) < epsilon);
        CHECK(std::abs(ouputs[29] - 0.) < epsilon);
        CHECK(std::abs(ouputs[30] - 0.) < epsilon);
        CHECK(std::abs(ouputs[31] - -0.4841229182759271) < epsilon);
        CHECK(std::abs(ouputs[32] - 0.) < epsilon);
        CHECK(std::abs(ouputs[33] - 0.5229125165837972) < epsilon);
        CHECK(std::abs(ouputs[34] - 0.) < epsilon);
        CHECK(std::abs(ouputs[35] - -0.701560760020114) < epsilon);
        
        CHECK(std::abs(ouputs[49] - 0.) < epsilon);
        CHECK(std::abs(ouputs[50] - 0.) < epsilon);
        CHECK(std::abs(ouputs[51] - 0.) < epsilon);
        CHECK(std::abs(ouputs[52] - 0.) < epsilon);
        CHECK(std::abs(ouputs[53] - 0.) < epsilon);
        CHECK(std::abs(ouputs[54] - 0.) < epsilon);
        CHECK(std::abs(ouputs[55] - 0.) < epsilon);
        CHECK(std::abs(ouputs[56] - 0.) < epsilon);
        CHECK(std::abs(ouputs[57] - 0.4133986423538423) < epsilon);
        CHECK(std::abs(ouputs[58] - 0.) < epsilon);
        CHECK(std::abs(ouputs[59] - -0.42961647140211) < epsilon);
        CHECK(std::abs(ouputs[60] - 0.) < epsilon);
        CHECK(std::abs(ouputs[61] - 0.4749588797990832) < epsilon);
        CHECK(std::abs(ouputs[62] - 0.) < epsilon);
        CHECK(std::abs(ouputs[63] - -0.6472598492877494) < epsilon);
    }
    
    SECTION("Azimuth 0.5pi Elevation 0")
    {
        encoder.setAzimuth(HOA_PI2);
        encoder.setElevation(0.);
        encoder.process(&input, ouputs.data());
        
        CHECK(std::abs(ouputs[0] - 1.) < epsilon);
        
        CHECK(std::abs(ouputs[1] - -1.) < epsilon);
        CHECK(std::abs(ouputs[2] - 0.) < epsilon);
        CHECK(std::abs(ouputs[3] - 0.) < epsilon);
        
        CHECK(std::abs(ouputs[4] - 0.) < epsilon);
        CHECK(std::abs(ouputs[5] - 0.) < epsilon);
        CHECK(std::abs(ouputs[6] - -0.5) < epsilon);
        CHECK(std::abs(ouputs[7] - 0.) < epsilon);
        CHECK(std::abs(ouputs[8] - -0.8660254037844386) < epsilon);
        
        CHECK(std::abs(ouputs[9] - 0.7905694150420948) < epsilon);
        CHECK(std::abs(ouputs[10] - 0.) < epsilon);
        CHECK(std::abs(ouputs[11] - 00.6123724356957945) < epsilon);
        CHECK(std::abs(ouputs[12] - 0.) < epsilon);
        CHECK(std::abs(ouputs[13] - 0.) < epsilon);
        CHECK(std::abs(ouputs[14] - 0.) < epsilon);
        CHECK(std::abs(ouputs[15] - 0.) < epsilon);
        
        CHECK(std::abs(ouputs[25] - -0.701560760020114) < epsilon);
        CHECK(std::abs(ouputs[26] - 0.) < epsilon);
        CHECK(std::abs(ouputs[27] - -0.5229125165837972) < epsilon);
        CHECK(std::abs(ouputs[28] - 0.) < epsilon);
        CHECK(std::abs(ouputs[29] - -0.4841229182759271) < epsilon);
        CHECK(std::abs(ouputs[30] - 0.) < epsilon);
        CHECK(std::abs(ouputs[31] - 0.) < epsilon);
        CHECK(std::abs(ouputs[32] - 0.) < epsilon);
        CHECK(std::abs(ouputs[33] - 0.) < epsilon);
        CHECK(std::abs(ouputs[34] - 0.) < epsilon);
        CHECK(std::abs(ouputs[35] - 0.) < epsilon);
        
        CHECK(std::abs(ouputs[49] - 0.6472598492877494) < epsilon);
        CHECK(std::abs(ouputs[50] - 0.) < epsilon);
        CHECK(std::abs(ouputs[51] - 0.4749588797990832) < epsilon);
        CHECK(std::abs(ouputs[52] - 0.) < epsilon);
        CHECK(std::abs(ouputs[53] - 0.42961647140211) < epsilon);
        CHECK(std::abs(ouputs[54] - 0.) < epsilon);
        CHECK(std::abs(ouputs[55] - 0.4133986423538423) < epsilon);
        CHECK(std::abs(ouputs[56] - 0.) < epsilon);
        CHECK(std::abs(ouputs[57] - 0.) < epsilon);
        CHECK(std::abs(ouputs[58] - 0.) < epsilon);
        CHECK(std::abs(ouputs[59] - 0.) < epsilon);
        CHECK(std::abs(ouputs[60] - 0.) < epsilon);
        CHECK(std::abs(ouputs[61] - 0.) < epsilon);
        CHECK(std::abs(ouputs[62] - 0.) < epsilon);
        CHECK(std::abs(ouputs[63] - 0.) < epsilon);
    }
    
    SECTION("Azimuth 1pi Elevation 0")
    {
        encoder.setAzimuth(HOA_PI);
        encoder.setElevation(0.);
        encoder.process(&input, ouputs.data());
        
        CHECK(std::abs(ouputs[0] - 1.) < epsilon);
        
        CHECK(std::abs(ouputs[1] - 0.) < epsilon);
        CHECK(std::abs(ouputs[2] - 0.) < epsilon);
        CHECK(std::abs(ouputs[3] - 1.) < epsilon);
        
        CHECK(std::abs(ouputs[4] - 0.) < epsilon);
        CHECK(std::abs(ouputs[5] - 0.) < epsilon);
        CHECK(std::abs(ouputs[6] - -0.5) < epsilon);
        CHECK(std::abs(ouputs[7] - 0.) < epsilon);
        CHECK(std::abs(ouputs[8] - 0.8660254037844386) < epsilon);
        
        CHECK(std::abs(ouputs[9] - 0.) < epsilon);
        CHECK(std::abs(ouputs[10] - 0.) < epsilon);
        CHECK(std::abs(ouputs[11] - 0.) < epsilon);
        CHECK(std::abs(ouputs[12] - 0.) < epsilon);
        CHECK(std::abs(ouputs[13] - -0.6123724356957945) < epsilon);
        CHECK(std::abs(ouputs[14] - 0.) < epsilon);
        CHECK(std::abs(ouputs[15] - 0.7905694150420948) < epsilon);
        
        CHECK(std::abs(ouputs[25] - 0.) < epsilon);
        CHECK(std::abs(ouputs[26] - 0.) < epsilon);
        CHECK(std::abs(ouputs[27] - 0.) < epsilon);
        CHECK(std::abs(ouputs[28] - 0.) < epsilon);
        CHECK(std::abs(ouputs[29] - 0.) < epsilon);
        CHECK(std::abs(ouputs[30] - 0.) < epsilon);
        CHECK(std::abs(ouputs[31] - 0.4841229182759271) < epsilon);
        CHECK(std::abs(ouputs[32] - 0.) < epsilon);
        CHECK(std::abs(ouputs[33] - -0.5229125165837972) < epsilon);
        CHECK(std::abs(ouputs[34] - 0.) < epsilon);
        CHECK(std::abs(ouputs[35] - 0.701560760020114) < epsilon);
        
        CHECK(std::abs(ouputs[49] - 0.) < epsilon);
        CHECK(std::abs(ouputs[50] - 0.) < epsilon);
        CHECK(std::abs(ouputs[51] - 0.) < epsilon);
        CHECK(std::abs(ouputs[52] - 0.) < epsilon);
        CHECK(std::abs(ouputs[53] - 0.) < epsilon);
        CHECK(std::abs(ouputs[54] - 0.) < epsilon);
        CHECK(std::abs(ouputs[55] - 0.) < epsilon);
        CHECK(std::abs(ouputs[56] - 0.) < epsilon);
        CHECK(std::abs(ouputs[57] - -0.4133986423538423) < epsilon);
        CHECK(std::abs(ouputs[58] - 0.) < epsilon);
        CHECK(std::abs(ouputs[59] - 0.42961647140211) < epsilon);
        CHECK(std::abs(ouputs[60] - 0.) < epsilon);
        CHECK(std::abs(ouputs[61] - -0.4749588797990832) < epsilon);
        CHECK(std::abs(ouputs[62] - 0.) < epsilon);
        CHECK(std::abs(ouputs[63] - 0.6472598492877494) < epsilon);
    }
    
    SECTION("Azimuth 0 pi Elevation pi/2")
    {
        encoder.setAzimuth(0);
        encoder.setElevation(HOA_PI2);
        encoder.process(&input, ouputs.data());
        
        CHECK(std::abs(ouputs[0] - 1.) < epsilon);
        
        CHECK(std::abs(ouputs[1] - 0.) < epsilon);
        CHECK(std::abs(ouputs[2] - 1.) < epsilon);
        CHECK(std::abs(ouputs[3] - 0.) < epsilon);
        
        CHECK(std::abs(ouputs[4] - 0.) < epsilon);
        CHECK(std::abs(ouputs[5] - 0.) < epsilon);
        CHECK(std::abs(ouputs[6] - 1.) < epsilon);
        CHECK(std::abs(ouputs[7] - 0.) < epsilon);
        CHECK(std::abs(ouputs[8] - 0.) < epsilon);
        
        CHECK(std::abs(ouputs[9] - 0.) < epsilon);
        CHECK(std::abs(ouputs[10] - 0.) < epsilon);
        CHECK(std::abs(ouputs[11] - 0.) < epsilon);
        CHECK(std::abs(ouputs[12] - 1.) < epsilon);
        CHECK(std::abs(ouputs[13] - 0.) < epsilon);
        CHECK(std::abs(ouputs[14] - 0.) < epsilon);
        CHECK(std::abs(ouputs[15] - 0.) < epsilon);
        
        CHECK(std::abs(ouputs[25] - 0.) < epsilon);
        CHECK(std::abs(ouputs[26] - 0.) < epsilon);
        CHECK(std::abs(ouputs[27] - 0.) < epsilon);
        CHECK(std::abs(ouputs[28] - 0.) < epsilon);
        CHECK(std::abs(ouputs[29] - 0.) < epsilon);
        CHECK(std::abs(ouputs[30] - 1.) < epsilon);
        CHECK(std::abs(ouputs[31] - 0.) < epsilon);
        CHECK(std::abs(ouputs[32] - 0.) < epsilon);
        CHECK(std::abs(ouputs[33] - 0.) < epsilon);
        CHECK(std::abs(ouputs[34] - 0.) < epsilon);
        CHECK(std::abs(ouputs[35] - 0.) < epsilon);
        
        CHECK(std::abs(ouputs[49] - 0.) < epsilon);
        CHECK(std::abs(ouputs[50] - 0.) < epsilon);
        CHECK(std::abs(ouputs[51] - 0.) < epsilon);
        CHECK(std::abs(ouputs[52] - 0.) < epsilon);
        CHECK(std::abs(ouputs[53] - 0.) < epsilon);
        CHECK(std::abs(ouputs[54] - 0.) < epsilon);
        CHECK(std::abs(ouputs[55] - 0.) < epsilon);
        CHECK(std::abs(ouputs[56] - 1.) < epsilon);
        CHECK(std::abs(ouputs[57] - 0.) < epsilon);
        CHECK(std::abs(ouputs[58] - 0.) < epsilon);
        CHECK(std::abs(ouputs[59] - 0.) < epsilon);
        CHECK(std::abs(ouputs[60] - 0.) < epsilon);
        CHECK(std::abs(ouputs[61] - 0.) < epsilon);
        CHECK(std::abs(ouputs[62] - 0.) < epsilon);
        CHECK(std::abs(ouputs[63] - 0.) < epsilon);
    }
    
    SECTION("Azimuth 0 pi Elevation pi/4")
    {
        encoder.setAzimuth(HOA_PI4);
        encoder.setElevation(HOA_PI4);
        encoder.process(&input, ouputs.data());
        
        CHECK(std::abs(ouputs[0] - 1.) < epsilon);
        
        CHECK(std::abs(ouputs[1] - -0.5) < epsilon);
        CHECK(std::abs(ouputs[2] - 0.7071067811865475) < epsilon);
        CHECK(std::abs(ouputs[3] - -0.5) < epsilon);
        
        CHECK(std::abs(ouputs[4] - 0.4330127018922194) < epsilon);
        CHECK(std::abs(ouputs[5] - -0.6123724356957945) < epsilon);
        CHECK(std::abs(ouputs[6] - 0.2499999999999999) < epsilon);
        CHECK(std::abs(ouputs[7] - -0.6123724356957946) < epsilon);
        CHECK(std::abs(ouputs[8] - 0.) < epsilon);
        
        CHECK(std::abs(ouputs[9] - -0.1976423537605238) < epsilon);
        CHECK(std::abs(ouputs[10] - 0.6846531968814576) < epsilon);
        CHECK(std::abs(ouputs[11] - -0.4592793267718457) < epsilon);
        CHECK(std::abs(ouputs[12] - -0.176776695296637) < epsilon);
        CHECK(std::abs(ouputs[13] - -0.4592793267718458) < epsilon);
        CHECK(std::abs(ouputs[14] - 0.) < epsilon);
        CHECK(std::abs(ouputs[15] - 0.1976423537605237) < epsilon);
        
        CHECK(std::abs(ouputs[25] - 0.08769509500251425) < epsilon);
        CHECK(std::abs(ouputs[26] - 0.) < epsilon);
        CHECK(std::abs(ouputs[27] - -0.4575484520108226) < epsilon);
        CHECK(std::abs(ouputs[28] - 0.4528555233184198) < epsilon);
        CHECK(std::abs(ouputs[29] - 0.1815460943534729) < epsilon);
        CHECK(std::abs(ouputs[30] - -0.3756504775053533) < epsilon);
        CHECK(std::abs(ouputs[31] - 0.1815460943534729) < epsilon);
        CHECK(std::abs(ouputs[32] - 0.) < epsilon);
        CHECK(std::abs(ouputs[33] - 0.4575484520108223) < epsilon);
        CHECK(std::abs(ouputs[34] - -0.3921843874378481) < epsilon);
        CHECK(std::abs(ouputs[35] - 0.08769509500251432) < epsilon);
        
        CHECK(std::abs(ouputs[49] - 0.04045374058048438) < epsilon);
        CHECK(std::abs(ouputs[50] - -0.2140610743565666) < epsilon);
        CHECK(std::abs(ouputs[51] - 0.3265342298618697) < epsilon);
        CHECK(std::abs(ouputs[52] - 0.) < epsilon);
        CHECK(std::abs(ouputs[53] - -0.2058578925468442) < epsilon);
        CHECK(std::abs(ouputs[54] - -0.3043116672431614) < epsilon);
        CHECK(std::abs(ouputs[55] - 0.3152164647948047) < epsilon);
        CHECK(std::abs(ouputs[56] - 0.1270582497444579) < epsilon);
        CHECK(std::abs(ouputs[57] - 0.3152164647948047) < epsilon);
        CHECK(std::abs(ouputs[58] - 0.) < epsilon);
        CHECK(std::abs(ouputs[59] - 0.2058578925468441) < epsilon);
        CHECK(std::abs(ouputs[60] - -0.5877316282087218) < epsilon);
        CHECK(std::abs(ouputs[61] - 0.3265342298618699) < epsilon);
        CHECK(std::abs(ouputs[62] - 0.) < epsilon);
        CHECK(std::abs(ouputs[63] - -0.04045374058048434) < epsilon);
    }
}

