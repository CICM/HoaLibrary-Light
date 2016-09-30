/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <Hoa.hpp>
#include "catch.hpp"
#include <cfloat>
#include <iostream>

using namespace hoa;
typedef double hoa_float_t;

TEST_CASE("Encoder 2D", "[Encoder] [2D]")
{
    Encoder<Hoa2d, hoa_float_t> encoder(7);
    const float epsilon = FLT_EPSILON;
    hoa_float_t              input(1.);
    std::vector<hoa_float_t> ouputs(15);
    
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
        
        CHECK(fabs(ouputs[0] - 1) < epsilon);
        CHECK(fabs(ouputs[1] - 0.707106812) < epsilon);
        CHECK(fabs(ouputs[2] - 0.707106812) < epsilon);
        CHECK(fabs(ouputs[3] - 1) < epsilon);
        CHECK(fabs(ouputs[4] - 0.) < epsilon);
        CHECK(fabs(ouputs[5] - 0.707106812) < epsilon);
        CHECK(fabs(ouputs[6] + 0.707106812) < epsilon);
        CHECK(fabs(ouputs[7] - 0.) < epsilon);
        CHECK(fabs(ouputs[8] + 1) < epsilon);
        CHECK(fabs(ouputs[9] + 0.707106812) < epsilon);
        CHECK(fabs(ouputs[10] + 0.707106812) < epsilon);
        CHECK(fabs(ouputs[11] + 1) < epsilon);
        CHECK(fabs(ouputs[12] - 0) < epsilon);
        CHECK(fabs(ouputs[13] + 0.707106812) < epsilon);
        CHECK(fabs(ouputs[14] - 0.7071067812) < epsilon);
    }
    
    SECTION("Azimuth 0.5pi")
    {
        encoder.setAzimuth(0.5*HOA_PI);
        encoder.process(&input, ouputs.data());
        
        CHECK(fabs(ouputs[0] - 1.) < epsilon);
        CHECK(fabs(ouputs[1] - 1.) < epsilon);
        CHECK(fabs(ouputs[2] - 0.) < epsilon);
        CHECK(fabs(ouputs[3] - 0.) < epsilon);
        CHECK(fabs(ouputs[4] + 1.) < epsilon);
        CHECK(fabs(ouputs[5] + 1.) < epsilon);
        CHECK(fabs(ouputs[6] + 0.) < epsilon);
        CHECK(fabs(ouputs[7] - 0.) < epsilon);
        CHECK(fabs(ouputs[8] - 1.) < epsilon);
        CHECK(fabs(ouputs[9] - 1.) < epsilon);
        CHECK(fabs(ouputs[10] + 0.) < epsilon);
        CHECK(fabs(ouputs[11] + 0.) < epsilon);
        CHECK(fabs(ouputs[12] + 1.) < epsilon);
        CHECK(fabs(ouputs[13] + 1.) < epsilon);
        CHECK(fabs(ouputs[14] - 0.) < epsilon);
    }
    
    SECTION("Azimuth 0.75pi")
    {
        encoder.setAzimuth(0.75*HOA_PI);
        encoder.process(&input, ouputs.data());
        
        CHECK(fabs(ouputs[0] - 1) < epsilon);
        CHECK(fabs(ouputs[1] - 0.707106812) < epsilon);
        CHECK(fabs(ouputs[2] + 0.707106812) < epsilon);
        CHECK(fabs(ouputs[3] + 1) < epsilon);
        CHECK(fabs(ouputs[4] - 0.) < epsilon);
        CHECK(fabs(ouputs[5] - 0.707106812) < epsilon);
        CHECK(fabs(ouputs[6] - 0.707106812) < epsilon);
        CHECK(fabs(ouputs[7] - 0.) < epsilon);
        CHECK(fabs(ouputs[8] + 1) < epsilon);
        CHECK(fabs(ouputs[9] + 0.707106812) < epsilon);
        CHECK(fabs(ouputs[10] - 0.707106812) < epsilon);
        CHECK(fabs(ouputs[11] - 1) < epsilon);
        CHECK(fabs(ouputs[12] - 0) < epsilon);
        CHECK(fabs(ouputs[13] + 0.707106812) < epsilon);
        CHECK(fabs(ouputs[14] + 0.7071067812) < epsilon);
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
        
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CHECK(sum_of_elems == 8);
    }
    
    SECTION("Radius 0.5")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(0.5);
        encoder.process(&input, ouputs.data());
        
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CHECK(sum_of_elems == 8);
    }
    
    SECTION("Radius 0.75")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(0.75);
        encoder.process(&input, ouputs.data());
        
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CHECK(sum_of_elems == 8);
    }
    
    SECTION("Radius 2.")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(2.);
        encoder.process(&input, ouputs.data());
        
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CHECK(sum_of_elems == 4);
    }
}


TEST_CASE("Encoder 3D", "[Encoder] [3D]")
{
    Encoder<Hoa3d, hoa_float_t> encoder(7);
    const float epsilon = FLT_EPSILON;
    hoa_float_t              input(1.);
    std::vector<hoa_float_t> ouputs(64);

    SECTION("Azimuth 0 Elevation 0")
    {
        encoder.setAzimuth(0.);
        encoder.setElevation(0.);
        encoder.process(&input, ouputs.data());
        std::cout.precision(16);
        std::cout << "\n";
        for(size_t  i = 0, k = 0; i < 7; ++i)
        {
            for(size_t j = 0; j < i * 2 + 1; ++j)
                std::cout << ouputs[k++] << " ";
            std::cout << "\n";
        }
        
        CHECK(fabs(ouputs[0] - 1.) < epsilon);
        
        CHECK(fabs(ouputs[1] - 0.) < epsilon);
        CHECK(fabs(ouputs[2] - 0.) < epsilon);
        CHECK(fabs(ouputs[3] - -1.) < epsilon);
        
        CHECK(fabs(ouputs[4] - 0.) < epsilon);
        CHECK(fabs(ouputs[5] - 0.) < epsilon);
        CHECK(fabs(ouputs[6] - -0.5) < epsilon);
        CHECK(fabs(ouputs[7] - 0.) < epsilon);
        CHECK(fabs(ouputs[8] - 0.8660254037844386) < epsilon);
        
        CHECK(fabs(ouputs[9] - 0.) < epsilon);
        CHECK(fabs(ouputs[10] - 0.) < epsilon);
        CHECK(fabs(ouputs[11] - 0.) < epsilon);
        CHECK(fabs(ouputs[12] - 0.) < epsilon);
        CHECK(fabs(ouputs[13] - 0.6123724356957945) < epsilon);
        CHECK(fabs(ouputs[14] - 0.) < epsilon);
        CHECK(fabs(ouputs[15] - 0.) < epsilon); // 0.790569
        
        CHECK(fabs(ouputs[26] - 0.) < epsilon);
        CHECK(fabs(ouputs[27] - 0.) < epsilon);
        CHECK(fabs(ouputs[28] - 0.) < epsilon);
        CHECK(fabs(ouputs[29] - 0.) < epsilon);
        CHECK(fabs(ouputs[30] - 0) < epsilon);
        CHECK(fabs(ouputs[31] - -0.4841229182759271) < epsilon);
        CHECK(fabs(ouputs[32] - 0.) < epsilon);
        CHECK(fabs(ouputs[33] - 0) < epsilon);
        CHECK(fabs(ouputs[34] - 0.) < epsilon);
        CHECK(fabs(ouputs[35] - 0.) < epsilon);
        CHECK(fabs(ouputs[36] - 0.) < epsilon); // -0.701561
    }
}

