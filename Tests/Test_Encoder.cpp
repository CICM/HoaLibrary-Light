/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <limits>
#include <cmath>
#include <iostream>

#include <Hoa.hpp>
using namespace hoa;
typedef float hoa_float_t;

#define CATCH_CONFIG_PREFIX_ALL
#include "catch.hpp"


CATCH_TEST_CASE("Encoder 2D", "[Encoder] [2D]")
{
    /*
    Encoder<Hoa2d, hoa_float_t> encoder(7);
    const hoa_float_t epsilon = std::numeric_limits<hoa_float_t>::epsilon() * hoa_float_t(10);
    hoa_float_t              input(1.);
    std::vector<hoa_float_t> ouputs(15);
    
    CATCH_SECTION("Azimuth 0")
    {
        encoder.setAzimuth(0.);
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(ouputs[0] == 1.);
        CATCH_CHECK(ouputs[1] == 0.);
        CATCH_CHECK(ouputs[2] == 1.);
        CATCH_CHECK(ouputs[3] == 0.);
        CATCH_CHECK(ouputs[4] == 1.);
        CATCH_CHECK(ouputs[5] == 0.);
        CATCH_CHECK(ouputs[6] == 1.);
        CATCH_CHECK(ouputs[7] == 0.);
        CATCH_CHECK(ouputs[8] == 1.);
        CATCH_CHECK(ouputs[9] == 0.);
        CATCH_CHECK(ouputs[10] == 1.);
        CATCH_CHECK(ouputs[11] == 0.);
        CATCH_CHECK(ouputs[12] == 1.);
        CATCH_CHECK(ouputs[13] == 0.);
        CATCH_CHECK(ouputs[14] == 1.);
    }
    
    CATCH_SECTION("Azimuth 0.25pi")
    {
        encoder.setAzimuth(0.25*HOA_PI);
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(std::abs(ouputs[0]  - hoa_float_t(1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[1]  - hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2]  - hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3]  - hoa_float_t(1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[4]  - hoa_float_t(0.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5]  - hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6]  - hoa_float_t(-0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7]  - hoa_float_t(0.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8]  - hoa_float_t(-1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[9]  - hoa_float_t(-0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - hoa_float_t(-0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - hoa_float_t(-1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - hoa_float_t(0.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - hoa_float_t(-0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - hoa_float_t(0.7071067690849304)) < epsilon);
    }
    
    CATCH_SECTION("Azimuth 0.5pi")
    {
        encoder.setAzimuth(0.5*HOA_PI);
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(std::abs(ouputs[0]  - hoa_float_t(1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[1]  - hoa_float_t(1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2]  - hoa_float_t(0.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3]  - hoa_float_t(0.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[4]  - hoa_float_t(-1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5]  - hoa_float_t(-1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6]  - hoa_float_t(-0.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7]  - hoa_float_t(0.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8]  - hoa_float_t(1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[9]  - hoa_float_t(1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - hoa_float_t(-0.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - hoa_float_t(-0.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - hoa_float_t(-1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - hoa_float_t(-1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - hoa_float_t(0.0000000000000000)) < epsilon);
    }
    
    CATCH_SECTION("Azimuth 0.75pi")
    {
        encoder.setAzimuth(0.75*HOA_PI);
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(std::abs(ouputs[0]  - hoa_float_t(1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[1]  - hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2]  - hoa_float_t(-0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3]  - hoa_float_t(-1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[4]  - hoa_float_t(0.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5]  - hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6]  - hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7]  - hoa_float_t(0.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8]  - hoa_float_t(-1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[9]  - hoa_float_t(-0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - hoa_float_t(1.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - hoa_float_t(0.0000000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - hoa_float_t(-0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - hoa_float_t(-0.7071067690849304)) < epsilon);
    }
    
    CATCH_SECTION("Radius 1")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(1.);
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(ouputs[0] == 1.);
        CATCH_CHECK(ouputs[1] == 0.);
        CATCH_CHECK(ouputs[2] == 1.);
        CATCH_CHECK(ouputs[3] == 0.);
        CATCH_CHECK(ouputs[4] == 1.);
        CATCH_CHECK(ouputs[5] == 0.);
        CATCH_CHECK(ouputs[6] == 1.);
        CATCH_CHECK(ouputs[7] == 0.);
        CATCH_CHECK(ouputs[8] == 1.);
        CATCH_CHECK(ouputs[9] == 0.);
        CATCH_CHECK(ouputs[10] == 1.);
        CATCH_CHECK(ouputs[11] == 0.);
        CATCH_CHECK(ouputs[12] == 1.);
        CATCH_CHECK(ouputs[13] == 0.);
        CATCH_CHECK(ouputs[14] == 1.);
    }
    
    CATCH_SECTION("Radius 0")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(0.);
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(ouputs[0] == 8);
        CATCH_CHECK(ouputs[1] == 0.);
        CATCH_CHECK(ouputs[2] == 0.);
        CATCH_CHECK(ouputs[3] == 0.);
        CATCH_CHECK(ouputs[4] == 0);
        CATCH_CHECK(ouputs[5] == 0.);
        CATCH_CHECK(ouputs[6] == 0.);
        CATCH_CHECK(ouputs[7] == 0.);
        CATCH_CHECK(ouputs[8] == 0.);
        CATCH_CHECK(ouputs[9] == 0.);
        CATCH_CHECK(ouputs[10] == 0.);
        CATCH_CHECK(ouputs[11] == 0.);
        CATCH_CHECK(ouputs[12] == 0.);
        CATCH_CHECK(ouputs[13] == 0.);
        CATCH_CHECK(ouputs[14] == 0.);
    }
    
    CATCH_SECTION("Radius 0.25")
    {
        encoder.setAzimuth(0.25);
        encoder.setRadius(0.);
        encoder.process(&input, ouputs.data());
        
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CATCH_CHECK(std::abs(sum_of_elems  - hoa_float_t(8.)) < epsilon);
    }
    
    CATCH_SECTION("Radius 0.5")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(0.5);
        encoder.process(&input, ouputs.data());
        
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CATCH_CHECK(std::abs(sum_of_elems  - hoa_float_t(8.)) < epsilon);
    }
    
    CATCH_SECTION("Radius 0.75")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(0.75);
        encoder.process(&input, ouputs.data());
        
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CATCH_CHECK(std::abs(sum_of_elems  - hoa_float_t(8.)) < epsilon);
    }
    
    CATCH_SECTION("Radius 2.")
    {
        encoder.setAzimuth(0.);
        encoder.setRadius(2.);
        encoder.process(&input, ouputs.data());
        
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CATCH_CHECK(std::abs(sum_of_elems  - hoa_float_t(4.)) < epsilon);
    }
     */
}

/*
CATCH_TEST_CASE("Encoder 3D", "[Encoder] [3D]")
{
    Encoder<Hoa3d, hoa_float_t> encoder(7);
    const hoa_float_t epsilon = std::numeric_limits<hoa_float_t>::epsilon() * hoa_float_t(10.f);
    hoa_float_t              input(1.);
    std::vector<hoa_float_t> ouputs(64);

    CATCH_SECTION("Azimuth 0 Elevation 0")
    {
        encoder.setAzimuth(0.);
        encoder.setElevation(0.);
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(std::abs(ouputs[0] - 1.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[1] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3] - -1.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[4] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6] - -0.5) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8] - 0.8660254037844386) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[9] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - 0.6123724356957945) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[15] - -0.7905694150420948) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[25] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[26] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[27] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[28] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[29] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[30] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[31] - -0.4841229182759271) < epsilon);
        CATCH_CHECK(std::abs(ouputs[32] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[33] - 0.5229125165837972) < epsilon);
        CATCH_CHECK(std::abs(ouputs[34] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[35] - -0.701560760020114) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[49] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[50] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[51] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[52] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[53] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[54] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[55] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[56] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[57] - 0.4133986423538423) < epsilon);
        CATCH_CHECK(std::abs(ouputs[58] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[59] - -0.42961647140211) < epsilon);
        CATCH_CHECK(std::abs(ouputs[60] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[61] - 0.4749588797990832) < epsilon);
        CATCH_CHECK(std::abs(ouputs[62] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[63] - -0.6472598492877494) < epsilon);
    }
    
    CATCH_SECTION("Azimuth 0.5pi Elevation 0")
    {
        encoder.setAzimuth(HOA_PI2);
        encoder.setElevation(0.);
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(std::abs(ouputs[0] - 1.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[1] - -1.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3] - 0.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[4] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6] - -0.5) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8] - -0.8660254037844386) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[9] - 0.7905694150420948) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - 00.6123724356957945) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[15] - 0.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[25] - -0.701560760020114) < epsilon);
        CATCH_CHECK(std::abs(ouputs[26] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[27] - -0.5229125165837972) < epsilon);
        CATCH_CHECK(std::abs(ouputs[28] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[29] - -0.4841229182759271) < epsilon);
        CATCH_CHECK(std::abs(ouputs[30] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[31] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[32] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[33] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[34] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[35] - 0.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[49] - 0.6472598492877494) < epsilon);
        CATCH_CHECK(std::abs(ouputs[50] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[51] - 0.4749588797990832) < epsilon);
        CATCH_CHECK(std::abs(ouputs[52] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[53] - 0.42961647140211) < epsilon);
        CATCH_CHECK(std::abs(ouputs[54] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[55] - 0.4133986423538423) < epsilon);
        CATCH_CHECK(std::abs(ouputs[56] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[57] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[58] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[59] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[60] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[61] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[62] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[63] - 0.) < epsilon);
    }
    
    CATCH_SECTION("Azimuth 1pi Elevation 0")
    {
        encoder.setAzimuth(HOA_PI);
        encoder.setElevation(0.);
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(std::abs(ouputs[0] - 1.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[1] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3] - 1.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[4] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6] - -0.5) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8] - 0.8660254037844386) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[9] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - -0.6123724356957945) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[15] - 0.7905694150420948) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[25] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[26] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[27] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[28] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[29] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[30] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[31] - 0.4841229182759271) < epsilon);
        CATCH_CHECK(std::abs(ouputs[32] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[33] - -0.5229125165837972) < epsilon);
        CATCH_CHECK(std::abs(ouputs[34] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[35] - 0.701560760020114) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[49] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[50] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[51] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[52] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[53] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[54] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[55] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[56] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[57] - -0.4133986423538423) < epsilon);
        CATCH_CHECK(std::abs(ouputs[58] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[59] - 0.42961647140211) < epsilon);
        CATCH_CHECK(std::abs(ouputs[60] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[61] - -0.4749588797990832) < epsilon);
        CATCH_CHECK(std::abs(ouputs[62] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[63] - 0.6472598492877494) < epsilon);
    }
    
    CATCH_SECTION("Azimuth 0 pi Elevation pi/2")
    {
        encoder.setAzimuth(0);
        encoder.setElevation(HOA_PI2);
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(std::abs(ouputs[0] - 1.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[1] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2] - 1.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3] - 0.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[4] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6] - 1.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8] - 0.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[9] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - 1.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[15] - 0.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[25] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[26] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[27] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[28] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[29] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[30] - 1.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[31] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[32] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[33] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[34] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[35] - 0.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[49] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[50] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[51] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[52] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[53] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[54] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[55] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[56] - 1.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[57] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[58] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[59] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[60] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[61] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[62] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[63] - 0.) < epsilon);
    }
    
    CATCH_SECTION("Azimuth 0 pi Elevation pi/4")
    {
        encoder.setAzimuth(HOA_PI4);
        encoder.setElevation(HOA_PI4);
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(std::abs(ouputs[0] - 1.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[1] - -0.5) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2] - 0.7071067811865475) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3] - -0.5) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[4] - 0.4330127018922194) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5] - -0.6123724356957945) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6] - 0.2499999999999999) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7] - -0.6123724356957946) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8] - 0.) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[9] - -0.1976423537605238) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - 0.6846531968814576) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - -0.4592793267718457) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - -0.176776695296637) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - -0.4592793267718458) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[15] - 0.1976423537605237) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[25] - 0.08769509500251425) < epsilon);
        CATCH_CHECK(std::abs(ouputs[26] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[27] - -0.4575484520108226) < epsilon);
        CATCH_CHECK(std::abs(ouputs[28] - 0.4528555233184198) < epsilon);
        CATCH_CHECK(std::abs(ouputs[29] - 0.1815460943534729) < epsilon);
        CATCH_CHECK(std::abs(ouputs[30] - -0.3756504775053533) < epsilon);
        CATCH_CHECK(std::abs(ouputs[31] - 0.1815460943534729) < epsilon);
        CATCH_CHECK(std::abs(ouputs[32] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[33] - 0.4575484520108223) < epsilon);
        CATCH_CHECK(std::abs(ouputs[34] - -0.3921843874378481) < epsilon);
        CATCH_CHECK(std::abs(ouputs[35] - 0.08769509500251432) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[49] - 0.04045374058048438) < epsilon);
        CATCH_CHECK(std::abs(ouputs[50] - -0.2140610743565666) < epsilon);
        CATCH_CHECK(std::abs(ouputs[51] - 0.3265342298618697) < epsilon);
        CATCH_CHECK(std::abs(ouputs[52] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[53] - -0.2058578925468442) < epsilon);
        CATCH_CHECK(std::abs(ouputs[54] - -0.3043116672431614) < epsilon);
        CATCH_CHECK(std::abs(ouputs[55] - 0.3152164647948047) < epsilon);
        CATCH_CHECK(std::abs(ouputs[56] - 0.1270582497444579) < epsilon);
        CATCH_CHECK(std::abs(ouputs[57] - 0.3152164647948047) < epsilon);
        CATCH_CHECK(std::abs(ouputs[58] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[59] - 0.2058578925468441) < epsilon);
        CATCH_CHECK(std::abs(ouputs[60] - -0.5877316282087218) < epsilon);
        CATCH_CHECK(std::abs(ouputs[61] - 0.3265342298618699) < epsilon);
        CATCH_CHECK(std::abs(ouputs[62] - 0.) < epsilon);
        CATCH_CHECK(std::abs(ouputs[63] - -0.04045374058048434) < epsilon);
    }
}
 */

