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
    Encoder<Hoa2d, hoa_float_t> encoder(7);
    const hoa_float_t epsilon = std::numeric_limits<hoa_float_t>::epsilon() * hoa_float_t(10);
    hoa_float_t              input(1.);
    std::vector<hoa_float_t> ouputs(15);
    
    CATCH_SECTION("Azimuth 0")
    {
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(ouputs[0] == hoa_float_t(1.));
        CATCH_CHECK(ouputs[1] == hoa_float_t(0.));
        CATCH_CHECK(ouputs[2] == hoa_float_t(1.));
        CATCH_CHECK(ouputs[3] == hoa_float_t(0.));
        CATCH_CHECK(ouputs[4] == hoa_float_t(1.));
        CATCH_CHECK(ouputs[5] == hoa_float_t(0.));
        CATCH_CHECK(ouputs[6] == hoa_float_t(1.));
        CATCH_CHECK(ouputs[7] == hoa_float_t(0.));
        CATCH_CHECK(ouputs[8] == hoa_float_t(1.));
        CATCH_CHECK(ouputs[9] == hoa_float_t(0.));
        CATCH_CHECK(ouputs[10] == hoa_float_t(1.));
        CATCH_CHECK(ouputs[11] == hoa_float_t(0.));
        CATCH_CHECK(ouputs[12] == hoa_float_t(1.));
        CATCH_CHECK(ouputs[13] == hoa_float_t(0.));
        CATCH_CHECK(ouputs[14] == hoa_float_t(1.));
    }
    
    CATCH_SECTION("Azimuth 0.25pi")
    {
        encoder.setAzimuth(hoa_float_t(0.25*HOA_PI));
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
        encoder.setAzimuth(hoa_float_t(0.5*HOA_PI));
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
        encoder.setAzimuth(hoa_float_t(0.75*HOA_PI));
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
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setRadius(hoa_float_t(1.));
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
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setRadius(hoa_float_t(0.));
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
        encoder.setAzimuth(hoa_float_t(0.25));
        encoder.setRadius(hoa_float_t(0.));
        encoder.process(&input, ouputs.data());
        
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CATCH_CHECK(std::abs(sum_of_elems  - hoa_float_t(8.)) < epsilon);
    }
    
    CATCH_SECTION("Radius 0.5")
    {
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setRadius(hoa_float_t(0.5));
        encoder.process(&input, ouputs.data());
        
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CATCH_CHECK(std::abs(sum_of_elems  - hoa_float_t(8.)) < epsilon);
    }
    
    CATCH_SECTION("Radius 0.75")
    {
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setRadius(hoa_float_t(0.75));
        encoder.process(&input, ouputs.data());
        
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CATCH_CHECK(std::abs(sum_of_elems  - hoa_float_t(8.)) < epsilon);
    }
    
    CATCH_SECTION("Radius 2.")
    {
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setRadius(hoa_float_t(2.));
        encoder.process(&input, ouputs.data());
        
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CATCH_CHECK(std::abs(sum_of_elems  - hoa_float_t(4.)) < epsilon);
    }
}



template<typename A, typename B> bool hoa_compare(const A val1, const B val2)
{
    return std::abs(static_cast<hoa_float_t>(val1) - static_cast<hoa_float_t>(val2)) < (std::numeric_limits<hoa_float_t>::epsilon() * static_cast<hoa_float_t>(10.));
}


CATCH_TEST_CASE("Encoder 3D", "[Encoder] [3D]")
{
    Encoder<Hoa3d, hoa_float_t> encoder(5);
    hoa_float_t                 input(1.);
    std::vector<hoa_float_t>    ouputs(36);
    const hoa_float_t epsilon = std::numeric_limits<hoa_float_t>::epsilon() * static_cast<hoa_float_t>(10.);
    

    CATCH_SECTION("Azimuth 0 Elevation 0")
    {
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setElevation(hoa_float_t(0.));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], 1.));
        
        CATCH_CHECK(hoa_compare(ouputs[1], 0.));
        CATCH_CHECK(hoa_compare(ouputs[2], 0.));
        CATCH_CHECK(hoa_compare(ouputs[3], 1.));
        
        CATCH_CHECK(hoa_compare(ouputs[4], 0.));
        CATCH_CHECK(hoa_compare(ouputs[5], 0.));
        CATCH_CHECK(hoa_compare(ouputs[6], -0.5));
        CATCH_CHECK(hoa_compare(ouputs[7], 0.));
        CATCH_CHECK(hoa_compare(ouputs[8], 0.866025));
        
        CATCH_CHECK(hoa_compare(ouputs[9], 0.));
        CATCH_CHECK(hoa_compare(ouputs[10], 0.));
        CATCH_CHECK(hoa_compare(ouputs[11], 0.));
        CATCH_CHECK(hoa_compare(ouputs[12], 0.));
        CATCH_CHECK(hoa_compare(ouputs[13], -0.612372));
        CATCH_CHECK(hoa_compare(ouputs[14], 0.));
        CATCH_CHECK(hoa_compare(ouputs[15], 0.790569));
        
        CATCH_CHECK(hoa_compare(ouputs[16], 0.));
        CATCH_CHECK(hoa_compare(ouputs[17], 0.));
        CATCH_CHECK(hoa_compare(ouputs[18], 0.));
        CATCH_CHECK(hoa_compare(ouputs[19], 0.));
        CATCH_CHECK(hoa_compare(ouputs[20], 0.375));
        CATCH_CHECK(hoa_compare(ouputs[21], 0.));
        CATCH_CHECK(hoa_compare(ouputs[22], -0.559017));
        CATCH_CHECK(hoa_compare(ouputs[23], 0.));
        CATCH_CHECK(hoa_compare(ouputs[24], 0.73951));
        
        CATCH_CHECK(hoa_compare(ouputs[25], 0.));
        CATCH_CHECK(hoa_compare(ouputs[26], 0.));
        CATCH_CHECK(hoa_compare(ouputs[27], 0.));
        CATCH_CHECK(hoa_compare(ouputs[28], 0.));
        CATCH_CHECK(hoa_compare(ouputs[29], 0.));
        CATCH_CHECK(hoa_compare(ouputs[30], 0.));
        CATCH_CHECK(hoa_compare(ouputs[31], 0.484123));
        CATCH_CHECK(hoa_compare(ouputs[32], 0.));
        CATCH_CHECK(hoa_compare(ouputs[33], -0.522913));
        CATCH_CHECK(hoa_compare(ouputs[34], 0.));
        CATCH_CHECK(hoa_compare(ouputs[35], 0.701561));
    }
    /*
    CATCH_SECTION("Azimuth 0.5pi Elevation 0")
    {
        encoder.setAzimuth(hoa_float_t(HOA_PI2));
        encoder.setElevation(hoa_float_t(0.));
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
    }
    
    CATCH_SECTION("Azimuth 1pi Elevation 0")
    {
        encoder.setAzimuth(hoa_float_t(HOA_PI));
        encoder.setElevation(hoa_float_t(0.));
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
    }
    
    CATCH_SECTION("Azimuth 0 pi Elevation pi/2")
    {
        encoder.setAzimuth(hoa_float_t(0));
        encoder.setElevation(hoa_float_t(HOA_PI2));
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
    }
    
    CATCH_SECTION("Azimuth 0 pi Elevation pi/4")
    {
        encoder.setAzimuth(hoa_float_t(HOA_PI4));
        encoder.setElevation(hoa_float_t(HOA_PI4));
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
    }
     */
}

