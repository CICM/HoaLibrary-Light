/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <limits>
#include <numeric>
#include <cmath>
#include <iostream>

#include <Hoa.hpp>
using namespace hoa;
typedef float hoa_float_t;
template<typename A, typename B> bool hoa_compare(const A val1, const B val2)
{
    return std::abs(static_cast<hoa_float_t>(val1) - static_cast<hoa_float_t>(val2)) < (std::numeric_limits<hoa_float_t>::epsilon() * static_cast<hoa_float_t>(10.));
}

#define CATCH_CONFIG_PREFIX_ALL
#include "catch.hpp"

CATCH_TEST_CASE("Encoder 2D", "[Encoder] [2D]")
{
    Encoder<Hoa2d, hoa_float_t> encoder(7);
    hoa_float_t                 input(1.);
    std::vector<hoa_float_t>    ouputs(15);
    
    CATCH_SECTION("Azimuth 0")
    {
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], 1.));
        CATCH_CHECK(hoa_compare(ouputs[1], 0.));
        CATCH_CHECK(hoa_compare(ouputs[2], 1.));
        CATCH_CHECK(hoa_compare(ouputs[3], 0.));
        CATCH_CHECK(hoa_compare(ouputs[4], 1.));
        CATCH_CHECK(hoa_compare(ouputs[5], 0.));
        CATCH_CHECK(hoa_compare(ouputs[6], 1.));
        CATCH_CHECK(hoa_compare(ouputs[7], 0.));
        CATCH_CHECK(hoa_compare(ouputs[8], 1.));
        CATCH_CHECK(hoa_compare(ouputs[9], 0.));
        CATCH_CHECK(hoa_compare(ouputs[10], 1.));
        CATCH_CHECK(hoa_compare(ouputs[11], 0.));
        CATCH_CHECK(hoa_compare(ouputs[12], 1.));
        CATCH_CHECK(hoa_compare(ouputs[13], 0.));
        CATCH_CHECK(hoa_compare(ouputs[14], 1.));
    }
    
    CATCH_SECTION("Azimuth 0.25pi")
    {
        encoder.setAzimuth(hoa_float_t(0.25*HOA_PI));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], 1.));
        CATCH_CHECK(hoa_compare(ouputs[1], 0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[2], 0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[3], 1.));
        CATCH_CHECK(hoa_compare(ouputs[4], 0.));
        CATCH_CHECK(hoa_compare(ouputs[5], 0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[6], -0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[7], 0.));
        CATCH_CHECK(hoa_compare(ouputs[8], -1.));
        CATCH_CHECK(hoa_compare(ouputs[9], -0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[10], -0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[11], -1.));
        CATCH_CHECK(hoa_compare(ouputs[12], 0.));
        CATCH_CHECK(hoa_compare(ouputs[13], -0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[14], 0.7071067690849304));
    }
    
    CATCH_SECTION("Azimuth 0.5pi")
    {
        encoder.setAzimuth(hoa_float_t(0.5*HOA_PI));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], 1.));
        CATCH_CHECK(hoa_compare(ouputs[1], 1.));
        CATCH_CHECK(hoa_compare(ouputs[2], 0.));
        CATCH_CHECK(hoa_compare(ouputs[3], 0.));
        CATCH_CHECK(hoa_compare(ouputs[4], -1.));
        CATCH_CHECK(hoa_compare(ouputs[5], -1.));
        CATCH_CHECK(hoa_compare(ouputs[6], -0.));
        CATCH_CHECK(hoa_compare(ouputs[7], 0.));
        CATCH_CHECK(hoa_compare(ouputs[8], 1.));
        CATCH_CHECK(hoa_compare(ouputs[9], 1.));
        CATCH_CHECK(hoa_compare(ouputs[10], -0.));
        CATCH_CHECK(hoa_compare(ouputs[11], -0.));
        CATCH_CHECK(hoa_compare(ouputs[12], -1.));
        CATCH_CHECK(hoa_compare(ouputs[13], -1.));
        CATCH_CHECK(hoa_compare(ouputs[14], 0.));
    }
    
    CATCH_SECTION("Azimuth 0.75pi")
    {
        encoder.setAzimuth(hoa_float_t(0.75*HOA_PI));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], 1.));
        CATCH_CHECK(hoa_compare(ouputs[1], 0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[2], -0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[3], -1.));
        CATCH_CHECK(hoa_compare(ouputs[4], 0.));
        CATCH_CHECK(hoa_compare(ouputs[5], 0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[6], 0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[7], 0.));
        CATCH_CHECK(hoa_compare(ouputs[8], -1.));
        CATCH_CHECK(hoa_compare(ouputs[9], -0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[10], 0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[11], 1.));
        CATCH_CHECK(hoa_compare(ouputs[12], 0.));
        CATCH_CHECK(hoa_compare(ouputs[13], -0.7071067690849304));
        CATCH_CHECK(hoa_compare(ouputs[14], -0.7071067690849304));
    }
    
    CATCH_SECTION("Radius 1")
    {
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setRadius(hoa_float_t(1.));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], 1.));
        CATCH_CHECK(hoa_compare(ouputs[1], 0.));
        CATCH_CHECK(hoa_compare(ouputs[2], 1.));
        CATCH_CHECK(hoa_compare(ouputs[3], 0.));
        CATCH_CHECK(hoa_compare(ouputs[4], 1.));
        CATCH_CHECK(hoa_compare(ouputs[5], 0.));
        CATCH_CHECK(hoa_compare(ouputs[6], 1.));
        CATCH_CHECK(hoa_compare(ouputs[7], 0.));
        CATCH_CHECK(hoa_compare(ouputs[8], 1.));
        CATCH_CHECK(hoa_compare(ouputs[9], 0.));
        CATCH_CHECK(hoa_compare(ouputs[10], 1.));
        CATCH_CHECK(hoa_compare(ouputs[11], 0.));
        CATCH_CHECK(hoa_compare(ouputs[12], 1.));
        CATCH_CHECK(hoa_compare(ouputs[13], 0.));
        CATCH_CHECK(hoa_compare(ouputs[14], 1.));
    }
    
    CATCH_SECTION("Radius 0")
    {
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setRadius(hoa_float_t(0.));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], 8));
        CATCH_CHECK(hoa_compare(ouputs[1], 0.));
        CATCH_CHECK(hoa_compare(ouputs[2], 0.));
        CATCH_CHECK(hoa_compare(ouputs[3], 0.));
        CATCH_CHECK(hoa_compare(ouputs[4], 0));
        CATCH_CHECK(hoa_compare(ouputs[5], 0.));
        CATCH_CHECK(hoa_compare(ouputs[6], 0.));
        CATCH_CHECK(hoa_compare(ouputs[7], 0.));
        CATCH_CHECK(hoa_compare(ouputs[8], 0.));
        CATCH_CHECK(hoa_compare(ouputs[9], 0.));
        CATCH_CHECK(hoa_compare(ouputs[10], 0.));
        CATCH_CHECK(hoa_compare(ouputs[11], 0.));
        CATCH_CHECK(hoa_compare(ouputs[12], 0.));
        CATCH_CHECK(hoa_compare(ouputs[13], 0.));
        CATCH_CHECK(hoa_compare(ouputs[14], 0.));
    }
    
    CATCH_SECTION("Radius 0.25")
    {
        encoder.setAzimuth(hoa_float_t(0.25));
        encoder.setRadius(hoa_float_t(0.));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(std::accumulate(ouputs.begin(), ouputs.end(), hoa_float_t(0)), 8.));
    }
    
    CATCH_SECTION("Radius 0.5")
    {
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setRadius(hoa_float_t(0.5));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(std::accumulate(ouputs.begin(), ouputs.end(), hoa_float_t(0)), 8.));
    }
    
    CATCH_SECTION("Radius 0.75")
    {
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setRadius(hoa_float_t(0.75));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(std::accumulate(ouputs.begin(), ouputs.end(), hoa_float_t(0)), 8.));
    }
    
    CATCH_SECTION("Radius 2.")
    {
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setRadius(hoa_float_t(2.));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(std::accumulate(ouputs.begin(), ouputs.end(), hoa_float_t(0)), 4.));
    }
}


CATCH_TEST_CASE("Encoder 3D", "[Encoder] [3D]")
{
    Encoder<Hoa3d, hoa_float_t> encoder(5);
    hoa_float_t                 input(1.);
    std::vector<hoa_float_t>    ouputs(36);
    
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
    
    CATCH_SECTION("Azimuth 0 Elevation 90")
    {
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setElevation(hoa_float_t(HOA_PI2));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], 1.));
        
        CATCH_CHECK(hoa_compare(ouputs[1], 0.));
        CATCH_CHECK(hoa_compare(ouputs[2], 1.));
        CATCH_CHECK(hoa_compare(ouputs[3], 0.));
        
        CATCH_CHECK(hoa_compare(ouputs[4], 0.));
        CATCH_CHECK(hoa_compare(ouputs[5], 0.));
        CATCH_CHECK(hoa_compare(ouputs[6], 1.));
        CATCH_CHECK(hoa_compare(ouputs[7], 0.));
        CATCH_CHECK(hoa_compare(ouputs[8], 0.));
        
        CATCH_CHECK(hoa_compare(ouputs[9], 0.));
        CATCH_CHECK(hoa_compare(ouputs[10], 0.));
        CATCH_CHECK(hoa_compare(ouputs[11], 0.));
        CATCH_CHECK(hoa_compare(ouputs[12], 1.));
        CATCH_CHECK(hoa_compare(ouputs[13], 0.));
        CATCH_CHECK(hoa_compare(ouputs[14], 0.));
        CATCH_CHECK(hoa_compare(ouputs[15], 0.));
        
        CATCH_CHECK(hoa_compare(ouputs[16], 0.));
        CATCH_CHECK(hoa_compare(ouputs[17], 0.));
        CATCH_CHECK(hoa_compare(ouputs[18], 0.));
        CATCH_CHECK(hoa_compare(ouputs[19], 0.));
        CATCH_CHECK(hoa_compare(ouputs[20], 1.));
        CATCH_CHECK(hoa_compare(ouputs[21], 0.));
        CATCH_CHECK(hoa_compare(ouputs[22], 0.));
        CATCH_CHECK(hoa_compare(ouputs[23], 0.));
        CATCH_CHECK(hoa_compare(ouputs[24], 0.));
        
        CATCH_CHECK(hoa_compare(ouputs[25], 0.));
        CATCH_CHECK(hoa_compare(ouputs[26], 0.));
        CATCH_CHECK(hoa_compare(ouputs[27], 0.));
        CATCH_CHECK(hoa_compare(ouputs[28], 0.));
        CATCH_CHECK(hoa_compare(ouputs[29], 0.));
        CATCH_CHECK(hoa_compare(ouputs[30], 1.));
        CATCH_CHECK(hoa_compare(ouputs[31], 0.));
        CATCH_CHECK(hoa_compare(ouputs[32], 0.));
        CATCH_CHECK(hoa_compare(ouputs[33], 0.));
        CATCH_CHECK(hoa_compare(ouputs[34], 0.));
        CATCH_CHECK(hoa_compare(ouputs[35], 0.));
    }
    
    CATCH_SECTION("Azimuth 0 Elevation -19.5")
    {
        encoder.setCoordinates(hoa_float_t(1.), hoa_float_t(0.), hoa_float_t(-19.5/180.*HOA_PI));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], 1.));
        
        CATCH_CHECK(hoa_compare(ouputs[1], 0.));
        CATCH_CHECK(hoa_compare(ouputs[2], -0.333807));
        CATCH_CHECK(hoa_compare(ouputs[3], 0.942642));
        
        CATCH_CHECK(hoa_compare(ouputs[4], 0.));
        CATCH_CHECK(hoa_compare(ouputs[5], 0.));
        CATCH_CHECK(hoa_compare(ouputs[6], -0.33286));
        CATCH_CHECK(hoa_compare(ouputs[7], -0.545007));
        CATCH_CHECK(hoa_compare(ouputs[8], 0.769527));
        
        CATCH_CHECK(hoa_compare(ouputs[9], 0.));
        CATCH_CHECK(hoa_compare(ouputs[10], 0.));
        CATCH_CHECK(hoa_compare(ouputs[11], 0.));
        CATCH_CHECK(hoa_compare(ouputs[12], 0.407722));
        CATCH_CHECK(hoa_compare(ouputs[13], -0.255643));
        CATCH_CHECK(hoa_compare(ouputs[14], -0.574386));
        CATCH_CHECK(hoa_compare(ouputs[15], 0.662186));
        
        CATCH_CHECK(hoa_compare(ouputs[16], 0.));
        CATCH_CHECK(hoa_compare(ouputs[17], 0.));
        CATCH_CHECK(hoa_compare(ouputs[18], 0.));
        CATCH_CHECK(hoa_compare(ouputs[19], 0.));
        CATCH_CHECK(hoa_compare(ouputs[20], 0.011469));
        CATCH_CHECK(hoa_compare(ouputs[21], 0.552251));
        CATCH_CHECK(hoa_compare(ouputs[22], -0.109286));
        CATCH_CHECK(hoa_compare(ouputs[23], -0.584822));
        CATCH_CHECK(hoa_compare(ouputs[24], 0.583889));
        
        CATCH_CHECK(hoa_compare(ouputs[25], 0.));
        CATCH_CHECK(hoa_compare(ouputs[26], 0.));
        CATCH_CHECK(hoa_compare(ouputs[27], 0.));
        CATCH_CHECK(hoa_compare(ouputs[28], 0.));
        CATCH_CHECK(hoa_compare(ouputs[29], 0.));
        CATCH_CHECK(hoa_compare(ouputs[30], -0.333069));
        CATCH_CHECK(hoa_compare(ouputs[31], -0.13656));
        CATCH_CHECK(hoa_compare(ouputs[32], 0.505841));
        CATCH_CHECK(hoa_compare(ouputs[33], 0.001245));
        CATCH_CHECK(hoa_compare(ouputs[34], -0.584718));
        CATCH_CHECK(hoa_compare(ouputs[35], 0.522153));
    }
    
    CATCH_SECTION("Azimuth 120 Elevation -19.5")
    {
        encoder.setCoordinates(hoa_float_t(1.), hoa_float_t(120./180.*HOA_PI), hoa_float_t(-19.5/180.*HOA_PI));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], 1.));
        
        CATCH_CHECK(hoa_compare(ouputs[1], 0.816351));
        CATCH_CHECK(hoa_compare(ouputs[2], -0.333807));
        CATCH_CHECK(hoa_compare(ouputs[3], -0.471321));
        
        CATCH_CHECK(hoa_compare(ouputs[4], -0.66643));
        CATCH_CHECK(hoa_compare(ouputs[5], -0.47199));
        CATCH_CHECK(hoa_compare(ouputs[6], -0.33286));
        CATCH_CHECK(hoa_compare(ouputs[7], 0.272504));
        CATCH_CHECK(hoa_compare(ouputs[8], -0.384763));
        
        CATCH_CHECK(hoa_compare(ouputs[9], 0.));
        CATCH_CHECK(hoa_compare(ouputs[10], 0.497433));
        CATCH_CHECK(hoa_compare(ouputs[11], -0.221393));
        CATCH_CHECK(hoa_compare(ouputs[12], 0.407722));
        CATCH_CHECK(hoa_compare(ouputs[13], 0.127822));
        CATCH_CHECK(hoa_compare(ouputs[14], 0.287193));
        CATCH_CHECK(hoa_compare(ouputs[15], 0.662186));
        
        CATCH_CHECK(hoa_compare(ouputs[16], 0.505663));
        CATCH_CHECK(hoa_compare(ouputs[17], 0.));
        CATCH_CHECK(hoa_compare(ouputs[18], 0.094644));
        CATCH_CHECK(hoa_compare(ouputs[19], 0.478264));
        CATCH_CHECK(hoa_compare(ouputs[20], 0.011469));
        CATCH_CHECK(hoa_compare(ouputs[21], -0.276126));
        CATCH_CHECK(hoa_compare(ouputs[22], 0.054643));
        CATCH_CHECK(hoa_compare(ouputs[23], -0.584822));
        CATCH_CHECK(hoa_compare(ouputs[24], -0.291945));
        
        CATCH_CHECK(hoa_compare(ouputs[25], -0.452198));
        CATCH_CHECK(hoa_compare(ouputs[26], -0.506381));
        CATCH_CHECK(hoa_compare(ouputs[27], 0.));
        CATCH_CHECK(hoa_compare(ouputs[28], -0.438071));
        CATCH_CHECK(hoa_compare(ouputs[29], -0.118265));
        CATCH_CHECK(hoa_compare(ouputs[30], -0.333069));
        CATCH_CHECK(hoa_compare(ouputs[31], 0.06828));
        CATCH_CHECK(hoa_compare(ouputs[32], -0.25292));
        CATCH_CHECK(hoa_compare(ouputs[33], 0.001245));
        CATCH_CHECK(hoa_compare(ouputs[34], 0.29236));
        CATCH_CHECK(hoa_compare(ouputs[35], -0.261076));
    }
    
    CATCH_SECTION("Azimuth 240 Elevation -19.5")
    {
        encoder.setCoordinates(hoa_float_t(1.), hoa_float_t(240./180.*HOA_PI), hoa_float_t(-19.5/180.*HOA_PI));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], 1.));
        
        CATCH_CHECK(hoa_compare(ouputs[1], -0.816352));
        CATCH_CHECK(hoa_compare(ouputs[2], -0.333807));
        CATCH_CHECK(hoa_compare(ouputs[3], -0.471321));
        
        CATCH_CHECK(hoa_compare(ouputs[4], 0.66643));
        CATCH_CHECK(hoa_compare(ouputs[5], 0.47199));
        CATCH_CHECK(hoa_compare(ouputs[6], -0.33286));
        CATCH_CHECK(hoa_compare(ouputs[7], 0.272504));
        CATCH_CHECK(hoa_compare(ouputs[8], -0.384763));
        
        CATCH_CHECK(hoa_compare(ouputs[9], 0.));
        CATCH_CHECK(hoa_compare(ouputs[10], -0.497433));
        CATCH_CHECK(hoa_compare(ouputs[11], 0.221393));
        CATCH_CHECK(hoa_compare(ouputs[12], 0.407722));
        CATCH_CHECK(hoa_compare(ouputs[13], 0.127822));
        CATCH_CHECK(hoa_compare(ouputs[14], 0.287193));
        CATCH_CHECK(hoa_compare(ouputs[15], 0.662186));
        
        CATCH_CHECK(hoa_compare(ouputs[16], -0.505663));
        CATCH_CHECK(hoa_compare(ouputs[17], 0.));
        CATCH_CHECK(hoa_compare(ouputs[18], -0.094644));
        CATCH_CHECK(hoa_compare(ouputs[19], -0.478264));
        CATCH_CHECK(hoa_compare(ouputs[20], 0.011469));
        CATCH_CHECK(hoa_compare(ouputs[21], -0.276126));
        CATCH_CHECK(hoa_compare(ouputs[22], 0.054643));
        CATCH_CHECK(hoa_compare(ouputs[23], -0.584822));
        CATCH_CHECK(hoa_compare(ouputs[24], -0.291945));
        
        CATCH_CHECK(hoa_compare(ouputs[25], 0.452198));
        CATCH_CHECK(hoa_compare(ouputs[26], 0.506381));
        CATCH_CHECK(hoa_compare(ouputs[27], 0.));
        CATCH_CHECK(hoa_compare(ouputs[28], 0.438071));
        CATCH_CHECK(hoa_compare(ouputs[29], 0.118265));
        CATCH_CHECK(hoa_compare(ouputs[30], -0.333069));
        CATCH_CHECK(hoa_compare(ouputs[31], 0.06828));
        CATCH_CHECK(hoa_compare(ouputs[32], -0.25292));
        CATCH_CHECK(hoa_compare(ouputs[33], 0.001245));
        CATCH_CHECK(hoa_compare(ouputs[34], 0.29236));
        CATCH_CHECK(hoa_compare(ouputs[35], -0.261076));
    }
    
    CATCH_SECTION("Radius 0")
    {
        encoder.setCoordinates(hoa_float_t(0.), hoa_float_t(0.), hoa_float_t(0.));
        encoder.process(&input, ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], 6.));
        
        CATCH_CHECK(hoa_compare(ouputs[1], 0.));
        CATCH_CHECK(hoa_compare(ouputs[2], 0.));
        CATCH_CHECK(hoa_compare(ouputs[3], 0.));
        
        CATCH_CHECK(hoa_compare(ouputs[4], 0.));
        CATCH_CHECK(hoa_compare(ouputs[5], 0.));
        CATCH_CHECK(hoa_compare(ouputs[6], 0.));
        CATCH_CHECK(hoa_compare(ouputs[7], 0.));
        CATCH_CHECK(hoa_compare(ouputs[8], 0.));
        
        CATCH_CHECK(hoa_compare(ouputs[9], 0.));
        CATCH_CHECK(hoa_compare(ouputs[10], 0.));
        CATCH_CHECK(hoa_compare(ouputs[11], 0.));
        CATCH_CHECK(hoa_compare(ouputs[12], 0.));
        CATCH_CHECK(hoa_compare(ouputs[13], 0.));
        CATCH_CHECK(hoa_compare(ouputs[14], 0.));
        CATCH_CHECK(hoa_compare(ouputs[15], 0.));
        
        CATCH_CHECK(hoa_compare(ouputs[16], 0.));
        CATCH_CHECK(hoa_compare(ouputs[17], 0.));
        CATCH_CHECK(hoa_compare(ouputs[18], 0.));
        CATCH_CHECK(hoa_compare(ouputs[19], 0.));
        CATCH_CHECK(hoa_compare(ouputs[20], 0.));
        CATCH_CHECK(hoa_compare(ouputs[21], 0.));
        CATCH_CHECK(hoa_compare(ouputs[22], 0.));
        CATCH_CHECK(hoa_compare(ouputs[23], 0.));
        CATCH_CHECK(hoa_compare(ouputs[24], 0.));
        
        CATCH_CHECK(hoa_compare(ouputs[25], 0.));
        CATCH_CHECK(hoa_compare(ouputs[26], 0.));
        CATCH_CHECK(hoa_compare(ouputs[27], 0.));
        CATCH_CHECK(hoa_compare(ouputs[28], 0.));
        CATCH_CHECK(hoa_compare(ouputs[29], 0.));
        CATCH_CHECK(hoa_compare(ouputs[30], 0.));
        CATCH_CHECK(hoa_compare(ouputs[31], 0.));
        CATCH_CHECK(hoa_compare(ouputs[32], 0.));
        CATCH_CHECK(hoa_compare(ouputs[33], 0.));
        CATCH_CHECK(hoa_compare(ouputs[34], 0.));
        CATCH_CHECK(hoa_compare(ouputs[35], 0.));
    }
}

