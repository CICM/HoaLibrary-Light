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

CATCH_TEST_CASE("Wider 2D", "[Wider] [2D]")
{
    Wider<Hoa2d, hoa_float_t>   wider(7);
    Encoder<Hoa2d, hoa_float_t> encoder(7);
    std::vector<hoa_float_t>    inputs(15);
    std::vector<hoa_float_t>    ouputs(15);
    hoa_float_t                 input(1.);
    
    encoder.setAzimuth(hoa_float_t(0.));
    encoder.process(&input, inputs.data());
    
    CATCH_SECTION("Widening factor 1")
    {
        wider.setWidening(1.);
        wider.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], inputs[0]));
        CATCH_CHECK(hoa_compare(ouputs[1], inputs[1]));
        CATCH_CHECK(hoa_compare(ouputs[2], inputs[2]));
        CATCH_CHECK(hoa_compare(ouputs[3], inputs[3]));
        CATCH_CHECK(hoa_compare(ouputs[4], inputs[4]));
        CATCH_CHECK(hoa_compare(ouputs[5], inputs[5]));
        CATCH_CHECK(hoa_compare(ouputs[6], inputs[6]));
        CATCH_CHECK(hoa_compare(ouputs[7], inputs[7]));
        CATCH_CHECK(hoa_compare(ouputs[8], inputs[8]));
        CATCH_CHECK(hoa_compare(ouputs[9], inputs[9]));
        CATCH_CHECK(hoa_compare(ouputs[10], inputs[10]));
        CATCH_CHECK(hoa_compare(ouputs[11], inputs[11]));
        CATCH_CHECK(hoa_compare(ouputs[12], inputs[12]));
        CATCH_CHECK(hoa_compare(ouputs[13], inputs[13]));
        CATCH_CHECK(hoa_compare(ouputs[14], inputs[14]));
    }
    
    CATCH_SECTION("Widening factor 0")
    {
        wider.setWidening(0.);
        wider.process(inputs.data(), ouputs.data());
        
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
    
    CATCH_SECTION("Widening factor 0.25")
    {
        wider.setWidening(0.25);
        wider.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(hoa_compare(std::accumulate(ouputs.begin(), ouputs.end(), hoa_float_t(0)), 8.));
    }
    
    CATCH_SECTION("Widening factor 0.5")
    {
        wider.setWidening(0.5);
        wider.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(hoa_compare(std::accumulate(ouputs.begin(), ouputs.end(), hoa_float_t(0)), 8.));
    }
    
    CATCH_SECTION("Widening factor 0.75")
    {
        wider.setWidening(0.75);
        wider.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(hoa_compare(std::accumulate(ouputs.begin(), ouputs.end(), hoa_float_t(0)), 8.));
    }
}


CATCH_TEST_CASE("Wider 3D", "[Wider] [3D]")
{
    Wider<Hoa3d, float> wider(7);
    Encoder<Hoa3d, hoa_float_t> encoder(7);

    std::vector<float> inputs(64);
    std::vector<float> ouputs(64);
    hoa_float_t                 input(1.);
    
    
    encoder.setCoordinates(hoa_float_t(1.), hoa_float_t(0.), hoa_float_t(0.));
    encoder.process(&input, inputs.data());
    
    CATCH_SECTION("Widening factor 1")
    {
        wider.setWidening(hoa_float_t(1.));
        wider.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(hoa_compare(ouputs[0], inputs[0]));
        CATCH_CHECK(hoa_compare(ouputs[1], inputs[1]));
        CATCH_CHECK(hoa_compare(ouputs[2], inputs[2]));
        CATCH_CHECK(hoa_compare(ouputs[3], inputs[3]));
        CATCH_CHECK(hoa_compare(ouputs[4], inputs[4]));
        CATCH_CHECK(hoa_compare(ouputs[5], inputs[5]));
        CATCH_CHECK(hoa_compare(ouputs[6], inputs[6]));
        CATCH_CHECK(hoa_compare(ouputs[7], inputs[7]));
        CATCH_CHECK(hoa_compare(ouputs[8], inputs[8]));
        CATCH_CHECK(hoa_compare(ouputs[9], inputs[9]));
        CATCH_CHECK(hoa_compare(ouputs[10], inputs[10]));
        CATCH_CHECK(hoa_compare(ouputs[11], inputs[11]));
        CATCH_CHECK(hoa_compare(ouputs[12], inputs[12]));
        CATCH_CHECK(hoa_compare(ouputs[13], inputs[13]));
        CATCH_CHECK(hoa_compare(ouputs[14], inputs[14]));
        CATCH_CHECK(hoa_compare(ouputs[15], inputs[15]));
        
        CATCH_CHECK(hoa_compare(ouputs[49], inputs[49]));
        CATCH_CHECK(hoa_compare(ouputs[63], inputs[63]));
    }
    
    CATCH_SECTION("Widening factor 0")
    {
        wider.setWidening(hoa_float_t(0.));
        wider.process(inputs.data(), ouputs.data());
        
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
        CATCH_CHECK(hoa_compare(ouputs[15], 0.));
        
        CATCH_CHECK(hoa_compare(ouputs[49], 0.));
        CATCH_CHECK(hoa_compare(ouputs[63], 0.));
    }
    
    
    CATCH_SECTION("Widening factor 0.5")
    {
        wider.setWidening(hoa_float_t(0.5));
        wider.process(inputs.data(), ouputs.data());
        
        //CATCH_CHECK(hoa_compare(std::accumulate(ouputs.begin(), ouputs.end(), hoa_float_t(0)), 8.));
    }
    
    CATCH_SECTION("Widening factor 0.75")
    {
        wider.setWidening(hoa_float_t(0.75));
        wider.process(inputs.data(), ouputs.data());
        
        //CATCH_CHECK(hoa_compare(std::accumulate(ouputs.begin(), ouputs.end(), hoa_float_t(0)), 8.));
    }
}
