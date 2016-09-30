/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <Hoa.hpp>
#include "catch.hpp"
#include <cfloat>

using namespace hoa;

TEST_CASE("MultiEncoder 2D", "[Encoder] [2D]")
{
    MultiEncoder<Hoa2d, float> encoder(7, 3);
    const float epsilon = FLT_EPSILON;
    std::vector<float> inputs(5);
    std::vector<float> ouputs(15);
    
    encoder.setRadius(0, 0.);
    encoder.setAzimuth(0, 0.);
    encoder.setElevation(0, 0.);
    encoder.setMute(0, false);
    
    encoder.setRadius(1, 1.);
    encoder.setAzimuth(1, 1.);
    encoder.setElevation(1, 1.);
    encoder.setMute(1, false);
    
    encoder.setRadius(2, 2.);
    encoder.setAzimuth(2, 2.);
    encoder.setElevation(2, 2.);
    encoder.setMute(2, false);
    
    CHECK(encoder.getRadius(0) == 0.);
    CHECK(encoder.getAzimuth(0) == 0.);
    CHECK(encoder.getElevation(0) == 0.);
    CHECK(encoder.getMute(0) == false);
    
    CHECK(encoder.getRadius(1) == 1.);
    CHECK(encoder.getAzimuth(1) == 1.);
    CHECK(encoder.getElevation(1) == 1.);
    CHECK(encoder.getMute(1) == false);
    
    CHECK(encoder.getRadius(2) == 2.);
    CHECK(encoder.getAzimuth(2) == 2.);
    CHECK(encoder.getElevation(2) == 2.);
    CHECK(encoder.getMute(2) == false);
    
    encoder.process(inputs.data(), ouputs.data());
}


TEST_CASE("MultiEncoder 3D", "[Encoder] [3D]")
{
    MultiEncoder<Hoa3d, float> encoder(7, 3);
    const float epsilon = FLT_EPSILON;
    std::vector<float> inputs(5);
    std::vector<float> ouputs(64);
    
    encoder.setRadius(0, 0.);
    encoder.setAzimuth(0, 0.);
    encoder.setElevation(0, 0.);
    encoder.setMute(0, false);
    
    encoder.setRadius(1, 1.);
    encoder.setAzimuth(1, 1.);
    encoder.setElevation(1, 1.);
    encoder.setMute(1, false);
    
    encoder.setRadius(2, 2.);
    encoder.setAzimuth(2, 2.);
    encoder.setElevation(2, 2.);
    encoder.setMute(2, false);
    
    CHECK(encoder.getRadius(0) == 0.);
    CHECK(encoder.getAzimuth(0) == 0.);
    CHECK(encoder.getElevation(0) == 0.);
    CHECK(encoder.getMute(0) == false);
    
    CHECK(encoder.getRadius(1) == 1.);
    CHECK(encoder.getAzimuth(1) == 1.);
    CHECK(encoder.getElevation(1) == 1.);
    CHECK(encoder.getMute(1) == false);
    
    CHECK(encoder.getRadius(2) == 2.);
    CHECK(encoder.getAzimuth(2) == 2.);
    CHECK(encoder.getElevation(2) == 2.);
    CHECK(encoder.getMute(2) == false);
    
    encoder.process(inputs.data(), ouputs.data());
    
}

