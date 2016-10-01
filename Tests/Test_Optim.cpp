/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <cfloat>
#include <cmath>

#include <Hoa.hpp>
using namespace hoa;

#define CATCH_CONFIG_PREFIX_ALL
#include "catch.hpp"

CATCH_TEST_CASE("Optim 2D", "[Optim] [2D]")
{
    Optim<Hoa2d, float> optim(7);
    const float epsilon = FLT_EPSILON;
    std::vector<float> inputs(15);
    std::vector<float> ouputs(15);
    CATCH_SECTION("Basic")
    {
        optim.setMode(Optim<Hoa2d, float>::Basic);
        std::fill(inputs.begin(), inputs.end(), 1.);
        optim.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(ouputs[0] == inputs[0]);
        CATCH_CHECK(ouputs[1] == inputs[1]);
        CATCH_CHECK(ouputs[2] == inputs[2]);
        CATCH_CHECK(ouputs[3] == inputs[3]);
        CATCH_CHECK(ouputs[4] == inputs[4]);
        CATCH_CHECK(ouputs[5] == inputs[5]);
        CATCH_CHECK(ouputs[6] == inputs[6]);
        CATCH_CHECK(ouputs[7] == inputs[7]);
        CATCH_CHECK(ouputs[8] == inputs[8]);
        CATCH_CHECK(ouputs[9] == inputs[9]);
        CATCH_CHECK(ouputs[10] == inputs[10]);
        CATCH_CHECK(ouputs[11] == inputs[11]);
        CATCH_CHECK(ouputs[12] == inputs[12]);
        CATCH_CHECK(ouputs[13] == inputs[13]);
        CATCH_CHECK(ouputs[14] == inputs[14]);
    }

    CATCH_SECTION("MaxRe")
    {
        optim.setMode(Optim<Hoa2d, float>::MaxRe);
        std::fill(inputs.begin(), inputs.end(), 1.);
        optim.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(ouputs[0] == inputs[0]);
        CATCH_CHECK(std::abs(ouputs[1] - inputs[1] * 0.9807852506) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2] - inputs[2] * 0.9807852506) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3] - inputs[3] * 0.9238795042) < epsilon);
        CATCH_CHECK(std::abs(ouputs[4] - inputs[4] * 0.9238795042) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5] - inputs[5] * 0.8314695954) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6] - inputs[6] * 0.8314695954) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7] - inputs[7] * 0.7071067690) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8] - inputs[8] * 0.7071067690) < epsilon);
        CATCH_CHECK(std::abs(ouputs[9] - inputs[9] * 0.5555702448) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - inputs[10] * 0.5555702448) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - inputs[11] * 0.3826834261) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - inputs[12] * 0.3826834261) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - inputs[13] * 0.1950903236) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - inputs[14] * 0.1950903236) < epsilon);
    }

    CATCH_SECTION("InPhase")
    {
        optim.setMode(Optim<Hoa2d, float>::InPhase);
        std::fill(inputs.begin(), inputs.end(), 1.);
        optim.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(ouputs[0] == inputs[0]);
        CATCH_CHECK(std::abs(ouputs[1] - inputs[1] * 0.8750000000) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2] - inputs[2] * 0.8750000000) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3] - inputs[3] * 0.5833333134) < epsilon);
        CATCH_CHECK(std::abs(ouputs[4] - inputs[4] * 0.5833333134) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5] - inputs[5] * 0.2916666567) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6] - inputs[6] * 0.2916666567) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7] - inputs[7] * 0.1060606092) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8] - inputs[8] * 0.1060606092) < epsilon);
        CATCH_CHECK(std::abs(ouputs[9] - inputs[9] * 0.0265151523) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - inputs[10] * 0.0265151523) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - inputs[11] * 0.0040792538) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - inputs[12] * 0.0040792538) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - inputs[13] * 0.0002913752) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - inputs[14] * 0.0002913752) < epsilon);
    }
}

CATCH_TEST_CASE("Optim 3D", "[Optim] [3D]")
{
    Optim<Hoa3d, float> optim(7);
    const float epsilon = FLT_EPSILON;
    std::vector<float> inputs(64);
    std::vector<float> ouputs(64);
    CATCH_SECTION("Basic")
    {
        optim.setMode(Optim<Hoa3d, float>::Basic);
        std::fill(inputs.begin(), inputs.end(), 1.);
        optim.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(ouputs[0] == inputs[0]);
        CATCH_CHECK(ouputs[1] == inputs[1]);
        CATCH_CHECK(ouputs[2] == inputs[2]);
        CATCH_CHECK(ouputs[3] == inputs[3]);
        CATCH_CHECK(ouputs[4] == inputs[4]);
        CATCH_CHECK(ouputs[5] == inputs[5]);
        CATCH_CHECK(ouputs[6] == inputs[6]);
        CATCH_CHECK(ouputs[7] == inputs[7]);
        CATCH_CHECK(ouputs[8] == inputs[8]);
        CATCH_CHECK(ouputs[9] == inputs[9]);
        CATCH_CHECK(ouputs[10] == inputs[10]);
        CATCH_CHECK(ouputs[11] == inputs[11]);
        CATCH_CHECK(ouputs[12] == inputs[12]);
        CATCH_CHECK(ouputs[13] == inputs[13]);
        CATCH_CHECK(ouputs[14] == inputs[14]);
        CATCH_CHECK(ouputs[15] == inputs[15]);
        
        CATCH_CHECK(ouputs[49] == inputs[49]);
        CATCH_CHECK(ouputs[63] == inputs[63]);
    }
    
    CATCH_SECTION("MaxRe")
    {
        optim.setMode(Optim<Hoa3d, float>::MaxRe);
        std::fill(inputs.begin(), inputs.end(), 1.);
        optim.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(ouputs[0] == inputs[0]);
        CATCH_CHECK(std::abs(ouputs[1] - inputs[1] * 0.9807852506) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2] - inputs[2] * 0.9807852506) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3] - inputs[3] * 0.9807852506) < epsilon);
        CATCH_CHECK(std::abs(ouputs[4] - inputs[4] * 0.9238795042) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5] - inputs[5] * 0.9238795042) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6] - inputs[6] * 0.9238795042) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7] - inputs[7] * 0.9238795042) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8] - inputs[8] * 0.9238795042) < epsilon);
        CATCH_CHECK(std::abs(ouputs[9] - inputs[9] * 0.8314695954) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - inputs[10] * 0.8314695954) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - inputs[11] * 0.8314695954) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - inputs[12] * 0.8314695954) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - inputs[13] * 0.8314695954) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - inputs[14] * 0.8314695954) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[49] - inputs[49] * 0.1950903236) < epsilon);
        CATCH_CHECK(std::abs(ouputs[63] - inputs[63] * 0.1950903236) < epsilon);
    }
    
    CATCH_SECTION("InPhase")
    {
        optim.setMode(Optim<Hoa3d, float>::InPhase);
        std::fill(inputs.begin(), inputs.end(), 1.);
        optim.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(ouputs[0] == inputs[0]);
        CATCH_CHECK(std::abs(ouputs[1] - inputs[1] * 0.8750000000) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2] - inputs[2] * 0.8750000000) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3] - inputs[3] * 0.8750000000) < epsilon);
        CATCH_CHECK(std::abs(ouputs[4] - inputs[4] * 0.5833333134) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5] - inputs[5] * 0.5833333134) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6] - inputs[6] * 0.5833333134) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7] - inputs[7] * 0.5833333134) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8] - inputs[8] * 0.5833333134) < epsilon);
        CATCH_CHECK(std::abs(ouputs[9] - inputs[9] * 0.2916666567) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - inputs[10] * 0.2916666567) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - inputs[11] * 0.2916666567) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - inputs[12] * 0.2916666567) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - inputs[13] * 0.2916666567) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - inputs[14] * 0.2916666567) < epsilon);
        CATCH_CHECK(std::abs(ouputs[15] - inputs[15] * 0.2916666567) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[49] - inputs[49] * 0.0002913752) < epsilon);
        CATCH_CHECK(std::abs(ouputs[63] - inputs[63] * 0.0002913752) < epsilon);
    }
}
