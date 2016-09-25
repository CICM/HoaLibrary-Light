/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <Hoa.hpp>
#include "catch.hpp"

#include <iostream>

using namespace hoa;
typedef float hoa_float_t;

TEST_CASE("Optim 2D", "[Optim] [2D]")
{
    const float epsilon = FLT_EPSILON;
    std::vector<hoa_float_t> inputs(15);
    std::vector<hoa_float_t> ouputs(15);
    SECTION("Basic")
    {
        OptimBasic<Hoa2d, hoa_float_t> basic(7);
        std::fill(inputs.begin(), inputs.end(), 1.);
        basic.process(inputs.data(), ouputs.data());
        
        CHECK(ouputs[0] == inputs[0]);
        CHECK(ouputs[1] == inputs[1]);
        CHECK(ouputs[2] == inputs[2]);
        CHECK(ouputs[3] == inputs[3]);
        CHECK(ouputs[4] == inputs[4]);
        CHECK(ouputs[5] == inputs[5]);
        CHECK(ouputs[6] == inputs[6]);
        CHECK(ouputs[7] == inputs[7]);
        CHECK(ouputs[8] == inputs[8]);
        CHECK(ouputs[9] == inputs[9]);
        CHECK(ouputs[10] == inputs[10]);
        CHECK(ouputs[11] == inputs[11]);
        CHECK(ouputs[12] == inputs[12]);
        CHECK(ouputs[13] == inputs[13]);
        CHECK(ouputs[14] == inputs[14]);
    }

    SECTION("MaxRe")
    {
        OptimMaxRe<Hoa2d, hoa_float_t> basic(7);
        std::fill(inputs.begin(), inputs.end(), 1.);
        basic.process(inputs.data(), ouputs.data());
        
        CHECK(ouputs[0] == inputs[0]);
        CHECK(fabs(ouputs[1] - inputs[1] * 0.9807852506) < epsilon);
        CHECK(fabs(ouputs[2] - inputs[2] * 0.9807852506) < epsilon);
        CHECK(fabs(ouputs[3] - inputs[3] * 0.9238795042) < epsilon);
        CHECK(fabs(ouputs[4] - inputs[4] * 0.9238795042) < epsilon);
        CHECK(fabs(ouputs[5] - inputs[5] * 0.8314695954) < epsilon);
        CHECK(fabs(ouputs[6] - inputs[6] * 0.8314695954) < epsilon);
        CHECK(fabs(ouputs[7] - inputs[7] * 0.7071067690) < epsilon);
        CHECK(fabs(ouputs[8] - inputs[8] * 0.7071067690) < epsilon);
        CHECK(fabs(ouputs[9] - inputs[9] * 0.5555702448) < epsilon);
        CHECK(fabs(ouputs[10] - inputs[10] * 0.5555702448) < epsilon);
        CHECK(fabs(ouputs[11] - inputs[11] * 0.3826834261) < epsilon);
        CHECK(fabs(ouputs[12] - inputs[12] * 0.3826834261) < epsilon);
        CHECK(fabs(ouputs[13] - inputs[13] * 0.1950903236) < epsilon);
        CHECK(fabs(ouputs[14] - inputs[14] * 0.1950903236) < epsilon);
    }

    SECTION("InPhase")
    {
        OptimInPhase<Hoa2d, hoa_float_t> basic(7);
        std::fill(inputs.begin(), inputs.end(), 1.);
        basic.process(inputs.data(), ouputs.data());
        
        CHECK(ouputs[0] == inputs[0]);
        CHECK(fabs(ouputs[1] - inputs[1] * 0.8750000000) < epsilon);
        CHECK(fabs(ouputs[2] - inputs[2] * 0.8750000000) < epsilon);
        CHECK(fabs(ouputs[3] - inputs[3] * 0.5833333134) < epsilon);
        CHECK(fabs(ouputs[4] - inputs[4] * 0.5833333134) < epsilon);
        CHECK(fabs(ouputs[5] - inputs[5] * 0.2916666567) < epsilon);
        CHECK(fabs(ouputs[6] - inputs[6] * 0.2916666567) < epsilon);
        CHECK(fabs(ouputs[7] - inputs[7] * 0.1060606092) < epsilon);
        CHECK(fabs(ouputs[8] - inputs[8] * 0.1060606092) < epsilon);
        CHECK(fabs(ouputs[9] - inputs[9] * 0.0265151523) < epsilon);
        CHECK(fabs(ouputs[10] - inputs[10] * 0.0265151523) < epsilon);
        CHECK(fabs(ouputs[11] - inputs[11] * 0.0040792538) < epsilon);
        CHECK(fabs(ouputs[12] - inputs[12] * 0.0040792538) < epsilon);
        CHECK(fabs(ouputs[13] - inputs[13] * 0.0002913752) < epsilon);
        CHECK(fabs(ouputs[14] - inputs[14] * 0.0002913752) < epsilon);
    }
}

TEST_CASE("Optim 3D", "[Optim] [3D]")
{
    const float epsilon = FLT_EPSILON;
    std::vector<hoa_float_t> inputs(64);
    std::vector<hoa_float_t> ouputs(64);
    SECTION("Basic")
    {
        OptimBasic<Hoa3d, hoa_float_t> basic(7);
        std::fill(inputs.begin(), inputs.end(), 1.);
        basic.process(inputs.data(), ouputs.data());
        
        CHECK(ouputs[0] == inputs[0]);
        CHECK(ouputs[1] == inputs[1]);
        CHECK(ouputs[2] == inputs[2]);
        CHECK(ouputs[3] == inputs[3]);
        CHECK(ouputs[4] == inputs[4]);
        CHECK(ouputs[5] == inputs[5]);
        CHECK(ouputs[6] == inputs[6]);
        CHECK(ouputs[7] == inputs[7]);
        CHECK(ouputs[8] == inputs[8]);
        CHECK(ouputs[9] == inputs[9]);
        CHECK(ouputs[10] == inputs[10]);
        CHECK(ouputs[11] == inputs[11]);
        CHECK(ouputs[12] == inputs[12]);
        CHECK(ouputs[13] == inputs[13]);
        CHECK(ouputs[14] == inputs[14]);
        CHECK(ouputs[15] == inputs[15]);
        
        CHECK(ouputs[49] == inputs[49]);
        CHECK(ouputs[63] == inputs[63]);
    }
    
    SECTION("MaxRe")
    {
        OptimMaxRe<Hoa3d, hoa_float_t> basic(7);
        std::fill(inputs.begin(), inputs.end(), 1.);
        basic.process(inputs.data(), ouputs.data());
        
        CHECK(ouputs[0] == inputs[0]);
        CHECK(fabs(ouputs[1] - inputs[1] * 0.9807852506) < epsilon);
        CHECK(fabs(ouputs[2] - inputs[2] * 0.9807852506) < epsilon);
        CHECK(fabs(ouputs[3] - inputs[3] * 0.9807852506) < epsilon);
        CHECK(fabs(ouputs[4] - inputs[4] * 0.9238795042) < epsilon);
        CHECK(fabs(ouputs[5] - inputs[5] * 0.9238795042) < epsilon);
        CHECK(fabs(ouputs[6] - inputs[6] * 0.9238795042) < epsilon);
        CHECK(fabs(ouputs[7] - inputs[7] * 0.9238795042) < epsilon);
        CHECK(fabs(ouputs[8] - inputs[8] * 0.9238795042) < epsilon);
        CHECK(fabs(ouputs[9] - inputs[9] * 0.8314695954) < epsilon);
        CHECK(fabs(ouputs[10] - inputs[10] * 0.8314695954) < epsilon);
        CHECK(fabs(ouputs[11] - inputs[11] * 0.8314695954) < epsilon);
        CHECK(fabs(ouputs[12] - inputs[12] * 0.8314695954) < epsilon);
        CHECK(fabs(ouputs[13] - inputs[13] * 0.8314695954) < epsilon);
        CHECK(fabs(ouputs[14] - inputs[14] * 0.8314695954) < epsilon);
        
        CHECK(fabs(ouputs[49] - inputs[49] * 0.1950903236) < epsilon);
        CHECK(fabs(ouputs[63] - inputs[63] * 0.1950903236) < epsilon);
    }
    
    SECTION("InPhase")
    {
        OptimInPhase<Hoa3d, hoa_float_t> basic(7);
        std::fill(inputs.begin(), inputs.end(), 1.);
        basic.process(inputs.data(), ouputs.data());
        
        CHECK(ouputs[0] == inputs[0]);
        CHECK(fabs(ouputs[1] - inputs[1] * 0.8750000000) < epsilon);
        CHECK(fabs(ouputs[2] - inputs[2] * 0.8750000000) < epsilon);
        CHECK(fabs(ouputs[3] - inputs[3] * 0.8750000000) < epsilon);
        CHECK(fabs(ouputs[4] - inputs[4] * 0.5833333134) < epsilon);
        CHECK(fabs(ouputs[5] - inputs[5] * 0.5833333134) < epsilon);
        CHECK(fabs(ouputs[6] - inputs[6] * 0.5833333134) < epsilon);
        CHECK(fabs(ouputs[7] - inputs[7] * 0.5833333134) < epsilon);
        CHECK(fabs(ouputs[8] - inputs[8] * 0.5833333134) < epsilon);
        CHECK(fabs(ouputs[9] - inputs[9] * 0.2916666567) < epsilon);
        CHECK(fabs(ouputs[10] - inputs[10] * 0.2916666567) < epsilon);
        CHECK(fabs(ouputs[11] - inputs[11] * 0.2916666567) < epsilon);
        CHECK(fabs(ouputs[12] - inputs[12] * 0.2916666567) < epsilon);
        CHECK(fabs(ouputs[13] - inputs[13] * 0.2916666567) < epsilon);
        CHECK(fabs(ouputs[14] - inputs[14] * 0.2916666567) < epsilon);
        CHECK(fabs(ouputs[15] - inputs[15] * 0.2916666567) < epsilon);
        
        CHECK(fabs(ouputs[49] - inputs[49] * 0.0002913752) < epsilon);
        CHECK(fabs(ouputs[63] - inputs[63] * 0.0002913752) < epsilon);
    }
}
