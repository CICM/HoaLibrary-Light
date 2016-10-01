/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <limits>
#include <cmath>

#include <Hoa.hpp>
using namespace hoa;
typedef float hoa_float_t;

#define CATCH_CONFIG_PREFIX_ALL
#include "catch.hpp"
#include <iostream>

CATCH_TEST_CASE("Optim 2D", "[Optim] [2D]")
{
    Optim<Hoa2d, hoa_float_t> optim(7);
    const hoa_float_t epsilon = std::numeric_limits<hoa_float_t>::epsilon();
    std::vector<hoa_float_t> inputs(15);
    std::vector<hoa_float_t> ouputs(15);
    CATCH_SECTION("Basic")
    {
        optim.setMode(Optim<Hoa2d, hoa_float_t>::Basic);
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
        optim.setMode(Optim<Hoa2d, hoa_float_t>::MaxRe);
        std::fill(inputs.begin(), inputs.end(), 1.);
        optim.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(std::abs(ouputs[0]  - inputs[0]) < epsilon);
        CATCH_CHECK(std::abs(ouputs[1]  - inputs[1]  * hoa_float_t(0.9807852506637573)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2]  - inputs[2]  * hoa_float_t(0.9807852506637573)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3]  - inputs[3]  * hoa_float_t(0.9238795042037964)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[4]  - inputs[4]  * hoa_float_t(0.9238795042037964)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5]  - inputs[5]  * hoa_float_t(0.8314695954322815)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6]  - inputs[6]  * hoa_float_t(0.8314695954322815)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7]  - inputs[7]  * hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8]  - inputs[8]  * hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[9]  - inputs[9]  * hoa_float_t(0.5555701851844788)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - inputs[10] * hoa_float_t(0.5555701851844788)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - inputs[11] * hoa_float_t(0.3826834261417389)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - inputs[12] * hoa_float_t(0.3826834261417389)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - inputs[13] * hoa_float_t(0.1950902342796326)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - inputs[14] * hoa_float_t(0.1950902342796326)) < epsilon);
    }

    CATCH_SECTION("InPhase")
    {
        optim.setMode(Optim<Hoa2d, hoa_float_t>::InPhase);
        std::fill(inputs.begin(), inputs.end(), 1.);
        optim.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(std::abs(ouputs[0]  - inputs[0]) < epsilon);
        CATCH_CHECK(std::abs(ouputs[1]  - inputs[1]  * hoa_float_t(0.8750000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2]  - inputs[2]  * hoa_float_t(0.8750000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3]  - inputs[3]  * hoa_float_t(0.5833333134651184)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[4]  - inputs[4]  * hoa_float_t(0.5833333134651184)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5]  - inputs[5]  * hoa_float_t(0.2916666567325592)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6]  - inputs[6]  * hoa_float_t(0.2916666567325592)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7]  - inputs[7]  * hoa_float_t(0.1060606092214584)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8]  - inputs[8]  * hoa_float_t(0.1060606092214584)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[9]  - inputs[9]  * hoa_float_t(0.0265151523053646)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - inputs[10] * hoa_float_t(0.0265151523053646)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - inputs[11] * hoa_float_t(0.0040792538784444)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - inputs[12] * hoa_float_t(0.0040792538784444)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - inputs[13] * hoa_float_t(0.0002913752978201)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - inputs[14] * hoa_float_t(0.0002913752978201)) < epsilon);
    }
}

CATCH_TEST_CASE("Optim 3D", "[Optim] [3D]")
{
    Optim<Hoa3d, hoa_float_t> optim(7);
    const hoa_float_t epsilon = std::numeric_limits<hoa_float_t>::epsilon();
    std::vector<hoa_float_t> inputs(64);
    std::vector<hoa_float_t> ouputs(64);
    
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
        optim.setMode(Optim<Hoa3d, hoa_float_t>::MaxRe);
        std::fill(inputs.begin(), inputs.end(), 1.);
        optim.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(ouputs[0] == inputs[0]);
        CATCH_CHECK(std::abs(ouputs[1]  - inputs[1]  * hoa_float_t(0.9807852506637573)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2]  - inputs[2]  * hoa_float_t(0.9807852506637573)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3]  - inputs[3]  * hoa_float_t(0.9807852506637573)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[4]  - inputs[4]  * hoa_float_t(0.9238795042037964)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5]  - inputs[5]  * hoa_float_t(0.9238795042037964)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6]  - inputs[6]  * hoa_float_t(0.9238795042037964)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7]  - inputs[7]  * hoa_float_t(0.9238795042037964)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8]  - inputs[8]  * hoa_float_t(0.9238795042037964)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[9]  - inputs[9]  * hoa_float_t(0.8314695954322815)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - inputs[10] * hoa_float_t(0.8314695954322815)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - inputs[11] * hoa_float_t(0.8314695954322815)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - inputs[12] * hoa_float_t(0.8314695954322815)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - inputs[13] * hoa_float_t(0.8314695954322815)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - inputs[14] * hoa_float_t(0.8314695954322815)) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[49] - inputs[49] * hoa_float_t(0.1950902342796326)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[63] - inputs[63] * hoa_float_t(0.1950902342796326)) < epsilon);
    }
    
    CATCH_SECTION("InPhase")
    {
        optim.setMode(Optim<Hoa3d, hoa_float_t>::InPhase);
        std::fill(inputs.begin(), inputs.end(), 1.);
        optim.process(inputs.data(), ouputs.data());
        
        CATCH_CHECK(ouputs[0] == inputs[0]);
        CATCH_CHECK(std::abs(ouputs[1]  - inputs[1]  * hoa_float_t(0.8750000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[2]  - inputs[2]  * hoa_float_t(0.8750000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[3]  - inputs[3]  * hoa_float_t(0.8750000000000000)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[4]  - inputs[4]  * hoa_float_t(0.5833333134651184)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[5]  - inputs[5]  * hoa_float_t(0.5833333134651184)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[6]  - inputs[6]  * hoa_float_t(0.5833333134651184)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[7]  - inputs[7]  * hoa_float_t(0.5833333134651184)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[8]  - inputs[8]  * hoa_float_t(0.5833333134651184)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[9]  - inputs[9]  * hoa_float_t(0.2916666567325592)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[10] - inputs[10] * hoa_float_t(0.2916666567325592)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[11] - inputs[11] * hoa_float_t(0.2916666567325592)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[12] - inputs[12] * hoa_float_t(0.2916666567325592)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[13] - inputs[13] * hoa_float_t(0.2916666567325592)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[14] - inputs[14] * hoa_float_t(0.2916666567325592)) < epsilon);
        
        CATCH_CHECK(std::abs(ouputs[49] - inputs[49] * hoa_float_t(0.0002913752978201956)) < epsilon);
        CATCH_CHECK(std::abs(ouputs[63] - inputs[63] * hoa_float_t(0.0002913752978201956)) < epsilon);
    }
}
