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

CATCH_TEST_CASE("Wider 2D", "[Wider] [2D]")
{
    Wider<Hoa2d, float> wider(7);
    const float epsilon = FLT_EPSILON;
    std::vector<float> inputs(15);
    std::vector<float> ouputs(15);
    
    for(size_t i = 0; i < 15; ++i){
        inputs[i] = float(i%2 == 0);
    }
    
    CATCH_SECTION("Widening factor 1")
    {
        wider.setWidening(1.);
        wider.process(inputs.data(), ouputs.data());
        
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
    
    CATCH_SECTION("Widening factor 0")
    {
        wider.setWidening(0.);
        wider.process(inputs.data(), ouputs.data());
        
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
    
    CATCH_SECTION("Widening factor 0.25")
    {
        wider.setWidening(0.25);
        wider.process(inputs.data(), ouputs.data());
        float sum_of_elems = 0.;
        for(std::vector<float>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CATCH_CHECK(sum_of_elems == 8);
    }
    
    CATCH_SECTION("Widening factor 0.5")
    {
        wider.setWidening(0.5);
        wider.process(inputs.data(), ouputs.data());
        float sum_of_elems = 0.;
        for(std::vector<float>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CATCH_CHECK(sum_of_elems == 8);
    }
    
    CATCH_SECTION("Widening factor 0.75")
    {
        wider.setWidening(0.75);
        wider.process(inputs.data(), ouputs.data());
        float sum_of_elems = 0.;
        for(std::vector<float>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CATCH_CHECK(sum_of_elems == 8);
    }
}


CATCH_TEST_CASE("Wider 3D", "[Wider] [3D]")
{
    Wider<Hoa3d, float> wider(7);
    const float epsilon = FLT_EPSILON;
    std::vector<float> inputs(64);
    std::vector<float> ouputs(64);
    
    for(size_t i = 0; i < 64; ++i){
        inputs[i] = 1.;
    }
    
    CATCH_SECTION("Widening factor 1")
    {
        wider.setWidening(1.);
        wider.process(inputs.data(), ouputs.data());
        
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
    
    CATCH_SECTION("Widening factor 0")
    {
        wider.setWidening(0.);
        wider.process(inputs.data(), ouputs.data());
        
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
        CATCH_CHECK(ouputs[15] == 0.);
        
        CATCH_CHECK(ouputs[49] == 0.);
        CATCH_CHECK(ouputs[63] == 0.);
    }
    
    
    CATCH_SECTION("Widening factor 0.5")
    {
        wider.setWidening(0.5);
        wider.process(inputs.data(), ouputs.data());
        float sum_of_elems = 0.;
        for(std::vector<float>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        int todo_after_encoder;
        //CATCH_CHECK(sum_of_elems == 8);
    }
    
    CATCH_SECTION("Widening factor 0.75")
    {
        wider.setWidening(0.75);
        wider.process(inputs.data(), ouputs.data());
        float sum_of_elems = 0.;
        for(std::vector<float>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        int todo_after_encoder;
        //CATCH_CHECK(sum_of_elems == 8);
    }
}
