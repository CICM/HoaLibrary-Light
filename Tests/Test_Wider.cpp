/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <Hoa.hpp>
#include "catch.hpp"
#include <cfloat>

using namespace hoa;
typedef float hoa_float_t;

TEST_CASE("Wider 2D", "[Wider] [2D]")
{
    Wider<Hoa2d, hoa_float_t> wider(7);
    const float epsilon = FLT_EPSILON;
    std::vector<hoa_float_t> inputs(15);
    std::vector<hoa_float_t> ouputs(15);
    
    for(size_t i = 0; i < 15; ++i){
        inputs[i] = hoa_float_t(i%2 == 0);
    }
    
    SECTION("Widening factor 1")
    {
        wider.setWidening(1.);
        wider.process(inputs.data(), ouputs.data());
        
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
    
    SECTION("Widening factor 0")
    {
        wider.setWidening(0.);
        wider.process(inputs.data(), ouputs.data());
        
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
    
    SECTION("Widening factor 0.25")
    {
        wider.setWidening(0.25);
        wider.process(inputs.data(), ouputs.data());
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CHECK(sum_of_elems == 8);
    }
    
    SECTION("Widening factor 0.5")
    {
        wider.setWidening(0.5);
        wider.process(inputs.data(), ouputs.data());
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CHECK(sum_of_elems == 8);
    }
    
    SECTION("Widening factor 0.75")
    {
        wider.setWidening(0.75);
        wider.process(inputs.data(), ouputs.data());
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        CHECK(sum_of_elems == 8);
    }
}


TEST_CASE("Wider 3D", "[Wider] [3D]")
{
    Wider<Hoa3d, hoa_float_t> wider(7);
    const float epsilon = FLT_EPSILON;
    std::vector<hoa_float_t> inputs(64);
    std::vector<hoa_float_t> ouputs(64);
    
    for(size_t i = 0; i < 64; ++i){
        inputs[i] = 1.;
    }
    
    SECTION("Widening factor 1")
    {
        wider.setWidening(1.);
        wider.process(inputs.data(), ouputs.data());
        
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
    
    SECTION("Widening factor 0")
    {
        wider.setWidening(0.);
        wider.process(inputs.data(), ouputs.data());
        
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
        CHECK(ouputs[15] == 0.);
        
        CHECK(ouputs[49] == 0.);
        CHECK(ouputs[63] == 0.);
    }
    
    
    SECTION("Widening factor 0.5")
    {
        wider.setWidening(0.5);
        wider.process(inputs.data(), ouputs.data());
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        int todo_after_encoder;
        //CHECK(sum_of_elems == 8);
    }
    
    SECTION("Widening factor 0.75")
    {
        wider.setWidening(0.75);
        wider.process(inputs.data(), ouputs.data());
        hoa_float_t sum_of_elems = 0.;
        for(std::vector<hoa_float_t>::iterator it = ouputs.begin(); it != ouputs.end(); ++it)
            sum_of_elems += *it;
        
        int todo_after_encoder;
        //CHECK(sum_of_elems == 8);
    }
}
