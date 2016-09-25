/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <Hoa.hpp>
#include "catch.hpp"

using namespace hoa;
typedef float hoa_float_t;

#define HOA_ORDER 7

TEST_CASE("Optim 2D", "[Optim] [2D]")
{
    std::vector<hoa_float_t> inputs(Harmonic<Hoa2d, hoa_float_t>::getNumberOfHarmonics(HOA_ORDER));
    std::vector<hoa_float_t> ouputs(Harmonic<Hoa2d, hoa_float_t>::getNumberOfHarmonics(HOA_ORDER));
    SECTION("Basic")
    {
        OptimBasic<Hoa2d, hoa_float_t> basic(HOA_ORDER);
        std::fill(inputs.begin(), inputs.end(), 1.);
        basic.process(inputs.data(), ouputs.data());
        for(int i = 0; i < Harmonic<Hoa2d, hoa_float_t>::getNumberOfHarmonics(HOA_ORDER); i++)
        {
            
        }
    }

    SECTION("MaxRe")
    {
        
    }

    SECTION("InPhase")
    {
        
    }
}
