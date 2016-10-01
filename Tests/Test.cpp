/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <cfloat>
#include <cmath>
#include <iostream>

#include <Hoa.hpp>
using namespace hoa;

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

int main(int argc, char* const argv[])
{
    std::cout << "HOA Library Testing...\n";
    return Catch::Session().run(argc, argv);
}