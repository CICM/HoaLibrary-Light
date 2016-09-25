/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <Hoa.hpp>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cassert>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace hoa;
typedef float hoa_float_t;

TEST_CASE("Harmonics 2D", "[Harmonics] [2D]")
{
    typedef Harmonic<Hoa2d, hoa_float_t> Harmonic;
    Harmonic Harmonic0(0);
    Harmonic Harmonic1(1);
    Harmonic Harmonic2(2);
    Harmonic Harmonic3(3);
    Harmonic Harmonic4(4);
    Harmonic Harmonic5(5);
    Harmonic Harmonic6(6);
    Harmonic Harmonic13(13);
    Harmonic Harmonic14(14);
    Harmonic Harmonic21(21);
    Harmonic Harmonic22(22);
    
    SECTION("Static Degree")
    {
        CHECK(Harmonic::getDegree(0) == 0);
        CHECK(Harmonic::getDegree(1) == 1);
        CHECK(Harmonic::getDegree(2) == 1);
        CHECK(Harmonic::getDegree(3) == 2);
        CHECK(Harmonic::getDegree(4) == 2);
        CHECK(Harmonic::getDegree(5) == 3);
        CHECK(Harmonic::getDegree(6) == 3);
        CHECK(Harmonic::getDegree(13) == 7);
        CHECK(Harmonic::getDegree(14) == 7);
        CHECK(Harmonic::getDegree(21) == 11);
        CHECK(Harmonic::getDegree(22) == 11);
    }
    
    SECTION("Static Order")
    {
        CHECK(Harmonic::getOrder(0) == 0);
        CHECK(Harmonic::getOrder(1) == -1);
        CHECK(Harmonic::getOrder(2) == 1);
        CHECK(Harmonic::getOrder(3) == -2);
        CHECK(Harmonic::getOrder(4) == 2);
        CHECK(Harmonic::getOrder(5) == -3);
        CHECK(Harmonic::getOrder(6) == 3);
        CHECK(Harmonic::getOrder(13) == -7);
        CHECK(Harmonic::getOrder(14) == 7);
        CHECK(Harmonic::getOrder(21) == -11);
        CHECK(Harmonic::getOrder(22) == 11);
    }

    SECTION("Static Index")
    {
        CHECK(Harmonic::getIndex(0, 0) == 0);
        CHECK(Harmonic::getIndex(1, -1) == 1);
        CHECK(Harmonic::getIndex(1, 1) == 2);
        CHECK(Harmonic::getIndex(2, -2) == 3);
        CHECK(Harmonic::getIndex(2, 2) == 4);
        CHECK(Harmonic::getIndex(3, -3) == 5);
        CHECK(Harmonic::getIndex(3, 3) == 6);
        CHECK(Harmonic::getIndex(7, -7) == 13);
        CHECK(Harmonic::getIndex(7, 7) == 14);
        CHECK(Harmonic::getIndex(11, -11) == 21);
        CHECK(Harmonic::getIndex(11, 11) == 22);
    }
    
    SECTION("Static Number of Harmonics")
    {
        CHECK(Harmonic::getNumberOfHarmonics(0) == 1);
        CHECK(Harmonic::getNumberOfHarmonics(1) == 3);
        CHECK(Harmonic::getNumberOfHarmonics(2) == 5);
        CHECK(Harmonic::getNumberOfHarmonics(3) == 7);
        CHECK(Harmonic::getNumberOfHarmonics(7) == 15);
        CHECK(Harmonic::getNumberOfHarmonics(11) == 23);
    }
    
    SECTION("Static Normalization")
    {
        CHECK(Harmonic::getSemiNormalization(0, 0) == 1.);
        CHECK(Harmonic::getSemiNormalization(1, -1) == 1.);
        CHECK(Harmonic::getSemiNormalization(1, 1) == 1.);
        CHECK(Harmonic::getSemiNormalization(2, -2) == 1.);
        CHECK(Harmonic::getSemiNormalization(2, 2) == 1.);
        CHECK(Harmonic::getSemiNormalization(3, -3) == 1.);
        CHECK(Harmonic::getSemiNormalization(3, 3) == 1.);
        CHECK(Harmonic::getSemiNormalization(7, -7) == 1.);
        CHECK(Harmonic::getSemiNormalization(7, 7) == 1.);
        CHECK(Harmonic::getSemiNormalization(11, -11) == 1.);
        CHECK(Harmonic::getSemiNormalization(11, 11) == 1.);
    }
    
    SECTION("Static Semi Normalization")
    {
        CHECK(Harmonic::getSemiNormalization(0, 0) == 1.);
        CHECK(Harmonic::getSemiNormalization(1, -1) == 1.);
        CHECK(Harmonic::getSemiNormalization(1, 1) == 1.);
        CHECK(Harmonic::getSemiNormalization(2, -2) == 1.);
        CHECK(Harmonic::getSemiNormalization(2, 2) == 1.);
        CHECK(Harmonic::getSemiNormalization(3, -3) == 1.);
        CHECK(Harmonic::getSemiNormalization(3, 3) == 1.);
        CHECK(Harmonic::getSemiNormalization(7, -7) == 1.);
        CHECK(Harmonic::getSemiNormalization(7, 7) == 1.);
        CHECK(Harmonic::getSemiNormalization(11, -11) == 1.);
        CHECK(Harmonic::getSemiNormalization(11, 11) == 1.);
    }
    
    
    
    
    
    
    
    SECTION("Local Degree")
    {
        CHECK(Harmonic0.getDegree() == 0);
        CHECK(Harmonic1.getDegree() == 1);
        CHECK(Harmonic2.getDegree() == 1);
        CHECK(Harmonic3.getDegree() == 2);
        CHECK(Harmonic4.getDegree() == 2);
        CHECK(Harmonic5.getDegree() == 3);
        CHECK(Harmonic6.getDegree() == 3);
        CHECK(Harmonic13.getDegree() == 7);
        CHECK(Harmonic14.getDegree() == 7);
        CHECK(Harmonic21.getDegree() == 11);
        CHECK(Harmonic22.getDegree() == 11);
    }
    
    SECTION("Static Order")
    {
        CHECK(Harmonic0.getOrder() == 0);
        CHECK(Harmonic1.getOrder() == -1);
        CHECK(Harmonic2.getOrder() == 1);
        CHECK(Harmonic3.getOrder() == -2);
        CHECK(Harmonic4.getOrder() == 2);
        CHECK(Harmonic5.getOrder() == -3);
        CHECK(Harmonic6.getOrder() == 3);
        CHECK(Harmonic13.getOrder() == -7);
        CHECK(Harmonic14.getOrder() == 7);
        CHECK(Harmonic21.getOrder() == -11);
        CHECK(Harmonic22.getOrder() == 11);
    }
    
    SECTION("Static Index")
    {
        CHECK(Harmonic0.getIndex() == 0);
        CHECK(Harmonic1.getIndex() == 1);
        CHECK(Harmonic2.getIndex() == 2);
        CHECK(Harmonic3.getIndex() == 3);
        CHECK(Harmonic4.getIndex() == 4);
        CHECK(Harmonic5.getIndex() == 5);
        CHECK(Harmonic6.getIndex() == 6);
        CHECK(Harmonic13.getIndex() == 13);
        CHECK(Harmonic14.getIndex() == 14);
        CHECK(Harmonic21.getIndex() == 21);
        CHECK(Harmonic22.getIndex() == 22);
    }
    
    SECTION("Static Normalization")
    {
        CHECK(Harmonic0.getNormalization() == 1.);
        CHECK(Harmonic1.getNormalization() == 1.);
        CHECK(Harmonic2.getNormalization() == 1.);
        CHECK(Harmonic3.getNormalization() == 1.);
        CHECK(Harmonic4.getNormalization() == 1.);
        CHECK(Harmonic5.getNormalization() == 1.);
        CHECK(Harmonic6.getNormalization() == 1.);
        CHECK(Harmonic13.getNormalization() == 1.);
        CHECK(Harmonic14.getNormalization() == 1.);
        CHECK(Harmonic21.getNormalization() == 1.);
        CHECK(Harmonic22.getNormalization() == 1.);
    }
    
    SECTION("Static Semi Normalization")
    {
        CHECK(Harmonic0.getNormalization() == 1.);
        CHECK(Harmonic1.getNormalization() == 1.);
        CHECK(Harmonic2.getNormalization() == 1.);
        CHECK(Harmonic3.getNormalization() == 1.);
        CHECK(Harmonic4.getNormalization() == 1.);
        CHECK(Harmonic5.getNormalization() == 1.);
        CHECK(Harmonic6.getNormalization() == 1.);
        CHECK(Harmonic13.getNormalization() == 1.);
        CHECK(Harmonic14.getNormalization() == 1.);
        CHECK(Harmonic21.getNormalization() == 1.);
        CHECK(Harmonic22.getNormalization() == 1.);
    }
}


