/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <cfloat>
#include <cmath>

#include <Hoa.hpp>
using namespace hoa;
typedef float hoa_float_t;

#define CATCH_CONFIG_PREFIX_ALL
#include "catch.hpp"
#include <iostream>

CATCH_TEST_CASE("Planewaves 2D", "[Planewaves] [2D]")
{
    typedef Planewave<Hoa2d, hoa_float_t> Planewave;
    hoa_float_t const epsilon = std::numeric_limits<hoa_float_t>::epsilon() * hoa_float_t(10.);
    hoa_float_t const zero(0.);
    hoa_float_t const one(1.);
    
    CATCH_SECTION("Polar Constructor")
    {
        Planewave p1(size_t(1), zero, zero);
        CATCH_CHECK(p1.getIndex() == size_t(1));
        CATCH_CHECK(std::abs(p1.getAzimuth() - zero) < epsilon);
        CATCH_CHECK(std::abs(p1.getElevation() - zero) < epsilon);
        CATCH_CHECK(std::abs(p1.getAbscissa() - zero) < epsilon);
        CATCH_CHECK(std::abs(p1.getOrdinate() - one) < epsilon);
        CATCH_CHECK(std::abs(p1.getHeight() - zero) < epsilon);
        
        Planewave p2(size_t(2), hoa_float_t(HOA_PI2), zero);
        CATCH_CHECK(p2.getIndex() == size_t(2));
        CATCH_CHECK(std::abs(p2.getAzimuth() - hoa_float_t(HOA_PI2)) < epsilon);
        CATCH_CHECK(std::abs(p2.getElevation() - zero) < epsilon);
        CATCH_CHECK(std::abs(p2.getAbscissa() - -one) < epsilon);
        CATCH_CHECK(std::abs(p2.getOrdinate() - zero) < epsilon);
        CATCH_CHECK(std::abs(p2.getHeight() - zero) < epsilon);
        
        Planewave p3(size_t(3), hoa_float_t(HOA_PI), zero);
        CATCH_CHECK(p3.getIndex() == size_t(3));
        CATCH_CHECK(std::abs(p3.getAzimuth() - hoa_float_t(HOA_PI)) < epsilon);
        CATCH_CHECK(std::abs(p3.getElevation() - zero) < epsilon);
        CATCH_CHECK(std::abs(p3.getAbscissa() - zero) < epsilon);
        CATCH_CHECK(std::abs(p3.getOrdinate() - -one) < epsilon);
        CATCH_CHECK(std::abs(p3.getHeight() - zero) < epsilon);
        
        Planewave p4(size_t(4), hoa_float_t(3. * HOA_PI2), zero);
        CATCH_CHECK(p4.getIndex() == size_t(4));
        CATCH_CHECK(std::abs(p4.getAzimuth() - hoa_float_t(3. * HOA_PI2)) < epsilon);
        CATCH_CHECK(std::abs(p4.getElevation() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getAbscissa() - one) < epsilon);
        CATCH_CHECK(std::abs(p4.getOrdinate() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getHeight() - zero) < epsilon);
        
        p4.setAzimuth(zero);
        CATCH_CHECK(std::abs(p4.getAzimuth() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getElevation() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getAbscissa() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getOrdinate() - one) < epsilon);
        CATCH_CHECK(std::abs(p4.getHeight() - zero) < epsilon);
    }
    
    CATCH_SECTION("cartesian Constructor")
    {
        Planewave p1(size_t(1), zero, one, zero);
        CATCH_CHECK(p1.getIndex() == size_t(1));
        CATCH_CHECK(std::abs(p1.getAzimuth() - zero) < epsilon);
        CATCH_CHECK(std::abs(p1.getElevation() - zero) < epsilon);
        CATCH_CHECK(std::abs(p1.getAbscissa() - zero) < epsilon);
        CATCH_CHECK(std::abs(p1.getOrdinate() - one) < epsilon);
        CATCH_CHECK(std::abs(p1.getHeight() - zero) < epsilon);
        
        Planewave p2(size_t(2), -one, zero, zero);
        CATCH_CHECK(p2.getIndex() == size_t(2));
        CATCH_CHECK(std::abs(p2.getAzimuth() - hoa_float_t(HOA_PI2)) < epsilon);
        CATCH_CHECK(std::abs(p2.getElevation() - zero) < epsilon);
        CATCH_CHECK(std::abs(p2.getAbscissa() - -one) < epsilon);
        CATCH_CHECK(std::abs(p2.getOrdinate() - zero) < epsilon);
        CATCH_CHECK(std::abs(p2.getHeight() - zero) < epsilon);
        
        Planewave p3(size_t(3), zero, -one, zero);
        CATCH_CHECK(p3.getIndex() == size_t(3));
        CATCH_CHECK(std::abs(p3.getAzimuth() - hoa_float_t(HOA_PI)) < epsilon);
        CATCH_CHECK(std::abs(p3.getElevation() - zero) < epsilon);
        CATCH_CHECK(std::abs(p3.getAbscissa() - zero) < epsilon);
        CATCH_CHECK(std::abs(p3.getOrdinate() - -one) < epsilon);
        CATCH_CHECK(std::abs(p3.getHeight() - zero) < epsilon);
        
        Planewave p4(size_t(4), one, zero, zero);
        CATCH_CHECK(p4.getIndex() == size_t(4));
        CATCH_CHECK(std::abs(p4.getAzimuth() - hoa_float_t(3. * HOA_PI2)) < epsilon);
        CATCH_CHECK(std::abs(p4.getElevation() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getAbscissa() - one) < epsilon);
        CATCH_CHECK(std::abs(p4.getOrdinate() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getHeight() - zero) < epsilon);
        
        p4.setAzimuth(zero);
        CATCH_CHECK(std::abs(p4.getAzimuth() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getElevation() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getAbscissa() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getOrdinate() - one) < epsilon);
        CATCH_CHECK(std::abs(p4.getHeight() - zero) < epsilon);
    }
}

CATCH_TEST_CASE("Planewaves 3D", "[Planewaves] [3D]")
{
    typedef Planewave<Hoa3d, hoa_float_t> Planewave;
    hoa_float_t const epsilon = std::numeric_limits<hoa_float_t>::epsilon() * hoa_float_t(10.);
    hoa_float_t const zero(0.);
    hoa_float_t const one(1.);
    
    CATCH_SECTION("Polar Constructor")
    {
        Planewave p1(size_t(1), zero, zero);
        CATCH_CHECK(p1.getIndex() == size_t(1));
        CATCH_CHECK(std::abs(p1.getAzimuth() - zero) < epsilon);
        CATCH_CHECK(std::abs(p1.getElevation() - zero) < epsilon);
        CATCH_CHECK(std::abs(p1.getAbscissa() - zero) < epsilon);
        CATCH_CHECK(std::abs(p1.getOrdinate() - one) < epsilon);
        CATCH_CHECK(std::abs(p1.getHeight() - zero) < epsilon);
        
        Planewave p2(size_t(2), hoa_float_t(HOA_PI2), hoa_float_t(HOA_PI4));
        CATCH_CHECK(p2.getIndex() == size_t(2));
        CATCH_CHECK(std::abs(p2.getAzimuth() - hoa_float_t(HOA_PI2)) < epsilon);
        CATCH_CHECK(std::abs(p2.getElevation() - hoa_float_t(HOA_PI4)) < epsilon);
        CATCH_CHECK(std::abs(p2.getAbscissa() - hoa_float_t(-0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(p2.getOrdinate() - zero) < epsilon);
        CATCH_CHECK(std::abs(p2.getHeight() - hoa_float_t(0.7071067690849304)) < epsilon);
        
        Planewave p3(size_t(3), hoa_float_t(HOA_PI), hoa_float_t(HOA_PI2));
        CATCH_CHECK(p3.getIndex() == size_t(3));
        CATCH_CHECK(std::abs(p3.getAzimuth() - hoa_float_t(HOA_PI)) < epsilon);
        CATCH_CHECK(std::abs(p3.getElevation() - hoa_float_t(HOA_PI2)) < epsilon);
        CATCH_CHECK(std::abs(p3.getAbscissa() - zero) < epsilon);
        CATCH_CHECK(std::abs(p3.getOrdinate() - zero) < epsilon);
        CATCH_CHECK(std::abs(p3.getHeight() - one) < epsilon);
        
        Planewave p4(size_t(4), hoa_float_t(3. * HOA_PI2), hoa_float_t(-HOA_PI4));
        CATCH_CHECK(p4.getIndex() == size_t(4));
        CATCH_CHECK(std::abs(p4.getAzimuth() - hoa_float_t(3. * HOA_PI2)) < epsilon);
        CATCH_CHECK(std::abs(p4.getElevation() - hoa_float_t(-HOA_PI4)) < epsilon);
        CATCH_CHECK(std::abs(p4.getAbscissa() - hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(p4.getOrdinate() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getHeight() - hoa_float_t(-0.7071067690849304)) < epsilon);
        
        p4.setAzimuth(zero);
        p4.setElevation(hoa_float_t(HOA_PI4));
        CATCH_CHECK(std::abs(p4.getAzimuth() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getElevation() - hoa_float_t(HOA_PI4)) < epsilon);
        CATCH_CHECK(std::abs(p4.getAbscissa() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getOrdinate() - hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(p4.getHeight() - hoa_float_t(0.7071067690849304)) < epsilon);
    }
    
    CATCH_SECTION("cartesian Constructor")
    {
        Planewave p1(size_t(1), zero, one, zero);
        CATCH_CHECK(p1.getIndex() == size_t(1));
        CATCH_CHECK(std::abs(p1.getAzimuth() - zero) < epsilon);
        CATCH_CHECK(std::abs(p1.getElevation() - zero) < epsilon);
        CATCH_CHECK(std::abs(p1.getAbscissa() - zero) < epsilon);
        CATCH_CHECK(std::abs(p1.getOrdinate() - one) < epsilon);
        CATCH_CHECK(std::abs(p1.getHeight() - zero) < epsilon);
        
        Planewave p2(size_t(2), hoa_float_t(-0.7071067690849304), zero, hoa_float_t(0.7071067690849304));
        CATCH_CHECK(p2.getIndex() == size_t(2));
        CATCH_CHECK(std::abs(p2.getAzimuth() - hoa_float_t(HOA_PI2)) < epsilon);
        CATCH_CHECK(std::abs(p2.getElevation() - hoa_float_t(HOA_PI4)) < epsilon);
        CATCH_CHECK(std::abs(p2.getAbscissa() - hoa_float_t(-0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(p2.getOrdinate() - zero) < epsilon);
        CATCH_CHECK(std::abs(p2.getHeight() - hoa_float_t(0.7071067690849304)) < epsilon);
        
        Planewave p3(size_t(3), zero, zero, one);
        CATCH_CHECK(p3.getIndex() == size_t(3));
        CATCH_CHECK(std::abs(p3.getAzimuth() - zero) < epsilon);
        CATCH_CHECK(std::abs(p3.getElevation() - hoa_float_t(HOA_PI2)) < epsilon);
        CATCH_CHECK(std::abs(p3.getAbscissa() - zero) < epsilon);
        CATCH_CHECK(std::abs(p3.getOrdinate() - zero) < epsilon);
        CATCH_CHECK(std::abs(p3.getHeight() - one) < epsilon);
        
        Planewave p4(size_t(4), hoa_float_t(0.7071067690849304), zero, hoa_float_t(-0.7071067690849304));
        CATCH_CHECK(p4.getIndex() == size_t(4));
        CATCH_CHECK(std::abs(p4.getAzimuth() - hoa_float_t(3. * HOA_PI2)) < epsilon);
        CATCH_CHECK(std::abs(p4.getElevation() - hoa_float_t(-HOA_PI4)) < epsilon);
        CATCH_CHECK(std::abs(p4.getAbscissa() - hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(p4.getOrdinate() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getHeight() - hoa_float_t(-0.7071067690849304)) < epsilon);
        
        p4.setAzimuth(zero);
        p4.setElevation(hoa_float_t(HOA_PI4));
        CATCH_CHECK(std::abs(p4.getAzimuth() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getElevation() - hoa_float_t(HOA_PI4)) < epsilon);
        CATCH_CHECK(std::abs(p4.getAbscissa() - zero) < epsilon);
        CATCH_CHECK(std::abs(p4.getOrdinate() - hoa_float_t(0.7071067690849304)) < epsilon);
        CATCH_CHECK(std::abs(p4.getHeight() - hoa_float_t(0.7071067690849304)) < epsilon);
    }
}

