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

CATCH_TEST_CASE("ProcessorPlanewaves 2D", "[Planewaves] [Processor] [2D]")
{
    typedef ProcessorPlanewaves<Hoa2d, hoa_float_t> Processor;
    Processor proc(8);
    hoa_float_t const epsilon = std::numeric_limits<hoa_float_t>::epsilon() * hoa_float_t(10.);
    hoa_float_t const zero(0.);
    hoa_float_t const one(1.);
    
    CATCH_CHECK(proc.getPlanewaveIndex(0) == 1);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationX()) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationY()) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationZ()) < epsilon);
    
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(0, false) - hoa_float_t(0. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(0, true) - hoa_float_t(0. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(1) == 2);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(1, false) - hoa_float_t(1. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(1, true) - hoa_float_t(1. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(2) == 3);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(2, false) - hoa_float_t(2. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(2, true) - hoa_float_t(2. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(3) == 4);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(3, false) - hoa_float_t(3. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(3, true) - hoa_float_t(3. * HOA_2PI / 8.)) < epsilon);
    
    CATCH_CHECK(proc.getPlanewaveIndex(4) == 5);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(4, false) - hoa_float_t(4. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(4, true) - hoa_float_t(4. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(5) == 6);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(5, false) - hoa_float_t(5. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(5, true) - hoa_float_t(5. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(6) == 7);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(6, false) - hoa_float_t(6. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(6, true) - hoa_float_t(6. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(7) == 8);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(7, false) - hoa_float_t(7. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(7, true) - hoa_float_t(7. * HOA_2PI / 8.)) < epsilon);
    
    proc.setPlanewavesRotation(zero, zero, hoa_float_t(HOA_PI));
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationX()) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationY()) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationZ() - hoa_float_t(HOA_PI)) < epsilon);
    
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(0, false) - hoa_float_t(0. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(0, true) - hoa_float_t(0. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(1) == 2);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(1, false) - hoa_float_t(1. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(1, true) - hoa_float_t(1. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(2) == 3);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(2, false) - hoa_float_t(2. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(2, true) - hoa_float_t(2. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(3) == 4);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(3, false) - hoa_float_t(3. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(3, true) - hoa_float_t(3. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    
    CATCH_CHECK(proc.getPlanewaveIndex(4) == 5);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(4, false) - hoa_float_t(4. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(4, true) - hoa_float_t(4. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(5) == 6);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(5, false) - hoa_float_t(5. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(5, true) - hoa_float_t(5. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(6) == 7);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(6, false) - hoa_float_t(6. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(6, true) - hoa_float_t(6. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(7) == 8);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(7, false) - hoa_float_t(7. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(7, true) - hoa_float_t(7. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    
    proc.setPlanewavesRotation(zero, zero, zero);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationX()) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationY()) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationZ()) < epsilon);
    
    proc.setPlanewaveAzimuth(0, hoa_float_t(0. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(1, hoa_float_t(1. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(2, hoa_float_t(2. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(3, hoa_float_t(3. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(4, hoa_float_t(4. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(5, hoa_float_t(5. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(6, hoa_float_t(6. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(7, hoa_float_t(7. * HOA_2PI / 8. + HOA_PI));
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(0, false) - hoa_float_t(0. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(0, true) - hoa_float_t(0. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(1) == 2);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(1, false) - hoa_float_t(1. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(1, true) - hoa_float_t(1. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(2) == 3);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(2, false) - hoa_float_t(2. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(2, true) - hoa_float_t(2. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(3) == 4);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(3, false) - hoa_float_t(3. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(3, true) - hoa_float_t(3. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    
    CATCH_CHECK(proc.getPlanewaveIndex(4) == 5);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(4, false) - hoa_float_t(4. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(4, true) - hoa_float_t(4. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(5) == 6);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(5, false) - hoa_float_t(5. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(5, true) - hoa_float_t(5. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(6) == 7);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(6, false) - hoa_float_t(6. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(6, true) - hoa_float_t(6. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(7) == 8);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(7, false) - hoa_float_t(7. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(7, true) - hoa_float_t(7. * HOA_2PI / 8. - HOA_PI)) < epsilon);
}

CATCH_TEST_CASE("ProcessorPlanewaves 3D", "[Planewaves] [Processor] [3D]")
{
    typedef ProcessorPlanewaves<Hoa3d, hoa_float_t> Processor;
    Processor proc(8);
    hoa_float_t const epsilon = std::numeric_limits<hoa_float_t>::epsilon() * hoa_float_t(10.);
    hoa_float_t const zero(0.);
    hoa_float_t const one(1.);
    
    CATCH_CHECK(proc.getPlanewaveIndex(0) == 1);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationX()) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationY()) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationZ()) < epsilon);
    
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(0, false) - hoa_float_t(0. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(0, true) - hoa_float_t(0. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(1) == 2);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(1, false) - hoa_float_t(1. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(1, true) - hoa_float_t(1. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(2) == 3);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(2, false) - hoa_float_t(2. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(2, true) - hoa_float_t(2. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(3) == 4);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(3, false) - hoa_float_t(3. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(3, true) - hoa_float_t(3. * HOA_2PI / 8.)) < epsilon);
    
    CATCH_CHECK(proc.getPlanewaveIndex(4) == 5);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(4, false) - hoa_float_t(4. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(4, true) - hoa_float_t(4. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(5) == 6);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(5, false) - hoa_float_t(5. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(5, true) - hoa_float_t(5. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(6) == 7);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(6, false) - hoa_float_t(6. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(6, true) - hoa_float_t(6. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(7) == 8);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(7, false) - hoa_float_t(7. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(7, true) - hoa_float_t(7. * HOA_2PI / 8.)) < epsilon);
    
    proc.setPlanewavesRotation(zero, zero, hoa_float_t(HOA_PI));
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationX()) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationY()) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationZ() - hoa_float_t(HOA_PI)) < epsilon);
    
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(0, false) - hoa_float_t(0. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(0, true) - hoa_float_t(0. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(1) == 2);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(1, false) - hoa_float_t(1. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(1, true) - hoa_float_t(1. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(2) == 3);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(2, false) - hoa_float_t(2. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(2, true) - hoa_float_t(2. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(3) == 4);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(3, false) - hoa_float_t(3. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(3, true) - hoa_float_t(3. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    
    CATCH_CHECK(proc.getPlanewaveIndex(4) == 5);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(4, false) - hoa_float_t(4. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(4, true) - hoa_float_t(4. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(5) == 6);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(5, false) - hoa_float_t(5. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(5, true) - hoa_float_t(5. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(6) == 7);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(6, false) - hoa_float_t(6. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(6, true) - hoa_float_t(6. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(7) == 8);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(7, false) - hoa_float_t(7. * HOA_2PI / 8.)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(7, true) - hoa_float_t(7. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    
    proc.setPlanewavesRotation(zero, zero, zero);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationX()) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationY()) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewavesRotationZ()) < epsilon);
    
    proc.setPlanewaveAzimuth(0, hoa_float_t(0. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(1, hoa_float_t(1. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(2, hoa_float_t(2. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(3, hoa_float_t(3. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(4, hoa_float_t(4. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(5, hoa_float_t(5. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(6, hoa_float_t(6. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(7, hoa_float_t(7. * HOA_2PI / 8. + HOA_PI));
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(0, false) - hoa_float_t(0. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(0, true) - hoa_float_t(0. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(1) == 2);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(1, false) - hoa_float_t(1. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(1, true) - hoa_float_t(1. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(2) == 3);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(2, false) - hoa_float_t(2. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(2, true) - hoa_float_t(2. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(3) == 4);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(3, false) - hoa_float_t(3. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(3, true) - hoa_float_t(3. * HOA_2PI / 8. + HOA_PI)) < epsilon);
    
    CATCH_CHECK(proc.getPlanewaveIndex(4) == 5);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(4, false) - hoa_float_t(4. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(4, true) - hoa_float_t(4. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(5) == 6);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(5, false) - hoa_float_t(5. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(5, true) - hoa_float_t(5. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(6) == 7);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(6, false) - hoa_float_t(6. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(6, true) - hoa_float_t(6. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(proc.getPlanewaveIndex(7) == 8);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(7, false) - hoa_float_t(7. * HOA_2PI / 8. - HOA_PI)) < epsilon);
    CATCH_CHECK(std::abs(proc.getPlanewaveAzimuth(7, true) - hoa_float_t(7. * HOA_2PI / 8. - HOA_PI)) < epsilon);
}

