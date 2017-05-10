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
template<typename A, typename B> bool hoa_compare(const A val1, const B val2)
{
    return std::abs(static_cast<hoa_float_t>(val1) - static_cast<hoa_float_t>(val2)) < (std::numeric_limits<hoa_float_t>::epsilon() * static_cast<hoa_float_t>(10.));
}

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
        CATCH_CHECK(hoa_compare(p1.getAzimuth(), zero));
        CATCH_CHECK(hoa_compare(p1.getElevation(), zero));
        CATCH_CHECK(hoa_compare(p1.getAbscissa(), zero));
        CATCH_CHECK(hoa_compare(p1.getOrdinate(), one));
        CATCH_CHECK(hoa_compare(p1.getHeight(), zero));
        
        Planewave p2(size_t(2), hoa_float_t(HOA_PI2), zero);
        CATCH_CHECK(p2.getIndex() == size_t(2));
        CATCH_CHECK(hoa_compare(p2.getAzimuth(), hoa_float_t(HOA_PI2)));
        CATCH_CHECK(hoa_compare(p2.getElevation(), zero));
        CATCH_CHECK(hoa_compare(p2.getAbscissa(), -one));
        CATCH_CHECK(hoa_compare(p2.getOrdinate(), zero));
        CATCH_CHECK(hoa_compare(p2.getHeight(), zero));
        
        Planewave p3(size_t(3), hoa_float_t(HOA_PI), zero);
        CATCH_CHECK(p3.getIndex() == size_t(3));
        CATCH_CHECK(hoa_compare(p3.getAzimuth(), hoa_float_t(HOA_PI)));
        CATCH_CHECK(hoa_compare(p3.getElevation(), zero));
        CATCH_CHECK(hoa_compare(p3.getAbscissa(), zero));
        CATCH_CHECK(hoa_compare(p3.getOrdinate(), -one));
        CATCH_CHECK(hoa_compare(p3.getHeight(), zero));
        
        Planewave p4(size_t(4), hoa_float_t(3. * HOA_PI2), zero);
        CATCH_CHECK(p4.getIndex() == size_t(4));
        CATCH_CHECK(hoa_compare(p4.getAzimuth(), hoa_float_t(3. * HOA_PI2)));
        CATCH_CHECK(hoa_compare(p4.getElevation(), zero));
        CATCH_CHECK(hoa_compare(p4.getAbscissa(), one));
        CATCH_CHECK(hoa_compare(p4.getOrdinate(), zero));
        CATCH_CHECK(hoa_compare(p4.getHeight(), zero));
        
        p4.setAzimuth(zero);
        CATCH_CHECK(hoa_compare(p4.getAzimuth(), zero));
        CATCH_CHECK(hoa_compare(p4.getElevation(), zero));
        CATCH_CHECK(hoa_compare(p4.getAbscissa(), zero));
        CATCH_CHECK(hoa_compare(p4.getOrdinate(), one));
        CATCH_CHECK(hoa_compare(p4.getHeight(), zero));
    }
    
    CATCH_SECTION("cartesian Constructor")
    {
        Planewave p1(size_t(1), zero, one, zero);
        CATCH_CHECK(p1.getIndex() == size_t(1));
        CATCH_CHECK(hoa_compare(p1.getAzimuth(), zero));
        CATCH_CHECK(hoa_compare(p1.getElevation(), zero));
        CATCH_CHECK(hoa_compare(p1.getAbscissa(), zero));
        CATCH_CHECK(hoa_compare(p1.getOrdinate(), one));
        CATCH_CHECK(hoa_compare(p1.getHeight(), zero));
        
        Planewave p2(size_t(2), -one, zero, zero);
        CATCH_CHECK(p2.getIndex() == size_t(2));
        CATCH_CHECK(hoa_compare(p2.getAzimuth(), hoa_float_t(HOA_PI2)));
        CATCH_CHECK(hoa_compare(p2.getElevation(), zero));
        CATCH_CHECK(hoa_compare(p2.getAbscissa(), -one));
        CATCH_CHECK(hoa_compare(p2.getOrdinate(), zero));
        CATCH_CHECK(hoa_compare(p2.getHeight(), zero));
        
        Planewave p3(size_t(3), zero, -one, zero);
        CATCH_CHECK(p3.getIndex() == size_t(3));
        CATCH_CHECK(hoa_compare(p3.getAzimuth(), hoa_float_t(HOA_PI)));
        CATCH_CHECK(hoa_compare(p3.getElevation(), zero));
        CATCH_CHECK(hoa_compare(p3.getAbscissa(), zero));
        CATCH_CHECK(hoa_compare(p3.getOrdinate(), -one));
        CATCH_CHECK(hoa_compare(p3.getHeight(), zero));
        
        Planewave p4(size_t(4), one, zero, zero);
        CATCH_CHECK(p4.getIndex() == size_t(4));
        CATCH_CHECK(hoa_compare(p4.getAzimuth(), hoa_float_t(3. * HOA_PI2)));
        CATCH_CHECK(hoa_compare(p4.getElevation(), zero));
        CATCH_CHECK(hoa_compare(p4.getAbscissa(), one));
        CATCH_CHECK(hoa_compare(p4.getOrdinate(), zero));
        CATCH_CHECK(hoa_compare(p4.getHeight(), zero));
        
        p4.setAzimuth(zero);
        CATCH_CHECK(hoa_compare(p4.getAzimuth(), zero));
        CATCH_CHECK(hoa_compare(p4.getElevation(), zero));
        CATCH_CHECK(hoa_compare(p4.getAbscissa(), zero));
        CATCH_CHECK(hoa_compare(p4.getOrdinate(), one));
        CATCH_CHECK(hoa_compare(p4.getHeight(), zero));
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
        CATCH_CHECK(hoa_compare(p1.getAzimuth(), zero));
        CATCH_CHECK(hoa_compare(p1.getElevation(), zero));
        CATCH_CHECK(hoa_compare(p1.getAbscissa(), zero));
        CATCH_CHECK(hoa_compare(p1.getOrdinate(), one));
        CATCH_CHECK(hoa_compare(p1.getHeight(), zero));
        
        Planewave p2(size_t(2), hoa_float_t(HOA_PI2), hoa_float_t(HOA_PI4));
        CATCH_CHECK(p2.getIndex() == size_t(2));
        CATCH_CHECK(hoa_compare(p2.getAzimuth(), hoa_float_t(HOA_PI2)));
        CATCH_CHECK(hoa_compare(p2.getElevation(), hoa_float_t(HOA_PI4)));
        CATCH_CHECK(hoa_compare(p2.getAbscissa(), hoa_float_t(-0.7071067690849304)));
        CATCH_CHECK(hoa_compare(p2.getOrdinate(), zero));
        CATCH_CHECK(hoa_compare(p2.getHeight(), hoa_float_t(0.7071067690849304)));
        
        Planewave p3(size_t(3), hoa_float_t(HOA_PI), hoa_float_t(HOA_PI2));
        CATCH_CHECK(p3.getIndex() == size_t(3));
        CATCH_CHECK(hoa_compare(p3.getAzimuth(), hoa_float_t(HOA_PI)));
        CATCH_CHECK(hoa_compare(p3.getElevation(), hoa_float_t(HOA_PI2)));
        CATCH_CHECK(hoa_compare(p3.getAbscissa(), zero));
        CATCH_CHECK(hoa_compare(p3.getOrdinate(), zero));
        CATCH_CHECK(hoa_compare(p3.getHeight(), one));
        
        Planewave p4(size_t(4), hoa_float_t(3. * HOA_PI2), hoa_float_t(-HOA_PI4));
        CATCH_CHECK(p4.getIndex() == size_t(4));
        CATCH_CHECK(hoa_compare(p4.getAzimuth(), hoa_float_t(3. * HOA_PI2)));
        CATCH_CHECK(hoa_compare(p4.getElevation(), hoa_float_t(-HOA_PI4)));
        CATCH_CHECK(hoa_compare(p4.getAbscissa(), hoa_float_t(0.7071067690849304)));
        CATCH_CHECK(hoa_compare(p4.getOrdinate(), zero));
        CATCH_CHECK(hoa_compare(p4.getHeight(), hoa_float_t(-0.7071067690849304)));
        
        p4.setAzimuth(zero);
        p4.setElevation(hoa_float_t(HOA_PI4));
        CATCH_CHECK(hoa_compare(p4.getAzimuth(), zero));
        CATCH_CHECK(hoa_compare(p4.getElevation(), hoa_float_t(HOA_PI4)));
        CATCH_CHECK(hoa_compare(p4.getAbscissa(), zero));
        CATCH_CHECK(hoa_compare(p4.getOrdinate(), hoa_float_t(0.7071067690849304)));
        CATCH_CHECK(hoa_compare(p4.getHeight(), hoa_float_t(0.7071067690849304)));
    }
    
    CATCH_SECTION("cartesian Constructor")
    {
        Planewave p1(size_t(1), zero, one, zero);
        CATCH_CHECK(p1.getIndex() == size_t(1));
        CATCH_CHECK(hoa_compare(p1.getAzimuth(), zero));
        CATCH_CHECK(hoa_compare(p1.getElevation(), zero));
        CATCH_CHECK(hoa_compare(p1.getAbscissa(), zero));
        CATCH_CHECK(hoa_compare(p1.getOrdinate(), one));
        CATCH_CHECK(hoa_compare(p1.getHeight(), zero));
        
        Planewave p2(size_t(2), hoa_float_t(-0.7071067690849304), zero, hoa_float_t(0.7071067690849304));
        CATCH_CHECK(p2.getIndex() == size_t(2));
        CATCH_CHECK(hoa_compare(p2.getAzimuth(), hoa_float_t(HOA_PI2)));
        CATCH_CHECK(hoa_compare(p2.getElevation(), hoa_float_t(HOA_PI4)));
        CATCH_CHECK(hoa_compare(p2.getAbscissa(), hoa_float_t(-0.7071067690849304)));
        CATCH_CHECK(hoa_compare(p2.getOrdinate(), zero));
        CATCH_CHECK(hoa_compare(p2.getHeight(), hoa_float_t(0.7071067690849304)));
        
        Planewave p3(size_t(3), zero, zero, one);
        CATCH_CHECK(p3.getIndex() == size_t(3));
        CATCH_CHECK(hoa_compare(p3.getAzimuth(), zero));
        CATCH_CHECK(hoa_compare(p3.getElevation(), hoa_float_t(HOA_PI2)));
        CATCH_CHECK(hoa_compare(p3.getAbscissa(), zero));
        CATCH_CHECK(hoa_compare(p3.getOrdinate(), zero));
        CATCH_CHECK(hoa_compare(p3.getHeight(), one));
        
        Planewave p4(size_t(4), hoa_float_t(0.7071067690849304), zero, hoa_float_t(-0.7071067690849304));
        CATCH_CHECK(p4.getIndex() == size_t(4));
        CATCH_CHECK(hoa_compare(p4.getAzimuth(), hoa_float_t(3. * HOA_PI2)));
        CATCH_CHECK(hoa_compare(p4.getElevation(), hoa_float_t(-HOA_PI4)));
        CATCH_CHECK(hoa_compare(p4.getAbscissa(), hoa_float_t(0.7071067690849304)));
        CATCH_CHECK(hoa_compare(p4.getOrdinate(), zero));
        CATCH_CHECK(hoa_compare(p4.getHeight(), hoa_float_t(-0.7071067690849304)));
        
        p4.setAzimuth(zero);
        p4.setElevation(hoa_float_t(HOA_PI4));
        CATCH_CHECK(hoa_compare(p4.getAzimuth(), zero));
        CATCH_CHECK(hoa_compare(p4.getElevation(), hoa_float_t(HOA_PI4)));
        CATCH_CHECK(hoa_compare(p4.getAbscissa(), zero));
        CATCH_CHECK(hoa_compare(p4.getOrdinate(), hoa_float_t(0.7071067690849304)));
        CATCH_CHECK(hoa_compare(p4.getHeight(), hoa_float_t(0.7071067690849304)));
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
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationX(), 0));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationY(), 0));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationZ(), 0));
    
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(0, false), hoa_float_t(0. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(0, true), hoa_float_t(0. * HOA_2PI / 8.)));
    CATCH_CHECK(proc.getPlanewaveIndex(1) == 2);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(1, false), hoa_float_t(1. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(1, true), hoa_float_t(1. * HOA_2PI / 8.)));
    CATCH_CHECK(proc.getPlanewaveIndex(2) == 3);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(2, false), hoa_float_t(2. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(2, true), hoa_float_t(2. * HOA_2PI / 8.)));
    CATCH_CHECK(proc.getPlanewaveIndex(3) == 4);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(3, false), hoa_float_t(3. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(3, true), hoa_float_t(3. * HOA_2PI / 8.)));
    
    CATCH_CHECK(proc.getPlanewaveIndex(4) == 5);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(4, false), hoa_float_t(4. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(4, true), hoa_float_t(4. * HOA_2PI / 8.)));
    CATCH_CHECK(proc.getPlanewaveIndex(5) == 6);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(5, false), hoa_float_t(5. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(5, true), hoa_float_t(5. * HOA_2PI / 8.)));
    CATCH_CHECK(proc.getPlanewaveIndex(6) == 7);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(6, false), hoa_float_t(6. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(6, true), hoa_float_t(6. * HOA_2PI / 8.)));
    CATCH_CHECK(proc.getPlanewaveIndex(7) == 8);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(7, false), hoa_float_t(7. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(7, true), hoa_float_t(7. * HOA_2PI / 8.)));
    
    proc.setPlanewavesRotation(zero, zero, hoa_float_t(HOA_PI));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationX(), 0));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationY(), 0));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationZ(), hoa_float_t(HOA_PI)));
    
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(0, false), hoa_float_t(0. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(0, true), hoa_float_t(0. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(1) == 2);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(1, false), hoa_float_t(1. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(1, true), hoa_float_t(1. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(2) == 3);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(2, false), hoa_float_t(2. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(2, true), hoa_float_t(2. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(3) == 4);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(3, false), hoa_float_t(3. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(3, true), hoa_float_t(3. * HOA_2PI / 8. + HOA_PI)));
    
    CATCH_CHECK(proc.getPlanewaveIndex(4) == 5);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(4, false), hoa_float_t(4. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(4, true), hoa_float_t(4. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(5) == 6);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(5, false), hoa_float_t(5. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(5, true), hoa_float_t(5. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(6) == 7);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(6, false), hoa_float_t(6. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(6, true), hoa_float_t(6. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(7) == 8);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(7, false), hoa_float_t(7. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(7, true), hoa_float_t(7. * HOA_2PI / 8. - HOA_PI)));
    
    proc.setPlanewavesRotation(zero, zero, zero);
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationX(), 0));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationY(), 0));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationZ(), 0));
    
    proc.setPlanewaveAzimuth(0, hoa_float_t(0. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(1, hoa_float_t(1. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(2, hoa_float_t(2. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(3, hoa_float_t(3. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(4, hoa_float_t(4. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(5, hoa_float_t(5. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(6, hoa_float_t(6. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(7, hoa_float_t(7. * HOA_2PI / 8. + HOA_PI));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(0, false), hoa_float_t(0. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(0, true), hoa_float_t(0. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(1) == 2);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(1, false), hoa_float_t(1. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(1, true), hoa_float_t(1. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(2) == 3);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(2, false), hoa_float_t(2. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(2, true), hoa_float_t(2. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(3) == 4);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(3, false), hoa_float_t(3. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(3, true), hoa_float_t(3. * HOA_2PI / 8. + HOA_PI)));
    
    CATCH_CHECK(proc.getPlanewaveIndex(4) == 5);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(4, false), hoa_float_t(4. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(4, true), hoa_float_t(4. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(5) == 6);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(5, false), hoa_float_t(5. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(5, true), hoa_float_t(5. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(6) == 7);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(6, false), hoa_float_t(6. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(6, true), hoa_float_t(6. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(7) == 8);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(7, false), hoa_float_t(7. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(7, true), hoa_float_t(7. * HOA_2PI / 8. - HOA_PI)));
}

CATCH_TEST_CASE("ProcessorPlanewaves 3D", "[Planewaves] [Processor] [3D]")
{
    typedef ProcessorPlanewaves<Hoa3d, hoa_float_t> Processor;
    Processor proc(8);
    hoa_float_t const epsilon = std::numeric_limits<hoa_float_t>::epsilon() * hoa_float_t(10.);
    hoa_float_t const zero(0.);
    hoa_float_t const one(1.);
    /*
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationX()));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationY()));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationZ()));
    
    CATCH_CHECK(proc.getPlanewaveIndex(0) == 1);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(0, false), hoa_float_t(0. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(0, true), hoa_float_t(0. * HOA_2PI / 8.)));
    CATCH_CHECK(proc.getPlanewaveIndex(1) == 2);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(1, false), hoa_float_t(1. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(1, true), hoa_float_t(1. * HOA_2PI / 8.)));
    CATCH_CHECK(proc.getPlanewaveIndex(2) == 3);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(2, false), hoa_float_t(2. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(2, true), hoa_float_t(2. * HOA_2PI / 8.)));
    CATCH_CHECK(proc.getPlanewaveIndex(3) == 4);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(3, false), hoa_float_t(3. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(3, true), hoa_float_t(3. * HOA_2PI / 8.)));
    
    CATCH_CHECK(proc.getPlanewaveIndex(4) == 5);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(4, false), hoa_float_t(4. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(4, true), hoa_float_t(4. * HOA_2PI / 8.)));
    CATCH_CHECK(proc.getPlanewaveIndex(5) == 6);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(5, false), hoa_float_t(5. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(5, true), hoa_float_t(5. * HOA_2PI / 8.)));
    CATCH_CHECK(proc.getPlanewaveIndex(6) == 7);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(6, false), hoa_float_t(6. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(6, true), hoa_float_t(6. * HOA_2PI / 8.)));
    CATCH_CHECK(proc.getPlanewaveIndex(7) == 8);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(7, false), hoa_float_t(7. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(7, true), hoa_float_t(7. * HOA_2PI / 8.)));
    
    proc.setPlanewavesRotation(zero, zero, hoa_float_t(HOA_PI));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationX()));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationY()));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationZ(), hoa_float_t(HOA_PI)));
    
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(0, false), hoa_float_t(0. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(0, true), hoa_float_t(0. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(1) == 2);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(1, false), hoa_float_t(1. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(1, true), hoa_float_t(1. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(2) == 3);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(2, false), hoa_float_t(2. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(2, true), hoa_float_t(2. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(3) == 4);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(3, false), hoa_float_t(3. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(3, true), hoa_float_t(3. * HOA_2PI / 8. + HOA_PI)));
    
    CATCH_CHECK(proc.getPlanewaveIndex(4) == 5);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(4, false), hoa_float_t(4. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(4, true), hoa_float_t(4. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(5) == 6);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(5, false), hoa_float_t(5. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(5, true), hoa_float_t(5. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(6) == 7);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(6, false), hoa_float_t(6. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(6, true), hoa_float_t(6. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(7) == 8);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(7, false), hoa_float_t(7. * HOA_2PI / 8.)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(7, true), hoa_float_t(7. * HOA_2PI / 8. - HOA_PI)));
    
    proc.setPlanewavesRotation(zero, zero, zero);
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationX()));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationY()));
    CATCH_CHECK(hoa_compare(proc.getPlanewavesRotationZ()));
    
    proc.setPlanewaveAzimuth(0, hoa_float_t(0. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(1, hoa_float_t(1. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(2, hoa_float_t(2. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(3, hoa_float_t(3. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(4, hoa_float_t(4. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(5, hoa_float_t(5. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(6, hoa_float_t(6. * HOA_2PI / 8. + HOA_PI));
    proc.setPlanewaveAzimuth(7, hoa_float_t(7. * HOA_2PI / 8. + HOA_PI));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(0, false), hoa_float_t(0. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(0, true), hoa_float_t(0. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(1) == 2);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(1, false), hoa_float_t(1. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(1, true), hoa_float_t(1. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(2) == 3);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(2, false), hoa_float_t(2. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(2, true), hoa_float_t(2. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(3) == 4);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(3, false), hoa_float_t(3. * HOA_2PI / 8. + HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(3, true), hoa_float_t(3. * HOA_2PI / 8. + HOA_PI)));
    
    CATCH_CHECK(proc.getPlanewaveIndex(4) == 5);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(4, false), hoa_float_t(4. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(4, true), hoa_float_t(4. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(5) == 6);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(5, false), hoa_float_t(5. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(5, true), hoa_float_t(5. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(6) == 7);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(6, false), hoa_float_t(6. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(6, true), hoa_float_t(6. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(proc.getPlanewaveIndex(7) == 8);
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(7, false), hoa_float_t(7. * HOA_2PI / 8. - HOA_PI)));
    CATCH_CHECK(hoa_compare(proc.getPlanewaveAzimuth(7, true), hoa_float_t(7. * HOA_2PI / 8. - HOA_PI)));
     */
}

