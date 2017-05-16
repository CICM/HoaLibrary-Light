/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <limits>
#include <cmath>
#include <iostream>

#include <Hoa.hpp>
using namespace hoa;
typedef float hoa_float_t;
template<typename A, typename B> bool hoa_compare(const A val1, const B val2)
{
    return std::abs(static_cast<hoa_float_t>(val1) - static_cast<hoa_float_t>(val2)) < (std::numeric_limits<hoa_float_t>::epsilon() * static_cast<hoa_float_t>(10.));
}

#define CATCH_CONFIG_PREFIX_ALL
#include "catch.hpp"


CATCH_TEST_CASE("Decoder 2D", "[Decoder] [2D]")
{
    Encoder<Hoa2d, hoa_float_t> encoder(7);
    hoa_float_t              input(1.);
    std::vector<hoa_float_t> harmonics(15);
    std::vector<hoa_float_t> outputs(16);
    
    CATCH_SECTION("Regular")
    {
        DecoderRegular<Hoa2d, hoa_float_t> decoder(7, 16);
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.process(&input, harmonics.data());
        decoder.process(harmonics.data(), outputs.data());
        
        CATCH_CHECK(hoa_compare(outputs[0], 0.5));
        CATCH_CHECK(hoa_compare(outputs[1], 0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[2], -0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[3], 0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[4], -0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[5], 0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[6], -0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[7], 0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[8], -0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[9], 0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[10], -0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[11], 0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[12], -0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[13], 0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[14], -0.03333324939012527));
        CATCH_CHECK(hoa_compare(outputs[15], 0.03333324939012527));
    }
}

CATCH_TEST_CASE("Decoder 3D", "[Decoder] [3D]")
{
    Encoder<Hoa3d, hoa_float_t> encoder(3);
    const hoa_float_t epsilon = std::numeric_limits<hoa_float_t>::epsilon() * hoa_float_t(10);
    hoa_float_t              input(1.);
    std::vector<hoa_float_t> harmonics(20);
    std::vector<hoa_float_t> outputs(20);
    
    CATCH_SECTION("Regular")
    {
        std::cout.precision(16);
        DecoderRegular<Hoa3d, hoa_float_t> decoder(1, 4);
        /*
        decoder.setPlanewaveAzimuth(0, 0.);
        decoder.setPlanewaveElevation(0, 0.);
        decoder.setPlanewaveAzimuth(1, 90.);
        decoder.setPlanewaveElevation(1, 0.);
        decoder.setPlanewaveAzimuth(3, 180.);
        decoder.setPlanewaveElevation(3, 0.);
        decoder.setPlanewaveAzimuth(4, 270.);
        decoder.setPlanewaveElevation(4, 0.);
        decoder.prepare();
         *
        
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setElevation(hoa_float_t(0.));
        encoder.process(&input, harmonics.data());
        for(size_t i = 0; i < 4; ++i)
        {
            std::cout << harmonics[i] << "\n";
        }
        harmonics[3] *= 3.;
         
         
        for(size_t i = 0; i < harmonics.size(); ++i) {
            harmonics[i] = hoa_float_t(0.);
        }
        harmonics[0] = 1.;
        
        decoder.process(harmonics.data(), outputs.data());
        
        
        for(size_t i = 0; i < outputs.size() && i < decoder.getNumberOfPlanewaves(); ++i)
        {
            std::cout << decoder.getPlanewaveName(i) << ": " << outputs[i] << "\n";
        }
        std::cout << "ratio 1 : " << (outputs[0] / outputs[1]) << "\n";
        std::cout << "ratio 2 : " << (outputs[1] / outputs[2]) << "\n";
        /
        CATCH_CHECK(hoa_compare(outputs[0], 0.93750000000000000);
        CATCH_CHECK(hoa_compare(outputs[1], 0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[2], -0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[3], 0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[4], -0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[5], 0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[6], -0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[7], 0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[8], -0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[9], 0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[10], -0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[11], 0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[12], -0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[13], 0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[14], -0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[15], 0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[16], -0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[17], 0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[18], -0.06249989569187164);
        CATCH_CHECK(hoa_compare(outputs[19], 0.06249989569187164);
         */
    }
    
    CATCH_SECTION("Regular2")
    {
        std::cout.precision(16);
        DecoderRegular<Hoa3d, hoa_float_t> decoder(3, 20);
        /*
         decoder.setPlanewaveAzimuth(0, 0.);
         decoder.setPlanewaveElevation(0, 0.);
         decoder.setPlanewaveAzimuth(1, 90.);
         decoder.setPlanewaveElevation(1, 0.);
         decoder.setPlanewaveAzimuth(3, 180.);
         decoder.setPlanewaveElevation(3, 0.);
         decoder.setPlanewaveAzimuth(4, 270.);
         decoder.setPlanewaveElevation(4, 0.);
         decoder.prepare();
         *
        
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.setElevation(hoa_float_t(0.));
        encoder.process(&input, harmonics.data());
        
         
         
         for(size_t i = 0; i < harmonics.size(); ++i) {
         harmonics[i] = hoa_float_t(0.);
         }
         harmonics[0] = 1.;
         *
        decoder.process(harmonics.data(), outputs.data());
        
        
        for(size_t i = 0; i < outputs.size() && i < decoder.getNumberOfPlanewaves(); ++i)
        {
            std::cout << decoder.getPlanewaveName(i) << ": " << outputs[i] << "\n";
        }
        std::cout << "ratio 1 : " << (outputs[0] / outputs[1]) << "\n";
        std::cout << "ratio 2 : " << (outputs[1] / outputs[2]) << "\n";
        
         
         
         CATCH_CHECK(hoa_compare(outputs[0], 0.93750000000000000);
         CATCH_CHECK(hoa_compare(outputs[1], 0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[2], -0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[3], 0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[4], -0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[5], 0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[6], -0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[7], 0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[8], -0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[9], 0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[10], -0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[11], 0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[12], -0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[13], 0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[14], -0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[15], 0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[16], -0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[17], 0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[18], -0.06249989569187164);
         CATCH_CHECK(hoa_compare(outputs[19], 0.06249989569187164);
         */
    }
}

