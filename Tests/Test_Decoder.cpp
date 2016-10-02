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

#define CATCH_CONFIG_PREFIX_ALL
#include "catch.hpp"


CATCH_TEST_CASE("Decoder 2D", "[Decoder] [2D]")
{
    Encoder<Hoa2d, hoa_float_t> encoder(7);
    const hoa_float_t epsilon = std::numeric_limits<hoa_float_t>::epsilon() * hoa_float_t(10);
    hoa_float_t              input(1.);
    std::vector<hoa_float_t> harmonics(15);
    std::vector<hoa_float_t> outputs(15);
    
    CATCH_SECTION("Regular")
    {
        /*
        DecoderRegular<Hoa2d, hoa_float_t> decoder(7, 16);
        encoder.setAzimuth(hoa_float_t(0.));
        encoder.process(&input, harmonics.data());
        decoder.process(harmonics.data(), outputs.data());
        
        CATCH_CHECK(std::abs(outputs[0]  - hoa_float_t(0.93750000000000000)) < epsilon);
        CATCH_CHECK(std::abs(outputs[1]  - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[2]  - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[3]  - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[4]  - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[5]  - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[6]  - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[7]  - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[8]  - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[9]  - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[10] - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[11] - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[12] - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[13] - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[14] - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[15] - hoa_float_t(0.06249989569187164)) < epsilon);
        */
        //test harmonic 0
        // test
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
        /*
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
        CATCH_CHECK(std::abs(outputs[0]  - hoa_float_t(0.93750000000000000)) < epsilon);
        CATCH_CHECK(std::abs(outputs[1]  - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[2]  - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[3]  - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[4]  - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[5]  - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[6]  - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[7]  - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[8]  - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[9]  - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[10] - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[11] - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[12] - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[13] - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[14] - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[15] - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[16] - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[17] - hoa_float_t(0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[18] - hoa_float_t(-0.06249989569187164)) < epsilon);
        CATCH_CHECK(std::abs(outputs[19] - hoa_float_t(0.06249989569187164)) < epsilon);
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
        /*
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
        /*
         CATCH_CHECK(std::abs(outputs[0]  - hoa_float_t(0.93750000000000000)) < epsilon);
         CATCH_CHECK(std::abs(outputs[1]  - hoa_float_t(0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[2]  - hoa_float_t(-0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[3]  - hoa_float_t(0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[4]  - hoa_float_t(-0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[5]  - hoa_float_t(0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[6]  - hoa_float_t(-0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[7]  - hoa_float_t(0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[8]  - hoa_float_t(-0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[9]  - hoa_float_t(0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[10] - hoa_float_t(-0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[11] - hoa_float_t(0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[12] - hoa_float_t(-0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[13] - hoa_float_t(0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[14] - hoa_float_t(-0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[15] - hoa_float_t(0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[16] - hoa_float_t(-0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[17] - hoa_float_t(0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[18] - hoa_float_t(-0.06249989569187164)) < epsilon);
         CATCH_CHECK(std::abs(outputs[19] - hoa_float_t(0.06249989569187164)) < epsilon);
         */
    }
}

