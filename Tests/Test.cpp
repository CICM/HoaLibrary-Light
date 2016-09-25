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
    /*
     CHECK(Harmonic::getOrder(0) == 0);
     CHECK(Harmonic::getIndex(0, 0) == 0);
     
     CHECK(Harmonic::getDegree(0) == 0);
     CHECK(Harmonic::getOrder(0) == 0);
     CHECK(Harmonic::getIndex(0, 0) == 0);
     
     for(size_t i = 1; i < HOA_MAX_ORDER; ++i)
     {
     const size_t index1 = Harmonic::getNumberOfHarmonics(i) - 1;
     const size_t index2 = Harmonic::getNumberOfHarmonics(i);
     CHECK(index1 == i * 2);
     CHECK(index2 == i * 2 + 1);
     
     CHECK(Harmonic::getDegree(index1) == i);
     CHECK(Harmonic::getDegree(index2) == i);
     
     CHECK(Harmonic::getOrder(index1) == -long(i));
     CHECK(Harmonic::getOrder(index2) == i);
     
     CHECK(Harmonic::getIndex(Harmonic::getDegree(index1), Harmonic::getOrder(index1)) == index1);
     CHECK(Harmonic::getIndex(Harmonic::getDegree(index2), Harmonic::getOrder(index2)) == index2);
     }
     */
}

static void test_binaural()
{
    /*
    const unsigned i_blck_size  = 128;
    const unsigned i_order      = 1;
    const unsigned i_input_nb   = 4;
    const unsigned i_output_nb  = 2;
    const unsigned i_nb_samples = 64000;

    srand(time(NULL));
    hoa::DecoderBinaural<hoa::Hoa3d, float> *decoder = new hoa::DecoderBinaural<hoa::Hoa3d, float>(i_order);
    float *p_src    = (float *)malloc(i_input_nb * i_nb_samples * sizeof(float));
    float *p_dest   = (float *)malloc(i_output_nb * i_nb_samples * sizeof(float));
    memset(p_src, 0, i_input_nb * i_nb_samples * sizeof(float));
    memset(p_dest, 0, i_output_nb * i_nb_samples * sizeof(float));
    for(unsigned i = 0; i < i_input_nb * i_nb_samples; i += 4)
    {
        p_src[i] = rand() / RAND_MAX + 1.f;
    }

    float **in_buf  = (float **)malloc(i_input_nb * sizeof(float *));
    for (unsigned i = 0; i < i_input_nb; ++i)
    {
        in_buf[i] = (float *)malloc(i_blck_size * sizeof(float));
    }

    float **out_buf = (float **)malloc(i_output_nb * sizeof(float *));
    for (unsigned i = 0; i < i_output_nb; ++i)
    {
        out_buf[i] = (float *)malloc(i_blck_size * sizeof(float));
    }

    decoder->computeRendering(i_blck_size);

    for(unsigned i = 0; i < i_nb_samples; i += i_blck_size)
    {
        for(unsigned j = 0; j < i_input_nb; ++j)
        {
            for(unsigned k = 0; k < i_blck_size; ++k)
            {
                in_buf[j][i] = p_src[(i + k) * i_input_nb + j];
            }
        }

        decoder->processBlock((const float**)in_buf, (float**)out_buf);

        for(unsigned j = 0; j < i_output_nb; ++j)
        {
            for(unsigned k = 0; k < i_blck_size; ++k)
            {
                p_dest[(i + k) * i_output_nb + j] = out_buf[j][k];
            }
        }
    }


    for(unsigned i = 0; i < i_output_nb * i_nb_samples; ++i)
    {
        assert(p_dest[i] != 0.f && "problem ici");
    }

    for(unsigned i = 0; i < i_output_nb; ++i)
    {
        free(out_buf[i]);
    }
    free(out_buf);

    for(unsigned i = 0; i < i_input_nb; ++i)
    {
        free(in_buf[i]);
    }
    free(in_buf);

    free(p_src);
    free(p_dest);

    delete decoder;
     */
}

