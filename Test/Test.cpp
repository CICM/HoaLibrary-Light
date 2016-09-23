/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "../Sources/Hoa.hpp"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cassert>

static void test_binaural()
{
    const unsigned i_blck_size  = 128;
    const unsigned i_order      = 1;
    const unsigned i_input_nb   = 4;
    const unsigned i_output_nb  = 2;
    const unsigned i_nb_samples = 64000;

    srand(time(NULL));
    hoa::Decoder<hoa::Hoa3d, float>::Binaural *decoder = new hoa::Decoder<hoa::Hoa3d, float>::Binaural(i_order);
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
}



int main(int argc, char** argv)
{
    std::cout << "binaural...";
    test_binaural();
    std::cout << "ok\n";
    return 0;
}
