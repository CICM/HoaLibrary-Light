/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_UTILS__
#define __DEF_HOA_UTILS__

#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string>
#include <assert.h>

namespace Hoa
{
	//! The int to string conversion
    /** The function converts a integer to a string.

	 @param     value   The value to convert.
	 @return    The function return value in a string format.
     */
    inline std::string int_to_string(int aValue)
    {
        char number[256];
        sprintf(number, "%i", aValue);
        return number;
    }

	class BinauralFilter
    {
    private:
        
        class Delay
        {
        private:
            
            float              m_buffer[128];
            int                m_ramp;
            int                m_delay;
            float              m_sample_rate;
        public:
            
            Delay()
            {
                for(int i = 0; i < 128; i++)
                    m_buffer[i] = 0;
                m_ramp      = 0;
                m_delay     = 0;
                m_sample_rate = 44100.;
            }
            
            void setSamplingRate(float sampleRate)
            {
                if(sampleRate != m_sample_rate)
                {
                    float delay = (float)m_delay / sampleRate;
                    m_sample_rate = sampleRate;
                    setDelay(delay);
                }
            }
            
            void setDelayMs(float delay)
            {
                setDelay(delay * m_sample_rate);
            }
            
            void setDelay(int delay)
            {
                m_delay = delay;
                if(m_delay > 128)
                    m_delay = 128;
                m_ramp      = 0;
                for(int i = 0; i < 128; i++)
                    m_buffer[i] = 0;
            }
            
            float process(float sample)
            {
                int delay = m_ramp - m_delay;
                m_buffer[m_ramp] = sample;
                if(++m_ramp >= 128)
                    m_ramp = 0;
                if(delay < 0)
                    delay += 128;
                return m_buffer[delay];
            }
            
            ~Delay(){}
        };
        
        class Biquad
        {
        protected:
            float  m_coeff_a0, m_coeff_a1, m_coeff_a2, m_coeff_b1, m_coeff_b2;
            float  m_delay_one, m_delay_two;
            float  m_weight;
            float  m_sample_rate;
            
        public:
            Biquad()
            {
                m_sample_rate = 44100.f;
                m_delay_one = 0.f;
                m_delay_two = 0.f;
                setCoefficients(0.f, 0.f, 0.f, 0.f, 0.f);
            };
            
            void setSamplingRate(float sampleRate)
            {
                if(m_sample_rate != sampleRate)
                {
                    m_sample_rate = sampleRate;
                    m_delay_one = 0.f;
                    m_delay_two = 0.f;
                    setCoefficients(m_coeff_a0, m_coeff_a1, m_coeff_a2, m_coeff_b1, m_coeff_b2);
                }
            }
            
            void setCoefficients(float a0, float a1, float a2, float b1, float b2)
            {
                m_coeff_a0 = a0 * m_sample_rate / 44100.f;
                m_coeff_a1 = a1 * m_sample_rate / 44100.f;
                m_coeff_a2 = a2 * m_sample_rate / 44100.f;
                m_coeff_b1 = b1 * m_sample_rate / 44100.f;
                m_coeff_b2 = b2 * m_sample_rate / 44100.f;
            }
            
            void setCoefficients(float* coeffs1, float* coeffs2, float factor, float* coeffs3, float factor2)
            {
                m_coeff_a0 = coeffs1[0] * factor + coeffs2[0] * (1.f - factor);
                m_coeff_a1 = coeffs1[1] * factor + coeffs2[1] * (1.f - factor);
                m_coeff_a2 = coeffs1[2] * factor + coeffs2[2] * (1.f - factor);
                m_coeff_b1 = coeffs1[3] * factor + coeffs2[3] * (1.f - factor);
                m_coeff_b2 = coeffs1[4] * factor + coeffs2[4] * (1.f - factor);
                m_coeff_a0 = coeffs3[0] * factor2 + m_coeff_a0 * (1.f - factor2);
                m_coeff_a1 = coeffs3[1] * factor2 + m_coeff_a1 * (1.f - factor2);
                m_coeff_a2 = coeffs3[2] * factor2 + m_coeff_a2 * (1.f - factor2);
                m_coeff_b1 = coeffs3[3] * factor2 + m_coeff_b1 * (1.f - factor2);
                m_coeff_b2 = coeffs3[4] * factor2 + m_coeff_b2 * (1.f - factor2);
                m_coeff_a0 *= m_sample_rate / 44100.f;
                m_coeff_a1 *= m_sample_rate / 44100.f;
                m_coeff_a2 *= m_sample_rate / 44100.f;
                m_coeff_b1 *= m_sample_rate / 44100.f;
                m_coeff_b2 *= m_sample_rate / 44100.f;
            }
            
            inline float process(float input)
            {
                float output = input * m_coeff_a0 + m_delay_one;
                m_delay_one = input * m_coeff_a1 + m_delay_two - m_coeff_b1 * output;
                m_delay_two = input * m_coeff_a2 - m_coeff_b2 * output;
                return output;
            }
            
            ~Biquad(){};
        };
        
        Delay              m_delay;
        Biquad             m_biquad[6];
        double             m_width;
        double             m_lenght_front;
        double             m_lenght_back;
        double             m_height;
        double             m_gain;
        
    public:
        
        BinauralFilter(double azimuth, double elevation)
        {
            float a0[] =
            {
                0.152034 ,0. ,-0.152034 ,-1.815326 ,0.914481,
                0.881248 ,-0.713042 ,0.854105 ,-0.713042 ,0.735353,
                1.090837 ,0.490579 ,0.72749 ,0.490579 ,0.818326,
                1.089543 ,0.491276 ,0.731369 ,0.491276 ,0.820913,
                1.002012 ,-1.27624 ,0.841009 ,-1.27624 ,0.843021,
                0.930155 ,1.706645 ,0.786791 ,1.693536 ,0.730055
            };
            
            float a90[] =
            {
                0.149582 ,0. ,-0.149582 ,-1.572131 ,0.682139,
                0.802786 ,-0.281769 ,0.792406 ,-0.281769 ,0.595192,
                1.200574 ,0.84932 ,0.398279 ,0.84932 ,0.598853,
                0.737594 ,0.925952 ,0.660416 ,0.925952 ,0.39801,
                0.908982 ,-1.236645 ,0.828076 ,-1.236645 ,0.737058,
                0.587956 ,1.464997 ,1. ,1.464997 ,0.587956
            };
            
            float a180[] =
            {
                0.135893 ,0. ,-0.135893 ,-1.622599 ,0.711227,
                0.920398 ,-0.295651 ,0.891452 ,-0.295651 ,0.81185,
                1.200574 ,0.84932 ,0.398279 ,0.84932 ,0.598853,
                0.737594 ,0.925952 ,0.660416 ,0.925952 ,0.39801,
                0.962813 ,-1.286172 ,0.843813 ,-1.286172 ,0.806626,
                0.667577 ,0.702615 ,0.307889 ,0.35535 ,0.322731
            };
            
            float a270[] =
            {
                0.140982 ,0. ,-0.140982 ,-1.683363 ,0.77531,
                1. ,-0.31358 ,0.921726 ,-0.31358 ,0.921726, //
                1.0196 ,0.791713 ,0.470808 ,0.791713 ,0.490408, //
                0.915346 ,0.977908 ,0.717821 ,0.977908 ,0.633167,
                0.962813 ,-1.286172 ,0.843813 ,-1.286172 ,0.806626, //
                0.246902 ,0.294791 ,0.081504 ,-0.343893 ,-0.032909
            };
            
            float e90[] =
            {
                /*
                0.863762, -1.121228, 0.711174, -1.121228, 0.574936,
                1.127982, -1.536144, 0.527789, -1.536144, 0.655771, //
                0.430404, 0.860807, 0.430404, 0.351578, 0.370037, //
                0.821834, 1.236388, 0.663464, 1.236388, 0.485298, //
                0.902467, -1.765776, 0.863522, -1.762164, 0.7696, //
                0.644363, 0.30387, 0.539764, 0.30387, 0.184127 //
                 */
                0.140982 ,0. ,-0.140982 ,-1.683363 ,0.77531,
               //0.337491, 0.674981, 0.337491, 0.147023, 0.202939,
                1.093091, -1.21119, 0.301819, -1.21119, 0.39491,
                0.992113, -0.28786, 0.930987, -0.28786, 0.923099,
                0.998075, 0.830722, 0.964389, 0.830722, 0.962464,
                1., -0.415965, 0.918568, -0.415965, 0.918568,
                1., 2., 1., 2., 1.
            };
            
            float e180[] =
            {
                /*
                0.878635, -1.14604, 0.851665, -1.14604, 0.730301,
                1.044199, -1.87162, 0.867404, -1.87162, 0.911603, //
                0.092146, 0.184291, 0.092146, -1.213448, 0.58203, //
                0.902024, 0.973359, 0.891432, 0.973359, 0.793455, //
                0.902467, -1.765776, 0.863522, -1.762164, 0.7696, //
                0.615178, -0.167428, 0.614407, -0.167428, 0.229585 //
                 */
                0.046686, 0.093371, 0.046686, -1.482523, 0.669265,
                0.903331, -0.920447, 0.868178, -0.920447, 0.771509,
                0.902094, -0.266711, 0.879716, -0.266711, 0.781811,
                0.951943, 0.801276, 0.940959, 0.801276, 0.892902,
                0.946868, -0.40153, 0.905122, -0.40153, 0.85199,
                0.683488, 1.120627, 0.48604, 0.906685, 0.38347
            };
            
            //int delay[] = {9, 0, 11, 27};
            
            m_width         = 0.5;
            m_lenght_front  = 0.3;
            m_lenght_back   = 0.4;
            m_height        = 0.4;
            azimuth = wrap_twopi(azimuth);
            elevation = clip_minmax(elevation, -HOA_PI2, HOA_PI2);
            
            float gain1 = (fabs(cos((azimuth - HOA_PI2) * 0.5)) + 0.1) / 1.1;
            float gain2 = cos(fabs(elevation) * 0.5);
            float gainfact = (HOA_PI2 - fabs(elevation)) / HOA_PI2;
            if(azimuth > HOA_PI2 && azimuth < 3. * HOA_PI2)
            {
                gain1 *= (m_lenght_front / m_lenght_back);
                m_delay.setDelayMs(distance(-m_width, 0, abscissa(m_width, azimuth, elevation), ordinate(m_lenght_back, azimuth, elevation), 0, height(m_height, azimuth, elevation)) / 340. * 0.2);
            }
            else
            {
                m_delay.setDelayMs(distance(-m_width, 0, abscissa(m_width, azimuth, elevation), ordinate(m_lenght_front, azimuth, elevation), 0, height(m_height, azimuth, elevation)) / 340. * 0.2);
            }
            m_gain = gain1 * gainfact + gain2 * (1. - gainfact);
            
            float f1 = 1, f2 = 0;
            float *tab1 = a0, *tab2 = a270, *tab3 = e90;
            // Between 7Pi/4 (-Pi/4) & Pi/4
            if(azimuth >= HOA_2PI - HOA_PI4)
            {
                f1 = (HOA_2PI - azimuth) / HOA_PI4;
                tab1 = a0; tab2 = a270;
            }
            else if(azimuth >= 0.f && azimuth <= HOA_PI4)
            {
                f1 = azimuth / HOA_PI4;
                tab1 = a0; tab2 = a90;
            }
            // Between Pi/4 & 3Pi/4
            else if(azimuth >= HOA_PI4 && azimuth <= HOA_PI2)
            {
                f1 = (HOA_PI2 - azimuth) / HOA_PI4;
                tab1 = a90; tab2 = a0;
            }
            else if(azimuth >= HOA_PI2 && azimuth <= HOA_PI2 + HOA_PI4)
            {
                f1 = (azimuth - HOA_PI2) / HOA_PI4;
                tab1 = a90; tab2 = a180;
            }
            // Between 3Pi/4 & 5Pi/4
            else if(azimuth >= HOA_PI2 + HOA_PI4 && azimuth <= HOA_PI)
            {
                f1 = (HOA_PI - azimuth) / HOA_PI4;
                tab1 = a180; tab2 = a90;
            }
            else if(azimuth >= HOA_PI && azimuth <= HOA_PI + HOA_PI4)
            {
                f1 = (azimuth - HOA_PI) / HOA_PI4;
                tab1 = a180; tab2 = a270;
            }
            // Between 5Pi/4 & 7Pi/4
            else if(azimuth >= HOA_PI + HOA_PI4 && azimuth <= HOA_PI + HOA_PI2)
            {
                f1 = (HOA_PI + HOA_PI2 - azimuth) / HOA_PI4;
                tab1 = a270; tab2 = a180;
            }
            else if(azimuth >= HOA_PI + HOA_PI2 && azimuth <= HOA_PI + HOA_PI2 + HOA_PI4)
            {
                f1 = (azimuth - (HOA_PI + HOA_PI2)) / HOA_PI4;
                tab1 = a270; tab2 = a0;
            }
            
            if(elevation >= 0)
            {
                f2 = elevation / HOA_PI2;
                tab3 = e90;
            }
            else
            {
                f2 = -elevation / HOA_PI2;
                tab3 = e180;
            }
            
            for(int i = 0; i < 6; i++)
                m_biquad[i].setCoefficients(tab1+i*5, tab2+i*5, f1, tab3+i*5, f2);
            
            //freq = (ordinate(1, azimuth) + 1) * 3000. + 100;
        }
        
        void setSampleRate(double samplerate)
        {
            for(int i = 0; i < 6; i++)
                m_biquad[i].setSamplingRate(samplerate);
            m_delay.setSamplingRate(samplerate);
        }
        
        inline float process(const float sample)
        {
            return m_biquad[0].process(
                                       m_biquad[1].process(
                                                           m_biquad[2].process(
                                                                               m_biquad[3].process(
                                                                                                   m_biquad[4].process(
                                                                                                                       m_biquad[5].process(
                                                                                                                       m_delay.process(sample * m_gain)))))));
        }
        
        ~BinauralFilter()
        {
            ;
        }
    };
}

#endif


