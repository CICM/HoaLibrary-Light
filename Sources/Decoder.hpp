/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_DECODER_LIGHT
#define DEF_HOA_DECODER_LIGHT

#include "Encoder.hpp"
#include "Hrtf.hpp"

namespace hoa
{
    //! The decoder class decodes a sound field in the harmonics domain through the planewave domain.
    /** The decoder should be used to decode a set the harmonics domain to a set of planewave for loudspeakers. There are three types of decoder. Regular for a perfect circle or sphere of loudspeakers. Irregular when the loudspeakers are not equally spaced over the circle or the sphere. Binaural for headphone restitution.
     */
    template <Dimension D, typename T> class Decoder : public Processor< Harmonic<D, T> >, public Processor< Planewave<D, T> >
    {
    public:
        //! The decoder constructor.
        /**	The decoder constructor allocates and initialize the base classes.
         @param     order                   The order
         @param     numberOfPlanewaves      The number of channels.
         */
        Decoder(const ulong order, const ulong numberOfPlanewaves) noexcept = 0;
        
        //! The destructor.
        /** The destructor free the memory.
         */
        virtual ~Decoder();
        
        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and performs the decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimym size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        virtual void process(const T* inputs, T* outputs) noexcept;
        
        //! This method computes the decoding matrice.
        /**	You should use this method after changing the position of the loudspeakers and/or calling the process method.
         @param vectorsize The vector size for binaural decoding.
         */
        virtual void computeRendering(const ulong vectorsize = 64);
        
        //! The regular decoder.
        /** The regular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if more or equal to the number of harmonics plus one and when the loudspeakers are equally spaced.
         */
        class Regular;
        
        //! The irregular decoder.
        /** The irregular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if less than the number of harmonics plus one or when the loudspeakers are not equally spaced.
         */
        class Irregular;
        
        //! The binaural decoder.
        /** The binaural decoder should be used to decode an ambisonic sound field for headphones.
         */
        class Binaural;
    };
    
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    
    template <typename T> class Decoder<Hoa2d, T> : public Processor< Harmonic<Hoa2d, T> >, public Processor< Planewave<Hoa2d, T> >
    {
    public:
        
        //! The decoder constructor.
        /**	The decoder constructor allocates and initialize the base classes.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        Decoder(const ulong order, const ulong numberOfPlanewaves) noexcept :
        Processor< Harmonic<Hoa2d, T> >(order),
        Processor< Planewave<Hoa2d, T> >(numberOfPlanewaves)
        {
            ;
        }
        
        //! The destructor.
        /** The destructor free the memory.
         */
        virtual ~Decoder()
        {
            ;
        }
        
        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and performs the decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimym size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        inline virtual void process(const T* inputs, T* outputs) noexcept = 0;
        
        //! This method computes the decoding matrice.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        virtual void computeRendering(const ulong vectorsize = 64) = 0;
        
        //! The ambisonic regular decoder.
        /** The regular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if more or equal to the number of harmonics plus one and when the loudspeakers are equally spaced.
         */
        class Regular;
        
        //! The ambisonic irregular decoder.
        /** The irregular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if less than the number of harmonics plus one or when the loudspeakers are not equally spaced.
         */
       class Irregular;
        
        //! The ambisonic binaural decoder.
        /** The binaural decoder should be used to decode an ambisonic sound field for headphones.
         */
        class Binaural;
    };
    
    template <typename T> class Decoder<Hoa2d, T>::Regular : public Decoder<Hoa2d, T>
    {
    private:
        T*  m_matrix;
    public:
        
        //! The regular constructor.
        /**	The regular constructor allocates and initialize the decoding matrix depending on a order of decomposition and a number of channels. The order must be at least 1 and the number of channels must be at least the number of harmonics.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        Regular(const ulong order, const ulong numberOfPlanewaves) noexcept : Decoder<Hoa2d, T>(order, numberOfPlanewaves)
        {
            m_matrix = new T[Decoder<Hoa2d, T>::getNumberOfPlanewaves() * Decoder<Hoa2d, T>::getNumberOfHarmonics()];
            computeRendering();
        }
        
        //! The destructor.
        /** The destructor free the memory.
         */
        ~Regular()
        {
            delete [] m_matrix;
        }
        
        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and performs the regular decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimym size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to channels.
         */
        inline void process(const T* inputs, T* outputs) noexcept override
        {
            Signal<T>::matrix_vector_mul(Decoder<Hoa2d, T>::getNumberOfHarmonics(), Decoder<Hoa2d, T>::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }
        
        //! This method computes the decoding matrice.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        void computeRendering(const ulong vectorsize = 64) override
        {
            typename Encoder<Hoa2d, T>::Basic encoder(Decoder<Hoa2d, T>::getDecompositionOrder());
            const T factor = 1. / (T)(Decoder<Hoa2d, T>::getDecompositionOrder() + 1.);
            for(ulong i = 0; i < Decoder<Hoa2d, T>::getNumberOfPlanewaves(); i++)
            {
                encoder.setAzimuth(Decoder<Hoa2d, T>::getPlanewaveAzimuth(i));
                encoder.process(&factor, m_matrix + i * Decoder<Hoa2d, T>::getNumberOfHarmonics());
                m_matrix[i * Encoder<Hoa2d, T>::Basic::getNumberOfHarmonics()] = factor * 0.5;
            }
        }
    };
    
    template <typename T> class Decoder<Hoa2d, T>::Irregular : public Decoder<Hoa2d, T>
    {
    private:
        T*  m_matrix;
    public:
        
        //! The irregular constructor.
        /**	The irregular constructor allocates and initialize the decoding matrix depending on a order of decomposition and a number of channels. The order must be at least 1.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        Irregular(const ulong order, const ulong numberOfPlanewaves) noexcept : Decoder<Hoa2d, T>(order, numberOfPlanewaves)
        {
            m_matrix = new T[Decoder<Hoa2d, T>::getNumberOfPlanewaves() * Decoder<Hoa2d, T>::getNumberOfHarmonics()];
            computeRendering();
        }
        
        //! The destructor.
        /** The destructor free the memory.
         */
        ~Irregular()
        {
            delete [] m_matrix;
        }
        
        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and performs the irregular decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimym size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        inline void process(const T* inputs, T* outputs) noexcept override
        {
            Signal<T>::matrix_vector_mul(Decoder<Hoa2d, T>::getNumberOfHarmonics(), Decoder<Hoa2d, T>::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }
        
        //! This method computes the decoding matrice.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        void computeRendering(const ulong vectorsize = 64)  override
        {
            typename Encoder<Hoa2d, T>::Basic encoder(Decoder<Hoa2d, T>::getDecompositionOrder());
            Signal<T>::vector_clear(Decoder<Hoa2d, T>::getNumberOfPlanewaves() * Decoder<Hoa2d, T>::getNumberOfHarmonics(), m_matrix);
            T vector_harmonics[Decoder<Hoa2d, T>::getNumberOfHarmonics()];
            
            if(Decoder<Hoa2d, T>::getNumberOfPlanewaves() == 1)
            {
                const ulong nls = ulong(Decoder<Hoa2d, T>::getDecompositionOrder() + 1.);
                const T factor = 1. / (T)(nls);
                for(ulong i = 0; i <nls; i++)
                {
                    encoder.setAzimuth(T(i) * HOA_2PI / T(nls));
                    encoder.process(&factor, vector_harmonics);
                    vector_harmonics[0] = factor * 0.5;
                    Signal<T>::vector_add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix);
                }
            }
            else
            {
                T smallest_distance = HOA_2PI;
                vector<Planewave<Hoa2d, T> > channels;
                for(ulong i = 0; i < Decoder<Hoa2d, T>::getNumberOfPlanewaves(); i++)
                {
                    channels.push_back(Planewave<Hoa2d, T>(i, Math<T>::wrap_twopi(Decoder<Hoa2d, T>::getPlanewaveAzimuth(i)), 0.));
                }
                
                std::sort(channels.begin(), channels.end(), Planewave<Hoa2d, T>::sort_azimuth);
                
                {
                    const T current_angle   = channels[0].getAzimuth(0., 0., 0.);
                    const T previous_angle  = channels[channels.size() - 1].getAzimuth(0., 0., 0.);
                    const T previous_portion= (HOA_2PI - previous_angle) + current_angle;
                    if(smallest_distance > previous_portion)
                    {
                        smallest_distance = previous_portion;
                    }
                    //post("channel %i : %f", (int)channels[0].getIndex(), (float)(channels[0].getAzimuth() / HOA_2PI * 360.f));
                }
                for(ulong i = 1; i < channels.size(); i++)
                {
                    const T current_angle   = channels[i].getAzimuth(0., 0., 0.);
                    const T previous_angle  = channels[i-1].getAzimuth(0., 0., 0.);
                    const T previous_portion= current_angle - previous_angle;
                    if(smallest_distance > previous_portion)
                    {
                        smallest_distance = previous_portion;
                    }
                    //post("channel %i : %f", (int)channels[i].getIndex(), (float)(channels[i].getAzimuth() / HOA_2PI * 360.f));
                }
                //post("");
                
                if(smallest_distance > HOA_2PI / T(Decoder<Hoa2d, T>::getNumberOfHarmonics() + 1.))
                {
                    smallest_distance = HOA_2PI / T(Decoder<Hoa2d, T>::getNumberOfHarmonics() + 1.);
                }
                const ulong nvirtual = ceil(HOA_2PI / smallest_distance);
                const T factor = 1. / (T)(nvirtual);
                
                //post("number of virtual %i", nvirtual);
                for(ulong i = 0; i < nvirtual; i++)
                {
                    const T angle = T(i) / T(nvirtual) * HOA_2PI;
                    //post("virtual %i :  %f", (int)i , (float)(angle / HOA_2PI * 360.f));
                    if(angle < channels[0].getAzimuth(0., 0., 0.))
                    {
                        const T portion = (HOA_2PI - channels[channels.size()-1].getAzimuth(0., 0., 0.)) + channels[0].getAzimuth(0., 0., 0.);
                        
                        const T factor1 = (1. - ((channels[0].getAzimuth(0., 0., 0.) - angle) / portion)) * factor;
                        encoder.setAzimuth(angle);
                        encoder.process(&factor1, vector_harmonics);
                        vector_harmonics[0] = factor1 * 0.5;
                        Signal<T>::vector_add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[0].getIndex() * Decoder<Hoa2d, T>::getNumberOfHarmonics());
                        
                        const T factor2 = ((channels[0].getAzimuth(0., 0., 0.) - angle) / portion) * factor;
                        encoder.process(&factor2, vector_harmonics);
                        vector_harmonics[0] = factor2 * 0.5;
                        Signal<T>::vector_add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[channels.size() - 1].getIndex() * Decoder<Hoa2d, T>::getNumberOfHarmonics());
                        
                        //post("portion : %f", (float)portion / HOA_2PI * 360.f);
                        //post("channel %i (%f) : %f", (int)channels[channels.size()-1].getIndex(),
                        //     (float)(channels[channels.size()-1].getAzimuth() / HOA_2PI * 360.f),
                        //     (float)(factor2 / factor) * 100.f);
                        //post("channel %i (%f) : %f", (int)channels[0].getIndex(),
                        //     (float)(channels[0].getAzimuth() / HOA_2PI * 360.f),
                        //     (float)(factor1 / factor) * 100.f);
                    }
                    else if(angle >= channels[channels.size() - 1].getAzimuth(0., 0., 0.))
                    {
                        const T portion = (HOA_2PI - channels[channels.size()-1].getAzimuth(0., 0., 0.)) + channels[0].getAzimuth(0., 0., 0.);
                        
                        const T factor1 = (1. - ((angle - channels[channels.size()-1].getAzimuth(0., 0., 0.)) / portion)) * factor;
                        encoder.setAzimuth(angle);
                        encoder.process(&factor1, vector_harmonics);
                        vector_harmonics[0] = factor1 * 0.5;
                        Signal<T>::vector_add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[channels.size()-1].getIndex() * Decoder<Hoa2d, T>::getNumberOfHarmonics());
                        
                        const T factor2 = ((angle - channels[channels.size()-1].getAzimuth(0., 0., 0.)) / portion) * factor;
                        encoder.process(&factor2, vector_harmonics);
                        vector_harmonics[0] = factor2 * 0.5;
                        Signal<T>::vector_add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[0].getIndex() * Decoder<Hoa2d, T>::getNumberOfHarmonics());
                        
                        //post("portion : %f", (float)portion / HOA_2PI * 360.f);
                        //post("channel %i (%f) : %f", (int)channels[channels.size()-1].getIndex(),
                        //     (float)(channels[channels.size()-1].getAzimuth() / HOA_2PI * 360.f),
                        //     (float)(factor1 / factor) * 100.f);
                        //post("channel %i (%f) : %f", (int)channels[0].getIndex(),
                        //     (float)(channels[0].getAzimuth() / HOA_2PI * 360.f),
                        //     (float)(factor2 / factor) * 100.f);
                    }
                    else
                    {
                        for(ulong j = 1; j < channels.size(); j++)
                        {
                            if(angle < channels[j].getAzimuth(0., 0., 0.) && angle >= channels[j-1].getAzimuth(0., 0., 0.))
                            {
                                const T portion = (channels[j].getAzimuth(0., 0., 0.) - channels[j-1].getAzimuth(0., 0., 0.));
                                
                                const T factor1 = (1. - ((channels[j].getAzimuth(0., 0., 0.) - angle) / portion)) * factor;
                                encoder.setAzimuth(angle);
                                encoder.process(&factor1, vector_harmonics);
                                vector_harmonics[0] = factor1 * 0.5;
                                Signal<T>::vector_add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[j].getIndex() * Decoder<Hoa2d, T>::getNumberOfHarmonics());
                                
                                const T factor2 = ((channels[j].getAzimuth(0., 0., 0.) - angle) / portion) * factor;
                                encoder.process(&factor2, vector_harmonics);
                                vector_harmonics[0] = factor2 * 0.5;
                                Signal<T>::vector_add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[j-1].getIndex() * Decoder<Hoa2d, T>::getNumberOfHarmonics());
                                
                                //post("portion : %f", (float)portion / HOA_2PI * 360.f);
                                //post("channel %i (%f) : %f", (int)channels[j-1].getIndex(),
                                //     (float)(channels[j-1].getAzimuth() / HOA_2PI * 360.f),
                                //     (float)(factor2 / factor) * 100.f);
                                //post("channel %i (%f) : %f", (int)channels[j].getIndex(),
                                //     (float)(channels[j].getAzimuth() / HOA_2PI * 360.f),
                                //     (float)(factor1 / factor) * 100.f);
                                
                                break;
                            }
                        }
                    }
                    //post("");
                }
                channels.clear();
            }
        }
    };
    
#define HOA_NBIN_I 512
#define HOA_NBIN_H 11

    template <typename T> class Decoder<Hoa2d, T>::Binaural : public Decoder<Hoa2d, T>
    {
    private:
        ulong       m_vector_size;
        ulong       m_counter;
        T*          m_inputs;
        T*          m_results;
        T*          m_result_matrix_left;
        T*          m_result_matrix_right;
        T*          m_linear_vector_left;
        T*          m_linear_vector_right;
        T*          m_output_left;
        T*          m_output_right;
    public:
        
        //! The binaural decoder constructor.
        /**	The binaural decoder constructor allocates and initialize the member values to the decoding matrix depending on a order of decomposition and a number of channels. The order and the number of channels must be at least 1.
         @param     order				The order
         */
        Binaural(const ulong order) : Decoder<Hoa2d, T>(order, 2),
        m_vector_size(0ul), m_counter(0ul), m_inputs(nullptr), m_results(nullptr), m_linear_vector_left(nullptr), m_linear_vector_right(nullptr), m_output_left(nullptr), m_output_right(nullptr)
        {
            Decoder<Hoa2d, T>::setPlanewaveAzimuth(0, HOA_PI2*3.);
            Decoder<Hoa2d, T>::setPlanewaveAzimuth(1, HOA_PI2);
        }
        
        //! The binaural decoder destructor.
        /**	The binaural decoder destructor free the memory.
         */
        ~Binaural()
        {
            if(m_inputs)
                delete [] m_inputs;
            if(m_results)
                delete [] m_results;
            if(m_linear_vector_left)
                delete [] m_linear_vector_left;
            if(m_linear_vector_right)
                delete [] m_linear_vector_right;
            if(m_output_left)
                delete [] m_output_left;
            if(m_output_right)
                delete [] m_output_right;
        }
        
        //! This method computes the decoding matrice.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        void computeRendering(const ulong vectorsize = 64)  override
        {
            m_counter     = 0;
            m_vector_size = vectorsize;
            if(m_inputs)
                delete [] m_inputs;
            if(m_results)
                delete [] m_results;
            if(m_linear_vector_left)
                delete [] m_linear_vector_left;
            if(m_linear_vector_right)
                delete [] m_linear_vector_right;
            if(m_output_left)
                delete [] m_output_left;
            if(m_output_right)
                delete [] m_output_right;
            
            m_output_left = new T[m_vector_size];
            m_output_right = new T[m_vector_size];
            Signal<T>::vector_clear(m_vector_size, m_output_left);
            Signal<T>::vector_clear(m_vector_size, m_output_right);
            
            m_inputs  = new T[HOA_NBIN_H * m_vector_size];
            Signal<T>::vector_clear(HOA_NBIN_H * m_vector_size, m_inputs);
            
            m_results = new T[HOA_NBIN_I * 2 * m_vector_size];
            Signal<T>::vector_clear(HOA_NBIN_I * 2 * m_vector_size, m_results);
            m_result_matrix_left    = m_results;
            m_result_matrix_right   = m_results + m_vector_size  * HOA_NBIN_I;
            
            m_linear_vector_left    = new T[m_vector_size + HOA_NBIN_I - 1];
            Signal<T>::vector_clear(m_vector_size + HOA_NBIN_I - 1, m_linear_vector_left);
            m_linear_vector_right   = new T[m_vector_size + HOA_NBIN_I - 1];
            Signal<T>::vector_clear(m_vector_size + HOA_NBIN_I - 1, m_linear_vector_right);
        }
        
        //! This method performs the binaural decoding.
        /**	You should use this method for not-in-place processing and performs the binaural decoding on block of samples. The inputs matrix contains the spherical harmonics samples : inputs[number of harmonics][vector size] and the outputs matrix contains the headphones samples : outputs[2][vector size].
         @param     inputs	The input samples.
         @param     outputs  The output array that contains samples destinated to channels.
         */
        void processBlock() noexcept
        {
   
            Signal<T>::matrix_matrix_mul(HOA_NBIN_I * 2, m_vector_size, HOA_NBIN_H, Hrtf<Hoa2d, T>::getImpulse(), m_inputs, m_results);
            
            for(ulong i = 0; i < m_vector_size; i++)
            {
                Signal<T>::vector_add(HOA_NBIN_I, m_results + i, m_vector_size, m_linear_vector_left + i, 1);
                m_output_left[i] = m_linear_vector_left[i];
            }
            
            for(ulong i = 0; i < m_vector_size; i++)
            {
                Signal<T>::vector_add(HOA_NBIN_I, m_results + i + m_vector_size * HOA_NBIN_I, m_vector_size, m_linear_vector_right + i, 1);
                m_output_right[i] = m_linear_vector_right[i];
            }

            Signal<T>::vector_copy(HOA_NBIN_I - 1, m_linear_vector_left + m_vector_size, m_linear_vector_left);
            Signal<T>::vector_copy(HOA_NBIN_I - 1, m_linear_vector_right + m_vector_size, m_linear_vector_right);
            
            Signal<T>::vector_clear(m_vector_size, m_linear_vector_left + HOA_NBIN_I - 1);
            Signal<T>::vector_clear(m_vector_size, m_linear_vector_right + HOA_NBIN_I - 1);
        }
        
        //! This method performs the binaural decoding.
        /**	You should use this method for not-in-place processing and performs the binaural decoding sample by sample. The inputs array contains the spherical harmonics samples : inputs[number of harmonics] and the outputs array contains the headphones samples : outputs[2].
         
         @param     inputs	The input samples.
         @param     outputs  The output array that contains samples destinated to channels.
         */
        inline void process(const T* inputs, T* outputs) noexcept override
        {
            for(ulong i = 0; i < Decoder<Hoa2d, T>::getNumberOfHarmonics() && i < 9; i++)
            {
                m_inputs[i*m_vector_size+m_counter] = inputs[i];
            }
       
            outputs[0] = m_output_left[m_counter];
            outputs[1] = m_output_right[m_counter];
            if(++m_counter >= m_vector_size)
            {
                processBlock();
                m_counter = 0;
            }
        }
        
    };
    
#undef HOA_NBIN_I
#undef HOA_NBIN_H
    
    template <typename T> class Decoder<Hoa3d, T> : public Encoder<Hoa3d, T>::Basic, public Planewave<Hoa3d, T>::Processor
    {
    public:
        
        //! The decoder constructor.
        /**	The decoder constructor allocates and initialize the base classes.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        Decoder(const ulong order, const ulong numberOfPlanewaves) noexcept :
        Encoder<Hoa3d, T>::Basic(order),
        Planewave<Hoa3d, T>::Processor(numberOfPlanewaves)
        {
            ;
        }
        
        //! The destructor.
        /** The destructor free the memory.
         */
        virtual ~Decoder()
        {
            ;
        }
        
        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and performs the decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimym size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        inline virtual void process(const T* inputs, T* outputs) noexcept = 0;
        
        //! This method computes the decoding matrice.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        virtual void computeRendering(const ulong vectorsize = 64) = 0;
        
        //! The ambisonic regular decoder.
        /** The regular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if more or equal to the number of harmonics plus one and when the loudspeakers are equally spaced.
         */
        class Regular;
        
        //! The ambisonic irregular decoder.
        /** The irregular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if less than the number of harmonics plus one or when the loudspeakers are not equally spaced.
         */
        class Irregular;
        
        //! The ambisonic binaural decoder.
        /** The binaural decoder should be used to decode an ambisonic sound field for headphones.
         */
        class Binaural;
    };
    
    template <typename T> class Decoder<Hoa3d, T>::Regular : public Decoder<Hoa3d, T>
    {
    private:
        T*  m_matrix;
    public:
        
        //! The regular constructor.
        /**	The regular constructor allocates and initialize the decoding matrix depending on a order of decomposition and a number of channels. The order must be at least 1 and the number of channels must be at least the number of harmonics.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        Regular(const ulong order, const ulong numberOfPlanewaves) noexcept : Decoder<Hoa3d, T>(order, numberOfPlanewaves)
        {
            m_matrix = new T[Planewave<Hoa3d, T>::Processor::getNumberOfPlanewaves()*Encoder<Hoa3d, T>::Basic::getNumberOfHarmonics()];
            computeRendering();
        }
        
        //! The regular destructor.
        /**	The regular destructor free the memory.
         */
        ~Regular()
        {
            delete [] m_matrix;
        }
        
        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and performs the regular decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimym size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to channels.
         */
        inline void process(const T* inputs, T* outputs) noexcept
        {
            Signal<T>::matrix_vector_mul(Encoder<Hoa3d, T>::Basic::getNumberOfHarmonics(), Planewave<Hoa3d, T>::Processor::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }
        
        void computeRendering(const ulong vectorsize = 64)
        {
            const T factor = 12.5 / (T)(Encoder<Hoa3d, T>::Basic::getNumberOfHarmonics());
            for(ulong i = 0; i < Planewave<Hoa3d, T>::Processor::getNumberOfPlanewaves(); i++)
            {
                Encoder<Hoa3d, T>::Basic::setAzimuth(Planewave<Hoa3d, T>::Processor::getPlanewaveAzimuth(i));
                Encoder<Hoa3d, T>::Basic::setElevation(Planewave<Hoa3d, T>::Processor::getPlanewaveElevation(i));
                Encoder<Hoa3d, T>::Basic::process(&factor, m_matrix + i * Encoder<Hoa3d, T>::Basic::getNumberOfHarmonics());
                for(ulong j = 0; j < Encoder<Hoa3d, T>::Basic::getNumberOfHarmonics(); j++)
                {
                    const ulong l = Encoder<Hoa3d, T>::Basic::getHarmonicDegree(j);
                    m_matrix[i * Encoder<Hoa3d, T>::Basic::getNumberOfHarmonics() + j] *= (2. * l + 1.) / (4. * HOA_PI);
                }
            }
        }
    };
    
#define HOA_NBIN_I 512
#define HOA_NBIN_H 16
    
    template <typename T> class Decoder<Hoa3d, T>::Binaural : public Decoder<Hoa3d, T>
    {
    private:
        ulong       m_vector_size;
        ulong       m_counter;
        T*          m_inputs;
        T*          m_results;
        T*          m_result_matrix_left;
        T*          m_result_matrix_right;
        T*          m_linear_vector_left;
        T*          m_linear_vector_right;
        T*          m_output_left;
        T*          m_output_right;
    public:
        
        //! The binaural decoder constructor.
        /**	The binaural decoder constructor allocates and initialize the member values to the decoding matrix depending on a order of decomposition and a number of channels. The order and the number of channels must be at least 1.
         @param     order				The order
         */
        Binaural(const ulong order) : Decoder<Hoa3d, T>(order, 2),
        m_vector_size(0ul), m_counter(0ul), m_inputs(nullptr), m_results(nullptr), m_linear_vector_left(nullptr), m_linear_vector_right(nullptr), m_output_left(nullptr), m_output_right(nullptr)
        {
            Decoder::setPlanewaveAzimuth(0, HOA_PI2*3);
            Decoder::setPlanewaveAzimuth(1, HOA_PI2);
        }
        
        //! The binaural decoder destructor.
        /**	The binaural decoder destructor free the memory.
         */
        ~Binaural()
        {
            if(m_inputs)
                delete [] m_inputs;
            if(m_results)
                delete [] m_results;
            if(m_linear_vector_left)
                delete [] m_linear_vector_left;
            if(m_linear_vector_right)
                delete [] m_linear_vector_right;
            if(m_output_left)
                delete [] m_output_left;
            if(m_output_right)
                delete [] m_output_right;
        }
        
        //! This method computes the decoding matrice.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        void computeRendering(const ulong vectorsize = 64)
        {
            m_counter     = 0;
            m_vector_size = vectorsize;
            if(m_inputs)
                delete [] m_inputs;
            if(m_results)
                delete [] m_results;
            if(m_linear_vector_left)
                delete [] m_linear_vector_left;
            if(m_linear_vector_right)
                delete [] m_linear_vector_right;
            if(m_output_left)
                delete [] m_output_left;
            if(m_output_right)
                delete [] m_output_right;
            
            m_output_left = new T[m_vector_size];
            m_output_right = new T[m_vector_size];
            Signal<T>::vector_clear(m_vector_size, m_output_left);
            Signal<T>::vector_clear(m_vector_size, m_output_right);
            
            m_inputs  = new T[HOA_NBIN_H * m_vector_size];
            Signal<T>::vector_clear(HOA_NBIN_H * m_vector_size, m_inputs);
            
            m_results = new T[HOA_NBIN_I * 2 * m_vector_size];
            Signal<T>::vector_clear(HOA_NBIN_I * 2 * m_vector_size, m_results);
            m_result_matrix_left    = m_results;
            m_result_matrix_right   = m_results + m_vector_size  * HOA_NBIN_I;
            
            m_linear_vector_left    = new T[m_vector_size + HOA_NBIN_I - 1];
            Signal<T>::vector_clear(m_vector_size + HOA_NBIN_I - 1, m_linear_vector_left);
            m_linear_vector_right   = new T[m_vector_size + HOA_NBIN_I - 1];
            Signal<T>::vector_clear(m_vector_size + HOA_NBIN_I - 1, m_linear_vector_right);
        }
        
        //! This method performs the binaural decoding.
        /**	You should use this method for not-in-place processing and performs the binaural decoding on block of samples. The inputs matrix contains the spherical harmonics samples : inputs[number of harmonics][vector size] and the outputs matrix contains the headphones samples : outputs[2][vector size].
         @param     inputs	The input samples.
         @param     outputs  The output array that contains samples destinated to channels.
         */
        void process() noexcept
        {
            Signal<T>::matrix_matrix_mul(HOA_NBIN_I * 2, m_vector_size, HOA_NBIN_H, Hrtf<Hoa3d, T>::getImpulse(), m_inputs, m_results);
            
            for(ulong i = 0; i < m_vector_size; i++)
            {
                Signal<T>::vector_add(HOA_NBIN_I, m_results + i, m_vector_size, m_linear_vector_left + i, 1);
                m_output_left[i] = m_linear_vector_left[i];
            }
            
            for(ulong i = 0; i < m_vector_size; i++)
            {
                Signal<T>::vector_add(HOA_NBIN_I, m_results + i + m_vector_size * HOA_NBIN_I, m_vector_size, m_linear_vector_right + i, 1);
                m_output_right[i] = m_linear_vector_right[i];
            }
            
            Signal<T>::vector_copy(HOA_NBIN_I - 1, m_linear_vector_left + m_vector_size, m_linear_vector_left);
            Signal<T>::vector_copy(HOA_NBIN_I - 1, m_linear_vector_right + m_vector_size, m_linear_vector_right);
            
            Signal<T>::vector_clear(m_vector_size, m_linear_vector_left + HOA_NBIN_I - 1);
            Signal<T>::vector_clear(m_vector_size, m_linear_vector_right + HOA_NBIN_I - 1);
        }
        
        //! This method performs the binaural decoding.
        /**	You should use this method for not-in-place processing and performs the binaural decoding sample by sample. The inputs array contains the spherical harmonics samples : inputs[number of harmonics] and the outputs array contains the headphones samples : outputs[2].
         
         @param     inputs	The input samples.
         @param     outputs  The output array that contains samples destinated to channels.
         */
        void process(const T* inputs, T* outputs) noexcept
        {
            
            for(ulong i = 0; i < Encoder<Hoa3d, T>::Basic::getNumberOfHarmonics() && i < HOA_NBIN_H; i++)
            {
                m_inputs[i*m_vector_size+m_counter] = inputs[i];
            }
            
            outputs[0] = m_output_left[m_counter] * 0.05;
            outputs[1] = m_output_right[m_counter] * 0.05;
            if(++m_counter >= m_vector_size)
            {
                process();
                m_counter = 0;
            }
        }
        
    };
    
#undef HOA_NBIN_I
#undef HOA_NBIN_H
    
#endif

}

#endif


