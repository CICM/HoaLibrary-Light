/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_DECODER_LIGHT
#define DEF_HOA_DECODER_LIGHT

#include "Encoder.hpp"

namespace hoa
{
    //! The ambisonic decoder.
    /** The decoder is a virtual class from witch inherit the regular, irregular and binaural decoders.
     */
    template <Dimension D, typename T> class Decoder;
    
    template <typename T> class Decoder<Hoa2d, T> : public Encoder<Hoa2d, T>, public Planewave<Hoa2d, T>::Processor
    {
    public:
        
        //! The decoder constructor.
        /**	The decoder constructor allocates and initialize the base classes.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        Decoder(const ulong order, const ulong numberOfPlanewaves) noexcept :
        Encoder<Hoa2d, T>(order),
        Planewave<Hoa2d, T>::Processor(numberOfPlanewaves)
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
        inline virtual void process(const T* inputs, T* outputs) const noexcept = 0;
        
        //! This method computes the decoding matrice.
        /**	You should use this method after changing the position of the loudspeakers.
         */
        virtual void computeMatrix() noexcept = 0;
        
        //! The ambisonic regular decoder.
        /** The regular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if more or equal to the number of harmonics plus one and when the loudspeakers are equally spaced.
         */
        class Regular;
        
        //! The ambisonic irregular decoder.
        /** The irregular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if less than the number of harmonics plus one or when the loudspeakers are not equally spaced.
         */
       class Irregular;
    };
    
    template <typename T> class Decoder<Hoa2d, T>::Regular : public Decoder<Hoa2d, T>
    {
    private:
        T*  m_matrix;
    public:
        
        //! The regular constructor.
        /**	The regular constructor allocates and initialize the decoding matrix depending of a decomposition order and a number of channels. The order must be at least 1 and the number of channels must be at least the number of harmonics.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        Regular(const ulong order, const ulong numberOfPlanewaves) noexcept : Decoder<Hoa2d, T>(order, numberOfPlanewaves)
        {
            m_matrix = new T[Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves()*Encoder<Hoa2d, T>::getNumberOfHarmonics()];
            computeMatrix();
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
        inline void process(const T* inputs, T* outputs) const noexcept
        {
            Signal<T>::matrix_vector_mul(Encoder<Hoa2d, T>::getNumberOfHarmonics(), Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }
        
        //! This method computes the decoding matrice.
        /**	You should use this method after changing the position of the loudspeakers.
         */
        void computeMatrix() noexcept
        {
            const T factor = 1. / (T)(Encoder<Hoa2d, T>::getDecompositionOrder() + 1.);
            for(ulong i = 0; i < Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves(); i++)
            {
                Encoder<Hoa2d, T>::setAzimuth(Planewave<Hoa2d, T>::Processor::getPlanewaveAzimuthRotated(i));
                Encoder<Hoa2d, T>::process(&factor, m_matrix + i * Encoder<Hoa2d, T>::getNumberOfHarmonics());
                m_matrix[i * Encoder<Hoa2d, T>::getNumberOfHarmonics()] = factor * 0.5;
            }
        }
    };
    
    template <typename T> class Decoder<Hoa2d, T>::Irregular : public Decoder<Hoa2d, T>
    {
    private:
        T*  m_matrix;
    public:
        
        //! The irregular constructor.
        /**	The irregular constructor allocates and initialize the decoding matrix depending of a decomposition order and a number of channels. The order must be at least 1.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        Irregular(const ulong order, const ulong numberOfPlanewaves) noexcept : Decoder<Hoa2d, T>(order, numberOfPlanewaves)
        {
            m_matrix = new T[Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves()*Encoder<Hoa2d, T>::getNumberOfHarmonics()];
            computeMatrix();
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
        inline void process(const T* inputs, T* outputs) const noexcept
        {
            Signal<T>::matrix_vector_mul(Encoder<Hoa2d, T>::getNumberOfHarmonics(), Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }
        
        //! This method computes the decoding matrice.
        /**	You should use this method after changing the position of the loudspeakers.
         */
        void computeMatrix() noexcept
        {
            Signal<T>::vector_clear(Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves()*Encoder<Hoa2d, T>::getNumberOfHarmonics(), m_matrix);
            T vector_harmonics[Encoder<Hoa2d, T>::getNumberOfHarmonics()];
            
            if(Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves() == 1)
            {
                const ulong nls = ulong(Encoder<Hoa2d, T>::getDecompositionOrder() + 1.);
                const T factor = 1. / (T)(nls);
                for(ulong i = 0; i <nls; i++)
                {
                    Encoder<Hoa2d, T>::setAzimuth(T(i) * HOA_2PI / T(nls));
                    Encoder<Hoa2d, T>::process(&factor, vector_harmonics);
                    vector_harmonics[0] = factor * 0.5;
                    Signal<T>::vector_add(Encoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix);
                }
            }
            else
            {
                T smallest_distance = HOA_2PI;
                vector<Planewave<Hoa2d, T> > channels;
                for(ulong i = 0; i < Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves(); i++)
                {
                    channels.push_back(Planewave<Hoa2d, T>(i, Math<T>::wrap_twopi(Planewave<Hoa2d, T>::Processor::getPlanewaveAzimuth(i) + Planewave<Hoa2d, T>::Processor::getPlanewavesRotation())));
                }
                
                std::sort(channels.begin(), channels.end(), Planewave<Hoa2d, T>::sort_azimuth);
                
                {
                    const T current_angle   = channels[0].getAzimuth();
                    const T previous_angle  = channels[channels.size() - 1].getAzimuth();
                    const T previous_portion= (HOA_2PI - previous_angle) + current_angle;
                    if(smallest_distance > previous_portion)
                    {
                        smallest_distance = previous_portion;
                    }
                    //post("channel %i : %f", (int)channels[0].getIndex(), (float)(channels[0].getAzimuth() / HOA_2PI * 360.f));
                }
                for(ulong i = 1; i < channels.size(); i++)
                {
                    const T current_angle   = channels[i].getAzimuth();
                    const T previous_angle  = channels[i-1].getAzimuth();
                    const T previous_portion= current_angle - previous_angle;
                    if(smallest_distance > previous_portion)
                    {
                        smallest_distance = previous_portion;
                    }
                    //post("channel %i : %f", (int)channels[i].getIndex(), (float)(channels[i].getAzimuth() / HOA_2PI * 360.f));
                }
                //post("");
                
                if(smallest_distance > HOA_2PI / T(Encoder<Hoa2d, T>::getNumberOfHarmonics() + 1.))
                {
                    smallest_distance = HOA_2PI / T(Encoder<Hoa2d, T>::getNumberOfHarmonics() + 1.);
                }
                const ulong nvirtual = ceil(HOA_2PI / smallest_distance);
                const T factor = 1. / (T)(nvirtual);
                
                //post("number of virtual %i", nvirtual);
                for(ulong i = 0; i < nvirtual; i++)
                {
                    const T angle = T(i) / T(nvirtual) * HOA_2PI;
                    //post("virtual %i :  %f", (int)i , (float)(angle / HOA_2PI * 360.f));
                    if(angle < channels[0].getAzimuth())
                    {
                        const T portion = (HOA_2PI - channels[channels.size()-1].getAzimuth()) + channels[0].getAzimuth();
                        
                        const T factor1 = (1. - ((channels[0].getAzimuth() - angle) / portion)) * factor;
                        Encoder<Hoa2d, T>::setAzimuth(angle);
                        Encoder<Hoa2d, T>::process(&factor1, vector_harmonics);
                        vector_harmonics[0] = factor1 * 0.5;
                        Signal<T>::vector_add(Encoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[0].getIndex() * Encoder<Hoa2d, T>::getNumberOfHarmonics());
                        
                        const T factor2 = ((channels[0].getAzimuth() - angle) / portion) * factor;
                        Encoder<Hoa2d, T>::process(&factor2, vector_harmonics);
                        vector_harmonics[0] = factor2 * 0.5;
                        Signal<T>::vector_add(Encoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[channels.size() - 1].getIndex() * Encoder<Hoa2d, T>::getNumberOfHarmonics());
                        
                        //post("portion : %f", (float)portion / HOA_2PI * 360.f);
                        //post("channel %i (%f) : %f", (int)channels[channels.size()-1].getIndex(),
                        //     (float)(channels[channels.size()-1].getAzimuth() / HOA_2PI * 360.f),
                        //     (float)(factor2 / factor) * 100.f);
                        //post("channel %i (%f) : %f", (int)channels[0].getIndex(),
                        //     (float)(channels[0].getAzimuth() / HOA_2PI * 360.f),
                        //     (float)(factor1 / factor) * 100.f);
                    }
                    else if(angle >= channels[channels.size() - 1].getAzimuth())
                    {
                        const T portion = (HOA_2PI - channels[channels.size()-1].getAzimuth()) + channels[0].getAzimuth();
                        
                        const T factor1 = (1. - ((angle - channels[channels.size()-1].getAzimuth()) / portion)) * factor;
                        Encoder<Hoa2d, T>::setAzimuth(angle);
                        Encoder<Hoa2d, T>::process(&factor1, vector_harmonics);
                        vector_harmonics[0] = factor1 * 0.5;
                        Signal<T>::vector_add(Encoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[channels.size()-1].getIndex() * Encoder<Hoa2d, T>::getNumberOfHarmonics());
                        
                        const T factor2 = ((angle - channels[channels.size()-1].getAzimuth()) / portion) * factor;
                        Encoder<Hoa2d, T>::process(&factor2, vector_harmonics);
                        vector_harmonics[0] = factor2 * 0.5;
                        Signal<T>::vector_add(Encoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[0].getIndex() * Encoder<Hoa2d, T>::getNumberOfHarmonics());
                        
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
                            if(angle < channels[j].getAzimuth() && angle >= channels[j-1].getAzimuth())
                            {
                                const T portion = (channels[j].getAzimuth() - channels[j-1].getAzimuth());
                                
                                const T factor1 = (1. - ((channels[j].getAzimuth() - angle) / portion)) * factor;
                                Encoder<Hoa2d, T>::setAzimuth(angle);
                                Encoder<Hoa2d, T>::process(&factor1, vector_harmonics);
                                vector_harmonics[0] = factor1 * 0.5;
                                Signal<T>::vector_add(Encoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[j].getIndex() * Encoder<Hoa2d, T>::getNumberOfHarmonics());
                                
                                const T factor2 = ((channels[j].getAzimuth() - angle) / portion) * factor;
                                Encoder<Hoa2d, T>::process(&factor2, vector_harmonics);
                                vector_harmonics[0] = factor2 * 0.5;
                                Signal<T>::vector_add(Encoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[j-1].getIndex() * Encoder<Hoa2d, T>::getNumberOfHarmonics());
                                
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
    
    template <typename T> class Decoder<Hoa3d, T> : public Encoder<Hoa3d, T>, public Planewave<Hoa3d, T>::Processor
    {
    private:
        T*  m_matrix;
    public:
        
        //! The regular constructor.
        /**	The regular constructor allocates and initialize the decoding matrix depending of a decomposition order and a number of channels. The order must be at least 1 and the number of channels must be at least the number of harmonics.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        Decoder(const ulong order, const ulong numberOfPlanewaves) noexcept :
        Encoder<Hoa3d, T>(order),
        Planewave<Hoa3d, T>::Processor(numberOfPlanewaves)
        {
            m_matrix = new T[Planewave<Hoa3d, T>::Processor::getNumberOfPlanewaves()*Encoder<Hoa3d, T>::getNumberOfHarmonics()];
            computeMatrix();
        }
        
        //! The Rotate destructor.
        /**	The Rotate destructor free the memory.
         */
        ~Decoder()
        {
            delete [] m_matrix;
        }
        
        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and performs the regular decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimym size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to channels.
         */
        inline void process(const T* inputs, T* outputs) const noexcept
        {
            Signal<T>::matrix_vector_mul(Encoder<Hoa3d, T>::getNumberOfHarmonics(), Planewave<Hoa3d, T>::Processor::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }
        
        void computeMatrix() noexcept
        {
            const T factor = 12.5 / (T)(Encoder<Hoa3d, T>::getNumberOfHarmonics());
            for(ulong i = 0; i < Planewave<Hoa3d, T>::Processor::getNumberOfPlanewaves(); i++)
            {
                Encoder<Hoa3d, T>::setAzimuth(Planewave<Hoa3d, T>::Processor::getPlanewaveAzimuthRotated(i));
                Encoder<Hoa3d, T>::setElevation(Planewave<Hoa3d, T>::Processor::getPlanewaveElevationRotated(i));
                Encoder<Hoa3d, T>::process(&factor, m_matrix + i * Encoder<Hoa3d, T>::getNumberOfHarmonics());
                for(ulong j = 0; j < Encoder<Hoa3d, T>::getNumberOfHarmonics(); j++)
                {
                    const ulong l = Encoder<Hoa3d, T>::getHarmonicDegree(j);
                    m_matrix[i * Encoder<Hoa3d, T>::getNumberOfHarmonics() + j] *= (2. * l + 1.) / (4. * HOA_PI);
                }
            }
        }
    };
}

#endif



