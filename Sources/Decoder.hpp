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
    /** The decoder should be used to decode an ambisonic sound field for a set a channels on a circle depending on a decomposition order or for headphones.
     */
    template <Dimension D, typename T> class Decoder;
    
    template <typename T> class Decoder<Hoa2d, T> : public Encoder<Hoa2d, T>, public Planewave<Hoa2d, T>::Processor
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
        Encoder<Hoa2d, T>(order),
        Planewave<Hoa2d, T>::Processor(numberOfPlanewaves)
        {
            m_matrix = new T[Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves()*Encoder<Hoa2d, T>::getNumberOfHarmonics()];
            computeMatrix();
        }
        
        //! The destructor.
        /** The destructor free the memory.
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
            Signal<T>::matrix_vector_mul(Encoder<Hoa2d, T>::getNumberOfHarmonics(), Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }
        
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



