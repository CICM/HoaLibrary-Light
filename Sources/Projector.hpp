/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_PROJECTOR_LIGHT
#define DEF_HOA_PROJECTOR_LIGHT

#include "Harmonics.hpp"
#include "Planewaves.hpp"

namespace hoa
{
    //! The ambisonic projector.
    /** The projector should be used to decode an ambisonic sound field for a set a channels on a circle depending on a decomposition order or for headphones.
     */
    template <Dimension D, typename T> class Projector;
    
    template <typename T> class Projector<Hoa2d, T> : public Encoder<Hoa2d, T>::Basic, public Planewave<Hoa2d, T>::Processor
    {
    private:
        T*  m_matrix;
    public:
        
        //! The regular constructor.
        /**	The regular constructor allocates and initialize the decoding matrix depending on a decomposition order and a number of channels. The order must be at least 1 and the number of channels must be at least the number of harmonics.
         @param     order				The order
         @param     numberOfPlanewaves     The number of planewaves.
         */
        Projector(const ulong order, const ulong numberOfPlanewaves) noexcept :
        Encoder<Hoa2d, T>::Basic(order),
        Planewave<Hoa2d, T>::Processor(numberOfPlanewaves)
        {
            m_matrix = new T[Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves() * Encoder<Hoa2d, T>::getNumberOfHarmonics()];
            const T factor = 1. / (T)(Encoder<Hoa2d, T>::getDecompositionOrder() + 1.);
            for(ulong i = 0; i < Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves(); i++)
            {
                Encoder<Hoa2d, T>::Basic::setAzimuth(Planewave<Hoa2d, T>::Processor::getPlanewaveAzimuth(i));
                Encoder<Hoa2d, T>::Basic::process(&factor, m_matrix + i * Encoder<Hoa2d, T>::getNumberOfHarmonics());
                m_matrix[i * Encoder<Hoa2d, T>::getNumberOfHarmonics()] = factor * 0.5;
            }
        }
        
        //! The Rotate destructor.
        /**	The Rotate destructor free the memory.
         */
        ~Projector()
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
    };
}

#endif



