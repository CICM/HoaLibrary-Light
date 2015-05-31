/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco, Thomas Le Meur & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_PROJECTOR_LIGHT
#define DEF_HOA_PROJECTOR_LIGHT

#include "Encoder.hpp"
#include "Planewaves.hpp"

namespace hoa
{
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    //! The ambisonic projector.
    /** The projector should be used to decode an ambisonic sound field for a set a channels on a circle depending on a order of decomposition or for headphones.
     */
    template <Dimension D, typename T> class Projector;

    template <typename T> class Projector<Hoa2d, T> : public Encoder<Hoa2d, T>::Basic, public Processor<Hoa2d, T>::Planewaves
    {
    private:
        T*  m_matrix;
    public:

        //! The regular constructor.
        /**	The regular constructor allocates and initialize the decoding matrix depending on a order of decomposition and a number of channels. The order must be at least 1 and the number of channels must be at least the number of harmonics.
         @param     order				The order
         @param     numberOfPlanewaves     The number of planewaves.
         */
        Projector(const ulong order, const ulong numberOfPlanewaves) noexcept :
        Encoder<Hoa2d, T>::Basic(order),
        Processor<Hoa2d, T>::Planewaves(numberOfPlanewaves)
        {
            m_matrix = new T[Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves() * Encoder<Hoa2d, T>::getNumberOfHarmonics()];
            const T factor = 1. / (T)(Encoder<Hoa2d, T>::getDecompositionOrder() + 1.);
            for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                Encoder<Hoa2d, T>::Basic::setAzimuth(Processor<Hoa2d, T>::Planewaves::getPlanewaveAzimuth(i));
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
        /**	You should use this method for in-place or not-in-place processing and performs the regular decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimum size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to channels.
         */
        inline void process(const T* inputs, T* outputs) noexcept override
        {
            Signal<T>::mul(Encoder<Hoa2d, T>::getNumberOfHarmonics(), Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }
    };
#endif
}

#endif



