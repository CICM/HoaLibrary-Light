/*
// Copyright (c) 2012-2016 CICM - Universite Paris 8 - Labex Arts H2H.
// Authors :
// 2012: Pierre Guillot, Eliott Paris & Julien Colafrancesco.
// 2012-2015: Pierre Guillot & Eliott Paris.
// 2015: Pierre Guillot & Eliott Paris & Thomas Le Meur (Light version)
// 2016: Pierre Guillot & Eliott Paris.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_PROJECTOR_LIGHT
#define DEF_HOA_PROJECTOR_LIGHT

#include "Hoa_Encoder.hpp"
#include "Hoa_Planewaves.hpp"

namespace hoa
{
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    //! The ambisonic projector.
    /** The projector should be used to decode an ambisonic sound field for a set a channels on a circle depending on a order of decomposition or for headphones.
     */
    template <Dimension D, typename T> class Projector;

    template <typename T> class Projector<Hoa2d, T> : public EncoderBasic<Hoa2d, T>, public ProcessorPlanewaves<Hoa2d, T>
    {
    private:
        T*  m_matrix;
    public:

        //! The regular constructor.
        /**	The regular constructor allocates and initialize the decoding matrix depending on a order of decomposition and a number of channels. The order must be at least 1 and the number of channels must be at least the number of harmonics.
         @param     order				The order
         @param     numberOfPlanewaves     The number of planewaves.
         */
        Projector(const size_t order, const size_t numberOfPlanewaves) hoa_noexcept :
        EncoderBasic<Hoa2d, T>(order),
        ProcessorPlanewaves<Hoa2d, T>(numberOfPlanewaves)
        {
            m_matrix = Signal<T>::alloc(ProcessorPlanewaves<Hoa2d, T>::getNumberOfPlanewaves() * Encoder<Hoa2d, T>::getNumberOfHarmonics());
            const T factor = 1. / (T)(Encoder<Hoa2d, T>::getDecompositionOrder() + 1.);
            for(size_t i = 0; i < ProcessorPlanewaves<Hoa2d, T>::getNumberOfPlanewaves(); i++)
            {
                EncoderBasic<Hoa2d, T>::setAzimuth(ProcessorPlanewaves<Hoa2d, T>::getPlanewaveAzimuth(i));
                EncoderBasic<Hoa2d, T>::process(&factor, m_matrix + i * Encoder<Hoa2d, T>::getNumberOfHarmonics());
                m_matrix[i * Encoder<Hoa2d, T>::getNumberOfHarmonics()] = factor * 0.5;
            }
        }

        //! The Rotate destructor.
        /**	The Rotate destructor free the memory.
         */
        ~Projector()
        {
            Signal<T>::free(m_matrix);
        }

        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and performs the regular decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimum size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to channels.
         */
        inline void process(const T* inputs, T* outputs) hoa_noexcept hoa_override
        {
            Signal<T>::mul(Encoder<Hoa2d, T>::getNumberOfHarmonics(), ProcessorPlanewaves<Hoa2d, T>::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }
    };
#endif
}

#endif



