/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
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
    template <Dimension D, typename T> class Decoder : public Encoder<D, T>, public Planewave<D, T>::Processor
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
        Encoder<D, T>(order),
        Planewave<D, T>::Processor(numberOfPlanewaves)
        {
            m_matrix = new T[Planewave<D, T>::Processor::getNumberOfPlanewaves()*Encoder<Hoa2d, T>::getNumberOfHarmonics()];
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
            Signal<T>::matrix_vector_mul(Encoder<Hoa2d, T>::getNumberOfHarmonics(), Planewave<D, T>::Processor::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }
        
        void computeMatrix() noexcept
        {
            const ulong nplanewaves = Planewave<D, T>::Processor::getNumberOfPlanewaves();
            const ulong nharmonics  = Encoder<Hoa2d, T>::getNumberOfHarmonics();
            if(nplanewaves == 1)
            {
                const T factor = 1. / (Encoder<D, T>::getDecompositionOrder() + 1.);
                for(ulong i = 0; i < nharmonics + 1; i++)
                {
                    Encoder<D, T>::setAzimuth((T)i / (T)(nharmonics + 1) * HOA_2PI);
                    Encoder<D, T>::process(&factor, m_matrix+i*(nharmonics));
                    m_matrix[i * nharmonics] = factor * 0.5f;
                }
            }
            else
            {
                T* vharmonics = new T[nharmonics];
                for(ulong i = 0; i < nharmonics * nplanewaves; i++)
                {
                    m_matrix[i] = 0.;
                }
                
                // Retrieve the planewaves and sort them
                vector<Planewave<Hoa2d, T>> planewaves;
                for(ulong i = 0; i < nplanewaves; i++)
                {
                    planewaves.push_back(Planewave<D, T>(i, Planewave<D, T>::Processor::getPlanewaveAzimuth(i)));
                }
                sort(planewaves.begin(), planewaves.end());
                
                // Retrieve the smallest angle and if the configuration is ambisonics at +/- 1°
                bool ambisonic = nplanewaves >= nharmonics;
                T refdist = HOA_2PI - planewaves[nplanewaves-1].getAzimuth() + planewaves[0].getAzimuth();
                T dist = refdist;
                for(ulong i = 1; i < nplanewaves; i++)
                {
                    T ndist = planewaves[i].getAzimuth() - planewaves[i-1].getAzimuth();
                    if(ndist > refdist + HOA_1DEG || ndist < refdist - HOA_1DEG)
                    {
                        ambisonic = false;
                    }
                    dist = min(ndist, dist);
                }
                dist = max(dist, (T)HOA_1DEG);
                
                //  Compute the matrix for the ambisonic decoding
                if(ambisonic)
                {
                    const T factor = 1. / (T)(Encoder<D, T>::getDecompositionOrder() + 1.);
                    for(ulong i = 0; i < planewaves.size(); i++)
                    {
                        Encoder<D, T>::setAzimuth(planewaves[i].getAzimuth() + Planewave<D, T>::Processor::getPlanewavesOffset());
                        Encoder<D, T>::process(&factor, m_matrix+planewaves[i].getIndex()*(nharmonics));
                        m_matrix[planewaves[i].getIndex() * nharmonics] = factor * 0.5;
                    }
                }
                // Compute the matrix for the irregular decoding
                else
                {
                    const ulong nvirtual= max((ulong)(HOA_2PI / dist), nharmonics + 1);
                    for(ulong i = 0; i < nvirtual; i++)
                    {
                        ulong  channel_index1 = 0, channel_index2 = 0;
                        T distance1 = HOA_2PI, distance2 = HOA_2PI;
                        T angle = (T)i / (T)nvirtual * HOA_2PI;
                        T factor_index1 = 0., factor_index2 = 0.;
                        
                        // Find the closest real planewave to the virtual one
                        for(ulong j = 0; j < nplanewaves; j++)
                        {
                            T distref = Math<T>::wrap_twopi(planewaves[j].getAzimuth() - angle);
                            if(distref > HOA_PI)
                            {
                                distref = HOA_2PI - distref;
                            }
                            if(distref < distance1)
                            {
                                distance1 = distref;
                                channel_index1 = j;
                            }
                        }
                        if(Math<T>::wrap_pi(planewaves[channel_index1].getAzimuth() - angle) >= 0)
                        {
                            channel_index2 = Math<T>::wrap_ptr(channel_index1 - 1, nplanewaves);
                        }
                        else
                        {
                            channel_index2 = Math<T>::wrap_ptr(channel_index1 + 1, nplanewaves);
                        }
                        distance2  = Math<T>::wrap_twopi(planewaves[channel_index2].getAzimuth() - angle);
                        if(distance2 > HOA_PI)
                        {
                            distance2 = HOA_2PI - distance2;
                        }
                        
                        if(distance1 - distance2 < HOA_2PI / (T)nvirtual)
                        {
                            factor_index1   = cos(distance1 / (distance1 + distance2) * HOA_PI2);
                            factor_index2   = cos(distance2 / (distance1 + distance2) * HOA_PI2);
                        }
                        else
                        {
                            factor_index1   = 1;
                            factor_index2   = 0;
                        }

                        const T factor = 1. / (T)(Encoder<D, T>::getDecompositionOrder() + 1.);
                        Encoder<D, T>::process(&factor, vharmonics);
                        Encoder<D, T>::setAzimuth(angle + Planewave<D, T>::Processor::getPlanewavesOffset());
                        
                        m_matrix[planewaves[channel_index1].getIndex() * nharmonics] += factor * 0.5 * factor_index1;
                        m_matrix[planewaves[channel_index2].getIndex() * nharmonics] += factor * 0.5 * factor_index2;
                        
                        for(ulong j = 1; j < nharmonics; j++)
                        {
                            m_matrix[channel_index1 * nharmonics + j] += vharmonics[j] * factor_index1;
                            m_matrix[channel_index2 * nharmonics + j] += vharmonics[j] * factor_index2;
                        }
                    }
                }
                delete [] vharmonics;
            }
        }
    };
}

#endif



