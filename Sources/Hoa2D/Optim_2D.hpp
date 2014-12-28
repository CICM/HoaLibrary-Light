/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_OPTIM
#define DEF_HOA_2D_OPTIM

#include "Ambisonic_2D.hpp"

namespace hoa
{

    //! The ambisonic optimization.
    /** The optimization should be used to optimize the ambisonic sound field. There are 3 optimization modes, Basic (no optimization), MaxRe (energy vector optimization) and InPhase (energy and velocity vector optimization). Basic has no effect, it should be used with a perfect ambisonic channels, arrengement where all the channels are to equal distance on a circle, and for a listener placed at the perfect center of the circle. MaxRe should be used for auditory confined to the center of the circle. InPhase should be used when the auditory covers the entire channels area and when the channels arragement is not a perfect circle or when the channels are not to equal distance. Note that the optimizations decrease the precision of the sound field restitution thus it can be compared to particular cases of the fractional orders.
     */
    template <typename T> class Optim : public Ambisonic2D<T>
    {
    private:
        Hoa::Mode   m_mode;
        T*       m_harmonics;
        
    public:
        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
            @param     order	The order.
         */
        Optim(unsigned long order, Hoa::Mode mode = Hoa::Basic) noexcept : Ambisonic2D<T>(order)
        {
            m_harmonics = new T[Ambisonic<T>::m_number_of_harmonics];
            setMode(mode);
        }
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~Optim()
        {
            delete [] m_harmonics;
        }
        
        //! Retrieve the optimization mode.
        /** Retrieve the optimization mode : Basic, MaxRe or InPhase.
         
            @return The method returns the optimization mode.
         */
        inline Hoa::Mode getMode() const noexcept
        {
            return m_mode;
        }
        
        //! This method set the optimization mode.
        /**	The mode should be one of the 3 optimization modes, Basic, MaxRe or InPhase.
         
            @param     mode The optimization mode.
         */
        void setMode(Hoa::Mode mode) noexcept
        {
            m_mode = mode;
            if(m_mode == Hoa::Basic)
            {
                for(unsigned long i = 0; i < Ambisonic<T>::m_number_of_harmonics; i++)
                {
                    m_harmonics[i] = 1.;
                }
            }
            else if(m_mode == Hoa::MaxRe)
            {
                for(unsigned long i = 0; i < Ambisonic<T>::m_number_of_harmonics; i++)
                {
                    m_harmonics[i] = std::cos(fabs((double)Ambisonic<T>::getHarmonicDegree(i)) * HOA_PI / (double)(2. * Ambisonic<T>::m_order_of_decomposition + 2));
                }
            }
            else
            {
                for(unsigned long i = 0; i < Ambisonic<T>::m_number_of_harmonics; i++)
                {
                    const long double temp1 = (long double)factorial(Ambisonic<T>::m_order_of_decomposition) / (long double)factorial(Ambisonic<T>::m_order_of_decomposition + Ambisonic<T>::getHarmonicDegree(i) + 1.);
                    const long double temp2 = (long double)factorial(Ambisonic<T>::m_order_of_decomposition + 1.) / (long double)factorial(Ambisonic<T>::m_order_of_decomposition - fabs((double)Ambisonic<T>::getHarmonicDegree(i)));
                    m_harmonics[i] = temp1 * temp2;
                }
            }
        }
        
        //! This method performs the optimization with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the optimization sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) const noexcept
        {
            T const* factor = m_harmonics;
            for(unsigned long i = 0; i < Ambisonic<T>::m_number_of_harmonics; i++)
            {
                *(outputs++) = (*inputs++) * (*factor++);
            }
        }
    };
}

#endif



