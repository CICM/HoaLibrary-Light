
/*
 // Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#ifndef DEF_HOA_HRIR_LIGHT
#define DEF_HOA_HRIR_LIGHT

#include "HrirIrc1002C.hpp"

namespace hoa
{
    //! The hrir class gives the impulse responses to decode in the binaural mode.
    /** The hrir class gives the impulse responses to decode in the binaural mode in 2d / 3d / simple and double precision.
     */
    template <Dimension D, typename T> class Hrir
    {
    public:
        //! Get the impulse response of the HRTFs
        /** The impulse response may be 2d/3d with simple /double precision.
         */
        virtual const float* getImpulse() noexcept = 0;
    };

    template<> class Hrir <Hoa2d, float>
    {
    public:
        
        //! Gets the order of decomposition used to compute the matrix.
        /** Gets the order of decomposition used to compute the matrix.
         @return The order of decomposition.
         */
        static ulong getOrderOfDecomposition() noexcept
        {
            return 5ul;
        }
        
        //! Gets the number rows of the matrix (or the harmonics used to compute the matrix).
        /** Gets the number rows of the matrix (or the harmonics used to compute the matrix).
         @return The number rows of the matrix.
         */
        static ulong getNumberOfRows() noexcept
        {
            return 11ul;
        }
        
        //! Gets the number columns of the matrix (or the size of the responses used to compute the matrix).
        /** Gets the number columns of the matrix (or the size of the responses used to compute the matrix).
         @return The number columns of the matrix.
         */
        static ulong getNumberOfcolumns() noexcept
        {
            return 256;
        }
        
        //! Get the HRIR matrix for the left ear.
        /** Get the HRIR matrix for the left ear.
         @return The HRIR matrix for the left ear.
         */
        static const float* getLeftMatrix() noexcept
        {
            return Irc1002C_left;
        }
        
        //! Get the HRIR matrix for the right ear.
        /** Get the HRIR matrix for the right ear.
         @return The HRIR matrix for the right ear.
         */
        static const float* getRightMatrix() noexcept
        {
            return Irc1002C_right;
        }
    };

    template<> class Hrir <Hoa2d, double>
    {
    public:
        //! Get the impulse response of the HRTFs
        /** Get the impulse response of the HRTFs
         */
        static const double* getImpulse() noexcept
        {
            return NULL;
        }
    };

    template<> class Hrir <Hoa3d, float>
    {
    public:
        //! Get the impulse response of the HRTFs
        /** Get the impulse response of the HRTFs
         */
        static const float* getImpulse() noexcept
        {
            return NULL;
        }
    };

    template<> class Hrir <Hoa3d, double>
    {
    public:
        //! Get the impulse response of the HRTFs
        /** Get the impulse response of the HRTFs
         */
        static const double* getImpulse() noexcept
        {
            return NULL;
        }
    };

}

#endif
