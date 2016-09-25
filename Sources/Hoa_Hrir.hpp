
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

#ifndef DEF_HOA_HRIR_LIGHT
#define DEF_HOA_HRIR_LIGHT

#include "Hoa_HrirIrc1002C2D.hpp"
#include "Hoa_HrirIrc1002C3D.hpp"

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
        virtual const float* getImpulse() hoa_noexcept = 0;
    };

    template<> class Hrir <Hoa2d, float>
    {
    public:
        
        //! Gets the order of decomposition used to compute the matrices.
        /** Gets the order of decomposition used to compute the matrices.
         @return The order of decomposition.
         */
        static size_t getOrderOfDecomposition() hoa_noexcept
        {
            return 5ul;
        }
        
        //! Gets the number rows of the matrices (or the size of the responses used to compute the matrices).
        /** Gets the number rows of the matrices (or the size of the responses used to compute the matrices).
         @return The number rows of the matrices.
         */
        static size_t getNumberOfRows() hoa_noexcept
        {
            return 512ul;
        }
        
        //! Gets the number columns of the matrices (or the number of harmonics used to compute the matrices).
        /** Gets the number columns of the matrices (or the number of harmonics used to compute the matrices).
         @return The number columns of the matrices.
         */
        static size_t getNumberOfColumns() hoa_noexcept
        {
            return 11ul;
        }
        
        //! Gets the size of the matrices (rows * columns).
        /** Gets the size of the matrices (rows * columns).
         @return The nsize of the matrices.
         */
        static size_t getMatricesSize() hoa_noexcept
        {
            return 5632ul;
        }
        
        //! Get the HRIR matrix for the left ear.
        /** Get the HRIR matrix for the left ear.
         @return The HRIR matrix for the left ear.
         */
        static const float* getLeftMatrix() hoa_noexcept
        {
            return Irc1002C_float_2d_left;
        }
        
        //! Get the HRIR matrix for the right ear.
        /** Get the HRIR matrix for the right ear.
         @return The HRIR matrix for the right ear.
         */
        static const float* getRightMatrix() hoa_noexcept
        {
            return Irc1002C_float_2d_right;
        }
    };

    template<> class Hrir <Hoa2d, double>
    {
    public:
        
        //! Gets the order of decomposition used to compute the matrices.
        /** Gets the order of decomposition used to compute the matrices.
         @return The order of decomposition.
         */
        static size_t getOrderOfDecomposition() hoa_noexcept
        {
            return 5ul;
        }
        
        //! Gets the number rows of the matrices (or the size of the responses used to compute the matrices).
        /** Gets the number rows of the matrices (or the size of the responses used to compute the matrices).
         @return The number rows of the matrices.
         */
        static size_t getNumberOfRows() hoa_noexcept
        {
            return 512ul;
        }
        
        //! Gets the number columns of the matrices (or the number of harmonics used to compute the matrices).
        /** Gets the number columns of the matrices (or the number of harmonics used to compute the matrices).
         @return The number columns of the matrices.
         */
        static size_t getNumberOfColumns() hoa_noexcept
        {
            return 11ul;
        }
        
        //! Gets the size of the matrices (rows * columns).
        /** Gets the size of the matrices (rows * columns).
         @return The nsize of the matrices.
         */
        static size_t getMatricesSize() hoa_noexcept
        {
            return 5632ul;
        }
        
        //! Get the HRIR matrix for the left ear.
        /** Get the HRIR matrix for the left ear.
         @return The HRIR matrix for the left ear.
         */
        static const double* getLeftMatrix() hoa_noexcept
        {
            return Irc1002C_double_2d_left;
        }
        
        //! Get the HRIR matrix for the right ear.
        /** Get the HRIR matrix for the right ear.
         @return The HRIR matrix for the right ear.
         */
        static const double* getRightMatrix() hoa_noexcept
        {
            return Irc1002C_double_2d_right;
        }
    };

    template<> class Hrir <Hoa3d, float>
    {
    public:
        
        //! Gets the order of decomposition used to compute the matrices.
        /** Gets the order of decomposition used to compute the matrices.
         @return The order of decomposition.
         */
        static size_t getOrderOfDecomposition() hoa_noexcept
        {
            return 3ul;
        }
        
        //! Gets the number rows of the matrices (or the size of the responses used to compute the matrices).
        /** Gets the number rows of the matrices (or the size of the responses used to compute the matrices).
         @return The number rows of the matrices.
         */
        static size_t getNumberOfRows() hoa_noexcept
        {
            return 512ul;
        }
        
        //! Gets the number columns of the matrices (or the number of harmonics used to compute the matrices).
        /** Gets the number columns of the matrices (or the number of harmonics used to compute the matrices).
         @return The number columns of the matrices.
         */
        static size_t getNumberOfColumns() hoa_noexcept
        {
            return 16ul;
        }
        
        //! Gets the size of the matrices (rows * columns).
        /** Gets the size of the matrices (rows * columns).
         @return The nsize of the matrices.
         */
        static size_t getMatricesSize() hoa_noexcept
        {
            return 8192ul;
        }
        
        //! Get the HRIR matrix for the left ear.
        /** Get the HRIR matrix for the left ear.
         @return The HRIR matrix for the left ear.
         */
        static const float* getLeftMatrix() hoa_noexcept
        {
            return Irc1002C_float_3d_left;
        }
        
        //! Get the HRIR matrix for the right ear.
        /** Get the HRIR matrix for the right ear.
         @return The HRIR matrix for the right ear.
         */
        static const float* getRightMatrix() hoa_noexcept
        {
            return Irc1002C_float_3d_right;
        }
    };

    template<> class Hrir <Hoa3d, double>
    {
    public:
        
        //! Gets the order of decomposition used to compute the matrices.
        /** Gets the order of decomposition used to compute the matrices.
         @return The order of decomposition.
         */
        static size_t getOrderOfDecomposition() hoa_noexcept
        {
            return 3ul;
        }
        
        //! Gets the number rows of the matrices (or the size of the responses used to compute the matrices).
        /** Gets the number rows of the matrices (or the size of the responses used to compute the matrices).
         @return The number rows of the matrices.
         */
        static size_t getNumberOfRows() hoa_noexcept
        {
            return 512ul;
        }
        
        //! Gets the number columns of the matrices (or the number of harmonics used to compute the matrices).
        /** Gets the number columns of the matrices (or the number of harmonics used to compute the matrices).
         @return The number columns of the matrices.
         */
        static size_t getNumberOfColumns() hoa_noexcept
        {
            return 16ul;
        }
        
        //! Gets the size of the matrices (rows * columns).
        /** Gets the size of the matrices (rows * columns).
         @return The nsize of the matrices.
         */
        static size_t getMatricesSize() hoa_noexcept
        {
            return 8192ul;
        }
        
        //! Get the HRIR matrix for the left ear.
        /** Get the HRIR matrix for the left ear.
         @return The HRIR matrix for the left ear.
         */
        static const double* getLeftMatrix() hoa_noexcept
        {
            return Irc1002C_double_3d_left;
        }
        
        //! Get the HRIR matrix for the right ear.
        /** Get the HRIR matrix for the right ear.
         @return The HRIR matrix for the right ear.
         */
        static const double* getRightMatrix() hoa_noexcept
        {
            return Irc1002C_double_3d_right;
        }
    };

}

#endif
