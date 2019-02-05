
/*
 // Copyright (c) 2012-2017 CICM - Universite Paris 8 - Labex Arts H2H.
// Authors :
// 2012: Pierre Guillot, Eliott Paris & Julien Colafrancesco.
// 2012-2015: Pierre Guillot & Eliott Paris.
// 2015: Pierre Guillot & Eliott Paris & Thomas Le Meur (Light version)
// 2016-2017: Pierre Guillot.
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#pragma once

#include "Hoa_HrirIrc1002C2D.hpp"
#include "Hoa_HrirIrc1002C3D.hpp"

namespace hoa
{
    // ================================================================================ //
    // HRIR //
    // ================================================================================ //
    
    //! @brief Manage an impulse responses to decode in binaural mode.
    //! @details Manage an impulse responses to decode in binaural mode
    //! in 2d / 3d / simple and double precision.
    template <Dimension D, typename T>
    class Hrir
    {
    public:
        
        //! @brief Get the impulse response of the HRTFs
        //! @details  The impulse response may be 2d/3d with simple /double precision.
        virtual const float* getImpulse() noexcept = 0;
    };

    template<>
    class Hrir <Hoa2d, float>
    {
    public:
        
        //! @brief Gets the order of decomposition used to compute the matrices.
        //! @return The order of decomposition.
        static size_t getOrderOfDecomposition() noexcept { return 5ul; }
        
        //! @brief Gets the number rows of the matrices.
        //! @details Or the size of the responses used to compute the matrices.
        //! @return The number rows of the matrices.
        static size_t getNumberOfRows() noexcept
        {
            return 512ul;
        }
        
        //! @brief Gets the number columns of the matrices.
        //! @details Or the number of harmonics used to compute the matrices.
        //! @return The number columns of the matrices.
        static size_t getNumberOfColumns() noexcept
        {
            return 11ul;
        }
        
        //! @brief Gets the size of the matrices (rows * columns).
        //! @return The size of the matrices.
        static size_t getMatricesSize() noexcept
        {
            return 5632ul;
        }
        
        //! @brief Get the HRIR matrix for the left ear.
        //! @return The HRIR matrix for the left ear.
        static const float* getLeftMatrix() noexcept
        {
            return Irc1002C_float_2d_left;
        }
        
        //! @brief Get the HRIR matrix for the right ear.
        //! @return The HRIR matrix for the right ear.
        static const float* getRightMatrix() noexcept
        {
            return Irc1002C_float_2d_right;
        }
    };
    
    // ================================================================================ //
    // HRIR 2D double //
    // ================================================================================ //

    template<>
    class Hrir <Hoa2d, double>
    {
    public:
        
        //! @brief Gets the order of decomposition used to compute the matrices.
        //! @return The order of decomposition.
        static size_t getOrderOfDecomposition() noexcept
        {
            return 5ul;
        }
        
        //! @brief Gets the number rows of the matrices.
        //! @details Or the size of the responses used to compute the matrices.
        //! @return The number rows of the matrices.
        static size_t getNumberOfRows() noexcept
        {
            return 512ul;
        }
        
        //! @brief Gets the number columns of the matrices.
        //! @details Or the number of harmonics used to compute the matrices.
        //! @return The number columns of the matrices.
        static size_t getNumberOfColumns() noexcept
        {
            return 11ul;
        }
        
        //! @brief Gets the size of the matrices (rows * columns).
        //! @return The nsize of the matrices.
        static size_t getMatricesSize() noexcept
        {
            return 5632ul;
        }
        
        //! @brief Get the HRIR matrix for the left ear.
        //! @return The HRIR matrix for the left ear.
        static const double* getLeftMatrix() noexcept
        {
            return Irc1002C_double_2d_left;
        }
        
        //! @brief Get the HRIR matrix for the right ear.
        //! @return The HRIR matrix for the right ear.
        static const double* getRightMatrix() noexcept
        {
            return Irc1002C_double_2d_right;
        }
    };
    
    // ================================================================================ //
    // HRIR 3D float //
    // ================================================================================ //

    template<> class Hrir <Hoa3d, float>
    {
    public:
        
        //! @brief Gets the order of decomposition used to compute the matrices.
        //! @return The order of decomposition.
        static size_t getOrderOfDecomposition() noexcept
        {
            return 3ul;
        }
        
        //! @brief Gets the number rows of the matrices.
        //! @details Or the size of the responses used to compute the matrices.
        //! @return The number rows of the matrices.
        static size_t getNumberOfRows() noexcept
        {
            return 512ul;
        }
        
        //! @brief Gets the number columns of the matrices.
        //! @details Or the number of harmonics used to compute the matrices.
        //! @return The number columns of the matrices.
        static size_t getNumberOfColumns() noexcept
        {
            return 16ul;
        }
        
        //! @brief Gets the size of the matrices (rows * columns).
        //! @return The nsize of the matrices.
        static size_t getMatricesSize() noexcept
        {
            return 8192ul;
        }
        
        //! @brief Get the HRIR matrix for the left ear.
        //! @return The HRIR matrix for the left ear.
        static const float* getLeftMatrix() noexcept
        {
            return Irc1002C_float_3d_left;
        }
        
        //! @brief Get the HRIR matrix for the right ear.
        //! @return The HRIR matrix for the right ear.
        static const float* getRightMatrix() noexcept
        {
            return Irc1002C_float_3d_right;
        }
    };

    // ================================================================================ //
    // HRIR 3D double //
    // ================================================================================ //
    
    template<>
    class Hrir <Hoa3d, double>
    {
    public:
        
        //! @brief Gets the order of decomposition used to compute the matrices.
        //! @return The order of decomposition.
        static size_t getOrderOfDecomposition() noexcept
        {
            return 3ul;
        }
        
        //! @brief Gets the number rows of the matrices.
        //! @details Or the size of the responses used to compute the matrices.
        //! @return The number rows of the matrices.
        static size_t getNumberOfRows() noexcept
        {
            return 512ul;
        }
        
        //! @brief Gets the number columns of the matrices.
        //! @details Or the number of harmonics used to compute the matrices.
        //! @return The number columns of the matrices.
        static size_t getNumberOfColumns() noexcept
        {
            return 16ul;
        }
        
        //! @brief Gets the size of the matrices (rows * columns).
        //! @return The nsize of the matrices.
        static size_t getMatricesSize() noexcept
        {
            return 8192ul;
        }
        
        //! @brief Get the HRIR matrix for the left ear.
        //! @return The HRIR matrix for the left ear.
        static const double* getLeftMatrix() noexcept
        {
            return Irc1002C_double_3d_left;
        }
        
        //! @brief Get the HRIR matrix for the right ear.
        //! @return The HRIR matrix for the right ear.
        static const double* getRightMatrix() noexcept
        {
            return Irc1002C_double_3d_right;
        }
    };

}
