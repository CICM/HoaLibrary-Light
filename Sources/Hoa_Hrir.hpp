// Copyright (c) 2012-2019 CICM - Universite Paris 8 - Labex Arts H2H.
// Authors :
// 2012: Pierre Guillot, Eliott Paris & Julien Colafrancesco.
// 2012-2015: Pierre Guillot & Eliott Paris.
// 2015: Pierre Guillot & Eliott Paris & Thomas Le Meur (Light version)
// 2016-2017: Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.

#pragma once

#include "Hoa_Hrir_Listen_1002C_2D.hpp"
#include "Hoa_Hrir_Listen_1002C_3D.hpp"
#include "Hoa_Hrir_Sadie_D2_2D.hpp"
#include "Hoa_Hrir_Sadie_D2_3D.hpp"

namespace hoa
{
    // ================================================================================ //
    // HRIR //
    // ================================================================================ //
    
    //! @brief Manage an impulse responses to decode in binaural mode.
    //! @details Manage an impulse responses to decode in binaural mode
    //! in 2d / 3d / simple and double precision.
    template<Dimension D, class HrirType>
    class Hrir
    {
    public:
        
        using hrir_t = HrirType;
        
        //! @brief Gets the order of decomposition used to compute the matrices.
        //! @return The order of decomposition.
        static constexpr size_t getOrderOfDecomposition() noexcept
        {
            return hrir_t::order;
        }
        
        //! @brief Gets the number rows of the matrices.
        //! @details Or the size of the responses used to compute the matrices.
        //! @return The number rows of the matrices.
        static constexpr size_t getNumberOfRows() noexcept
        {
            return hrir_t::responses_size;
        }
        
        //! @brief Gets the number columns of the matrices.
        //! @details Or the number of harmonics used to compute the matrices.
        //! @return The number columns of the matrices.
        static constexpr size_t getNumberOfColumns() noexcept
        {
            return hrir_t::number_of_harmonics;
        }
        
        //! @brief Gets the size of the matrices (rows * columns).
        //! @return The nsize of the matrices.
        static constexpr size_t getMatricesSize() noexcept
        {
            return getNumberOfRows() * getNumberOfColumns();
        }
        
        //! @brief Get the HRIR matrix for the left ear.
        //! @return The HRIR matrix for the left ear.
        template<typename FloatType>
        static FloatType const* getLeftMatrix()
        {
            return GetMatrix<FloatType>::left();
        }
        
        //! @brief Get the HRIR matrix for the right ear.
        //! @return The HRIR matrix for the right ear.
        template<typename FloatType>
        static FloatType const* getRightMatrix()
        {
            return GetMatrix<FloatType>::right();
        }
        
    private:
        
        template <typename FloatType, typename Unused = void>
        struct GetMatrix {
            static FloatType const* left() {}
            static FloatType const* right() {}
        };
        
        template <typename Unused>
        struct GetMatrix<float, Unused> {
            static float const* left() { return hrir_t::get_float_left(); }
            static float const* right() { return hrir_t::get_float_right(); }
        };
        
        template <typename Unused>
        struct GetMatrix<double, Unused> {
            static double const* left() { return hrir_t::get_double_left(); }
            static double const* right() { return hrir_t::get_double_right(); }
        };
    };
}
