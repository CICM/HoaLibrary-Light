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

namespace hoa
{
    //! Math utilities and funtions.
    template<typename T>
    class math
    {
    public:
        
        math() = delete;
        ~math() = delete;
        
        //! @brief Returns π constant
        static constexpr T pi() { return 3.14159265358979323846264338327950288; }
        
        //! @brief Returns 2π
        static constexpr T two_pi() { return pi() * 2.; }
        
        //! @brief Returns π/2
        static constexpr T pi_over_two() { return pi() * 0.5; }
        
        //! @brief Returns π/4
        static constexpr T pi_over_four() { return pi() * 0.25; }
        
        //! @brief Wraps the value between 0 and π
        static T wrap_pi(T value)
        {
            while(value < -pi()) { value += two_pi(); }
            while(value >= pi()) { value -= two_pi(); }
            return value;
        }
        
        //! @brief Wraps the value between 0 and 2π
        static T wrap_two_pi(T value)
        {
            while(value < 0.) { value += two_pi(); }
            while(value >= two_pi()) { value -= two_pi(); }
            return value;
        }
        
    };
}
