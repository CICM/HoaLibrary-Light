/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_AMBISONIC
#define DEF_HOA_2D_AMBISONIC

#include "../Hoa.hpp"

namespace hoa
{
    //! The 2D ambisonic class.
    /**
     Most of the 2D ambisonic classes inherit from this classe. It computes the number of harmonics, their degrees and their orders depending of the decomposition order. The harmonics are sorted by their orders, from 0 to the decomposition order and, in each order, there are the 2 harmonics with the orders -order and order. For the first orders, the harmonics arrangement is h[0] h[-1] h[1] h[-2] h[2] h[-3] h[3]etc. with h[order].
     */
    template <typename T> class Ambisonic2D : public hoa::Ambisonic<T>
    {
    public:
        //! The ambisonic constructor.
        /** The ambisonic constructor allocates and initializes the generale member values depending of a decomposition order.
         @param     order	The order, must be at least 1.
         */
        Ambisonic2D(unsigned long order) noexcept : hoa::Ambisonic<T>(order, order * 2 + 1)
        {
            ;
        }
        
        //! The ambisonic destructor.
        /** The ambisonic destructor.
         */
        virtual ~Ambisonic2D()
        {
            ;
        }
        
        //! Retrieve the order of an harmonic.
        /** The orders of the harmonic are in the range -order to order.
         @param     index	The index of an harmonic.
         @return    The method returns the order of the harmonic.
         @see       getHarmonicDegree()
         @see       getHarmonicName()
         */
        inline long getHarmonicOrder(const unsigned long index) const noexcept override
        {
            if(index % 2)
            {
                return -0.5 * (long)(index + 1);
            }
            else
            {
                return 0.5 * index;
            }
        }
        
        //! Retrieve the degree of an harmonic.
        /** The degrees of the harmonics are in the range 0 to decomposition order.
         @param     index	The index of an harmonic.
         @return    The method returns the order of the harmonic.
         @see       getHarmonicOrder()
         @see       getHarmonicName()
         */
        inline unsigned long getHarmonicDegree(unsigned long index) const noexcept override
        {
            return (index + index % 2) * 0.5;
        }
        
        //! Retrieve the index of an harmonic.
        /** The indices of the harmonics are in the range 0 to number of harmonics - 1.
         @param     order	The order an harmonic.
         @return    The method returns the index of the harmonic.
         @see       getHarmonicOrder()
         @see       getHarmonicName()
         */
        inline unsigned long getHarmonicIndex(long order) const noexcept
        {
            if(order < 0)
            {
                return -order * 2 - 1;
            }
            else
            {
                return order * 2;
            }
        };
        
        //! Retrieve the index of an harmonic.
        /** The indices of the harmonics are in the range 0 to number of harmonics - 1.
         @param     order	The order an harmonic.
         @return    The method returns the index of the harmonic.
         @see       getHarmonicOrder()
         @see       getHarmonicName()
         */
        inline unsigned long getHarmonicIndex(unsigned long degree, long order) const noexcept override
        {
            return getHarmonicIndex(order);
        }
    };
}


#endif


