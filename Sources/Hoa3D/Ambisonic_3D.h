/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_3D_AMBISONIC__
#define __DEF_HOA_3D_AMBISONIC__

#include "../Hoa.h"

namespace Hoa3D
{
    //! The ambisonic class.
    /** The ambisonics classes inherit from this classe. It computes the number of harmonics depending of the decomposition order. The harmonics are sorted by their degrees, from 0 to the decomposition order. Each degree contains 2 * degree + 1 harmonics with orders in the range -degree to degree. The harmonic input and output arrays in process method of ambisonic classes must have this configuration. For the first degrees, the harmonics arrangement is (0,0) (1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) with (degree,order).
     */
    class Ambisonic : public Hoa::Ambisonic
    {      
    public:
        
        //! The ambisonic constructor.
        /** The ambisonic constructor allocates and initializes the generale member values depending of a decomposition order. The order must be at least 1.
            @param     order	The order.
         */
        Ambisonic(unsigned long order) noexcept  : Hoa::Ambisonic(order, (order + 1) * (order + 1))
        {
            ;
        }
        
        //! The ambisonic destructor.
        /**	The ambisonic destructor free the memory.
         */
        ~Ambisonic()
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
        inline long getHarmonicOrder(unsigned long index) const noexcept override
        {
            const unsigned long order = sqrt(index);
            return index - (order * (order + 1));
        };
        
        //! Retrieve the degree of an harmonic.
        /** The degrees of the harmonics are in the range 0 to decomposition order.
         @param     index	The index of an harmonic.
         @return    The method returns the order of the harmonic.
         @see       getHarmonicOrder()
         @see       getHarmonicName()
         */
        inline unsigned long getHarmonicDegree(unsigned long index) const noexcept override
        {
            return sqrt(index);
        };
        
        //! Retrieve the index of an harmonic.
        /** The indices of the harmonics are in the range 0 to number of harmonics - 1.
         @param     order	The order an harmonic.
         @return    The method returns the index of the harmonic.
         @see       getHarmonicOrder()
         @see       getHarmonicName()
         */
        inline unsigned long getHarmonicIndex(unsigned long degree, long order) const noexcept
        {
            return degree * (degree + 1) + order;
        };
    };
}

#endif


