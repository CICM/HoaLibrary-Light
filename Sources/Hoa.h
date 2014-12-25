/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_LIBRARY__
#define __DEF_HOA_LIBRARY__

namespace Hoa{};

#include "HoaDefs.h"
#include "HoaMath.h"
#include "HoaUtils.h"

using namespace Hoa;
using namespace std;

namespace Hoa
{
    //! The ambisonic class.
    /**
     Most of the ambisonic classes inherit from this classe. It computes the number of harmonics, their degrees and their orders depending of the decomposition order.
     */
    class Ambisonic
    {
    protected:
        const unsigned long	m_order_of_decomposition;
        const unsigned long	m_number_of_harmonics;
        
    public:
        //! The ambisonic constructor.
        /** The ambisonic constructor allocates and initializes the generale member values depending of a decomposition order.
         @param orderOfDecomposition    The order of decomposition, must be at least 1.
         @param numberOfHarmonics       The number of harmonics.
         */
        Ambisonic(unsigned long orderOfDecomposition, unsigned long numberOfHarmonics) noexcept :
        m_order_of_decomposition(orderOfDecomposition),
        m_number_of_harmonics(numberOfHarmonics)
        {
            ;
        }
        
        //! The ambisonic destructor.
        /** The ambisonic destructor.
         */
        virtual ~Ambisonic()
        {
            ;
        }
        
        //! Retrieve the decomposition order.
        /** Retrieve the decomposition order.
         @return The order.
         */
        inline unsigned long getDecompositionOrder() const noexcept
        {
            return m_order_of_decomposition;
        }
        
        //! Retrieve the number of harmonics.
        /** Retrieve the number of harmonics.
         @return The number of harmonics.
         */
        inline unsigned long getNumberOfHarmonics() const noexcept
        {
            return m_number_of_harmonics;
        }
        
        //! Retrieve the order of an harmonic.
        /** The orders of the harmonic are in the range -order to order.
         @param     index	The index of an harmonic.
         @return    The method returns the order of the harmonic.
         @see       getHarmonicDegree()
         @see       getHarmonicName()
         */
        virtual inline long getHarmonicOrder(unsigned long index) const noexcept = 0;
        
        //! Retrieve the degree of an harmonic.
        /** The degrees of the harmonics are in the range 0 to decomposition order.
         @param     index	The index of an harmonic.
         @return    The method returns the order of the harmonic.
         @see       getHarmonicOrder()
         @see       getHarmonicName()
         */
        virtual inline unsigned long getHarmonicDegree(unsigned long index) const noexcept = 0;
        
        //! Retrieve the index of an harmonic.
        /** The indices of the harmonics are in the range 0 to number of harmonics - 1.
         @param     order	The order an harmonic.
         @return    The method returns the index of the harmonic.
         @see       getHarmonicOrder()
         @see       getHarmonicName()
         */
        virtual inline unsigned long getHarmonicIndex(unsigned long degree, long order) const noexcept = 0;
        
        //! Retrieve a name for an harmonic.
        /** Retrieve a name for an harmonic in a string format.
         @param     index	The index of an harmonic.
         @return    The method returns a name for the harmonic that contains its degree and its order.
         @see       getHarmonicDegree()
         @see       getHarmonicOrder()
         */
        string getHarmonicName(unsigned long index) const noexcept
        {
            return "Harmonic " + num_to_string(getHarmonicDegree(index)) + " " + num_to_string(getHarmonicOrder(index));
        }
    };
}

#endif


