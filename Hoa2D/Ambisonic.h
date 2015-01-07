/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_AMBISONIC
#define DEF_HOA_2D_AMBISONIC

#include "../Hoa.h"

namespace Hoa2D
{
    //! The ambisonic class.
    /**
     Most of the ambisonic classes inherit from this classe. It computes the number of harmonics, their degrees and their orders depending of the decomposition order. etc...
     */
    class Ambisonic
    {
    protected:
        unsigned int	m_order;
        unsigned int	m_number_of_harmonics;
        long*           m_harmonics_orders;
    
public:
        //! The ambisonic constructor.
        /** The ambisonic constructor allocates and initializes the generale member values depending of a decomposition order.
         @param     order	The order, must be at least 1.
         */
        Ambisonic(unsigned int order);
        
        //! The ambisonic destructor.
        /** The ambisonic destructor.
         */
        ~Ambisonic();
        
        //! Retrieve the decomposition order.
        /** Retrieve the decomposition order.
            
            @return The order.
         */
        unsigned int getDecompositionOrder() const
        {
            return m_order;
        }
        
        //! Retrieve the number of harmonics.
        /** Retrieve the number of harmonics.
         
            @return The number of harmonics.
         */
        unsigned int getNumberOfHarmonics() const
        {
            return m_number_of_harmonics;
        }
        
        //! Retrieve the order of an harmonic.
        /** The order of an harmonic is in the range -order to order. The harmonics are sorted by their orders, from 0 to the decomposition order and, in each order, there are the 2 harmonics with the orders -order and order. For the first orders, the harmonics arrangement is h[0] h[-1] h[1] h[-2] h[2] h[-3] h[3]etc. with h[order].
         
            @param     index	The index of an harmonic.
            @return    The method returns the order of the harmonic if the harmonic exists, otherwise the function generates an error.
            @see       getHarmonicDegree()
            @see       getHarmonicName()
         */
        long getHarmonicOrder(unsigned int index) const
        {
            assert(index < m_number_of_harmonics);
            return m_harmonics_orders[index];
        }
        
        //! Retrieve the order of an harmonic.
        /** The orders of the harmonics are in the range 0 to the decomposition order. Each order contains 2 harmonics with the orders -order and order. For the first orders, the harmonics arrangement is h[0] h[-1] h[1] h[-2] h[2] h[-3] h[3], etc. with h[order].
         
            @param     index	The index of an harmonic.
            @return    The method returns the order of the harmonic if the harmonic exists, otherwise the function generates an error.
            @see       getHarmonicOrder()
            @see       getHarmonicName()
         */
        long getHarmonicDegree(unsigned int index) const
        {
            assert(index < m_number_of_harmonics);
            return abs(m_harmonics_orders[index]);
        }
        
        //! Retrieve the index of an harmonic.
        /** The orders of the harmonics are in the range 0 to the decomposition order. Each order contains 2 harmonics with the orders -order and order. For the first orders, the harmonics arrangement is h[0] h[-1] h[1] h[-2] h[2] h[-3] h[3], etc. with h[order].
         
         @param     order	The order an harmonic.
         @return    The method returns the index of the harmonic if the harmonic exists, otherwise the function generates an error.
         @see       getHarmonicOrder()
         @see       getHarmonicName()
         */
        inline unsigned int getHarmonicIndex(const int harmOrder) const
        {
            assert(abs(harmOrder) <= getDecompositionOrder());
            if(harmOrder < 0)
                return -harmOrder * 2 - 1;
            else
                return harmOrder * 2;
        };
        
        //! Retrieve a name for an harmonic.
        /** Retrieve a name for an harmonic in a std::string format that will be "harmonic order".
         
            @param     index	The index of an harmonic.
            @return    The method returns a name for the harmonic that contains its order if the harmonic exists, otherwise the function generates an error.
         
            @see       getHarmonicDegree()
            @see       getHarmonicOrder()
         */
        std::string getHarmonicName(unsigned int index) const
        {
            assert(index < m_number_of_harmonics);
            return "Harmonic " + int_to_string(getHarmonicOrder(index));
        }
    };
}


#endif


