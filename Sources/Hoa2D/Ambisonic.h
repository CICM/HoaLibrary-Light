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
        const unsigned long	m_order;
        const unsigned long	m_number_of_harmonics;
    
    public:
        //! The ambisonic constructor.
        /** The ambisonic constructor allocates and initializes the generale member values depending of a decomposition order.
         @param     order	The order, must be at least 1.
         */
        Ambisonic(unsigned long order) noexcept :
        m_order(order),
        m_number_of_harmonics(order * 2 + 1)
        {
            ;
        }
        
        //! The ambisonic destructor.
        /** The ambisonic destructor.
         */
        ~Ambisonic()
        {
            ;
        }
        
        //! Retrieve the decomposition order.
        /** Retrieve the decomposition order.
            
            @return The order.
         */
        inline unsigned long getDecompositionOrder() const noexcept
        {
            return m_order;
        }
        
        //! Retrieve the number of harmonics.
        /** Retrieve the number of harmonics.
         
            @return The number of harmonics.
         */
        unsigned long getNumberOfHarmonics() const noexcept
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
        long getHarmonicOrder(unsigned long index) const noexcept
        {
            if(index % 2)
            {
                return index / 2;
            }
            else
            {
                return -(index + 1) / 2;
            }
        }
        
        //! Retrieve the degree of an harmonic.
        /** The orders of the harmonics are in the range 0 to the decomposition order. Each order contains 2 harmonics with the orders -order and order. For the first orders, the harmonics arrangement is h[0] h[-1] h[1] h[-2] h[2] h[-3] h[3], etc. with h[order].
         
            @param     index	The index of an harmonic.
            @return    The method returns the order of the harmonic if the harmonic exists, otherwise the function generates an error.
            @see       getHarmonicOrder()
            @see       getHarmonicName()
         */
        unsigned long getHarmonicDegree(unsigned long index) const noexcept
        {
            if(index % 2)
            {
                return index / 2;
            }
            else
            {
                return (index + 1) / 2;
            }
        }
        
        //! Retrieve the index of an harmonic.
        /** The orders of the harmonics are in the range 0 to the decomposition order. Each order contains 2 harmonics with the orders -order and order. For the first orders, the harmonics arrangement is h[0] h[-1] h[1] h[-2] h[2] h[-3] h[3], etc. with h[order].
         
         @param     order	The order an harmonic.
         @return    The method returns the index of the harmonic if the harmonic exists, otherwise the function generates an error.
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
        
        //! Retrieve a name for an harmonic.
        /** Retrieve a name for an harmonic in a std::string format that will be "harmonic order".
         
            @param     index	The index of an harmonic.
            @return    The method returns a name for the harmonic that contains its order if the harmonic exists, otherwise the function generates an error.
         
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


