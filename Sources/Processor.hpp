/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_PROCESSOR_LIGHT
#define DEF_HOA_PROCESSOR_LIGHT

#include "Harmonics.hpp"
#include "Planewaves.hpp"

namespace hoa
{
    //! The processor.
    /** The processor owns a set of harmonics or planewaves and performs digital signal processing.
     */
    template <typename T> class Processor
    {
    public:
        
        //! This method performs the processing.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The input array and the outputs array depends of the template and the processing.
         @param     input  The input array.
         @param     outputs The outputs array.
         */
        virtual void process(const T* input, T* outputs) const = 0;
    };
    
    //! The harmonic processor.
    /** The harmonic processor owns a set of harmonics depending on the order of decomposition.
     */
    template <Dimension D, typename T> class Processor< Harmonic<D, T> >
    {
    private:
        
        const ulong                 m_order_of_decomposition;
        const ulong                 m_number_of_harmonics;
        vector< Harmonic<D, T> >    m_harmonics;
    public:
        
        //! The harmonic processor constructor.
        /** The harmonic processor constructor allocates and initializes the generale member values depending on a order of decomposition .
         @param order    The order of decomposition, must be at least 1.
         */
        Processor(const ulong order) noexcept :
        m_order_of_decomposition(order),
        m_number_of_harmonics(Harmonic<D, T>::getNumberOfHarmonics(order))
        {
            for(ulong i = 0; i < m_number_of_harmonics; i++)
            {
                m_harmonics.push_back(Harmonic<D, T>(i));
            }
        }
        
        //! The ambisonic destructor.
        /** The ambisonic destructor.
         */
        ~Processor() noexcept
        {
            m_harmonics.clear();
        }
        
        //! Retrieve the order of decomposition.
        /** Retrieve the order of decomposition.
         @return The order.
         */
        inline ulong getDecompositionOrder() const noexcept
        {
            return m_order_of_decomposition;
        }
        
        //! Retrieve the number of harmonics.
        /** Retrieve the number of harmonics.
         @return The number of harmonics.
         */
        inline ulong getNumberOfHarmonics() const noexcept
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
        inline long getHarmonicOrder(const ulong index) const noexcept
        {
            return m_harmonics[index].getOrder();
        }
        
        //! Retrieve the degree of an harmonic.
        /** The degrees of the harmonics are in the range 0 to order of decomposition.
         @param     index	The index of an harmonic.
         @return    The method returns the degree of the harmonic.
         @see       getHarmonicOrder()
         @see       getHarmonicName()
         */
        inline ulong getHarmonicDegree(const ulong index) const noexcept
        {
            return m_harmonics[index].getDegree();
        }
        
        //! Retrieve the index of an harmonic.
        /** The indices of the harmonics are in the range 0 to number of harmonics - 1.
         @param     order	The order an harmonic.
         @return    The method returns the index of the harmonic.
         @see       getHarmonicOrder()
         @see       getHarmonicName()
         */
        inline ulong getHarmonicIndex(const ulong degree, long order) const noexcept
        {
            return Harmonic<D, T>::getHarmonicIndex(order, order);
        }
        
        //! Retrieve the name of an harmonic.
        /** Retrieve the name of an harmonic.
         @param     index	The index of an harmonic.
         @return    The method returns the name of the harmonic that contains its degree and its order.
         @see       getHarmonicDegree()
         @see       getHarmonicOrder()
         */
        string getHarmonicName(const ulong index) const noexcept
        {
            return m_harmonics[index].getName();
        }
    };
}

#endif


