/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_HARMONICS_LIGHT
#define DEF_HOA_HARMONICS_LIGHT

#include "Math.hpp"
#include "Signal.hpp"

namespace hoa
{
    //! The harmonic class owns basic ordering informations.
    /** The harmonic allows to retrieves informations about its ACN ordering, the degree and the order.
     */
    template <Dimension D, typename T> class Harmonic
    {
    public:
        
        //! The harmonic constructor.
        /** The harmonic constructor allocates and initializes the generale member values depending of an index.
         @param index    The index must be at least 1.
         */
        Harmonic(const ulong _index) noexcept = 0;
        
        //! The harmonic destructor.
        /** The harmonic destructor free the memory.
         */
        virtual ~Harmonic() noexcept = 0;
        
        //! Get the index of the harmonic.
        /** The method returns the index of the harmonic.
         @return     The index.
         */
        virtual ulong getIndex() const noexcept = 0;

        //! Get the degree of the harmonic.
        /** The method returns the degree of the harmonic.
         @return     The degree.
         */
        virtual ulong getDegree() const noexcept = 0;
        
        //! Get the order of the harmonic.
        /** The method returns the order of the harmonic.
         @return     The order.
         */
        virtual long getOrder() const noexcept = 0;
        
        //! Get the name of the harmonic.
        /** The method returns the name of the harmonic.
         @return     The name.
         */
        virtual string getName() const noexcept = 0;
    };
    
    template<typename T> class Harmonic<Hoa2d, T>
    {
    private:
        ulong m_index;
        ulong m_degree;
        long  m_order;
    public:
        
        //! The harmonic constructor.
        /** The harmonic constructor allocates and initializes the generale member values depending of an index.
         @param index    The index must be at least 1.
         */
        Harmonic(const ulong _index) noexcept :
        m_index(_index),
        m_degree((_index + _index % 2) * 0.5),
        m_order(m_degree * (1 - (_index % 2) * 2))
        {
            ;
        }
        
        //! The harmonic destructor.
        /** The harmonic destructor free the memory.
         */
        ~Harmonic() noexcept
        {
            ;
        }
        
        //! Get the index of the harmonic.
        /** The method returns the index of the harmonic.
         @return     The index.
         */
        inline ulong getIndex() const noexcept
        {
            return m_index;
        }
        
        //! Get the degree of the harmonic.
        /** The method returns the degree of the harmonic.
         @return     The degree.
         */
        inline ulong getDegree() const noexcept
        {
            return m_degree;
        }
        
        //! Get the order of the harmonic.
        /** The method returns the order of the harmonic.
         @return     The order.
         */
        inline long getOrder() const noexcept
        {
            return m_order;
        }
        
        //! Get the name of the harmonic.
        /** The method returns the name of the harmonic.
         @return     The name.
         */
        inline string getName() const noexcept
        {
            return "Harmonic " + to_string(getDegree()) + " " + to_string(getOrder());
        }
        
        //! The harmonics class.
        /**
         All the classes that perform on the harmonics inherit from this classe. It computes the number of harmonics, their degrees and their orders depending of the order of decomposition.
         */
        class Processor
        {
        private:
            const ulong                 m_order_of_decomposition;
            const ulong                 m_number_of_harmonics;
            vector<Harmonic<Hoa2d, T> > m_harmonics;
        public:
            
            //! The ambisonic constructor.
            /** The ambisonic constructor allocates and initializes the generale member values depending of a decomposition order.
             @param order    The order of decomposition, must be at least 1.
             */
            Processor(const ulong order) noexcept :
            m_order_of_decomposition(order),
            m_number_of_harmonics(2 *  order + 1)
            {
                for(ulong i = 0; i < m_number_of_harmonics; i++)
                {
                    m_harmonics.push_back(Harmonic<Hoa2d, T>(i));
                }
            }
            
            //! The ambisonic destructor.
            /** The ambisonic destructor.
             */
            ~Processor()
            {
                m_harmonics.clear();
            }
            
            //! Retrieve the decomposition order.
            /** Retrieve the decomposition order.
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
            /** The degrees of the harmonics are in the range 0 to decomposition order.
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
                return abs(order) *  2 - (order < 0);
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
    };
        
    template <typename T> class Harmonic<Hoa3d, T>
    {
    private:
        ulong m_index;
        ulong m_degree;
        long  m_order;
    public:
        
        //! The harmonic constructor.
        /** The harmonic constructor allocates and initializes the generale member values depending of an index.
         @param index    The index must be at least 1.
         */
        Harmonic(const ulong _index) noexcept :
        m_index(_index),
        m_degree(sqrt(_index)),
        m_order(_index - (m_degree * (m_degree + 1)))
        {
            
        }
        
        //! The harmonic destructor.
        /** The harmonic destructor free the memory.
         */
        ~Harmonic() noexcept
        {
            ;
        }
        
        //! Get the index of the harmonic.
        /** The method returns the index of the harmonic.
         @return     The index.
         */
        inline ulong getIndex() const noexcept
        {
            return m_index;
        }
        
        //! Get the degree of the harmonic.
        /** The method returns the degree of the harmonic.
         @return     The degree.
         */
        inline ulong getDegree() const noexcept
        {
            return m_degree;
        }
        
        //! Get the order of the harmonic.
        /** The method returns the order of the harmonic.
         @return     The order.
         */
        inline long getOrder() const noexcept
        {
            return m_order;
        }
        
        //! Get the name of the harmonic.
        /** The method returns the name of the harmonic.
         @return     The name.
         */
        inline string getName() const noexcept
        {
            return "Harmonic " + to_string(getDegree()) + " " + to_string(getOrder());
        }
        
        //! The harmonics class.
        /**
         All the classes that perform on the harmonics inherit from this classe. It computes the number of harmonics, their degrees and their orders depending of the order of decomposition.
         */
        class Processor
        {
        protected:
            const ulong                 m_order_of_decomposition;
            const ulong                 m_number_of_harmonics;
            vector<Harmonic<Hoa3d, T> > m_harmonics;
        public:
            
            //! The ambisonic constructor.
            /** The ambisonic constructor allocates and initializes the generale member values depending of a decomposition order.
             @param order    The order of decomposition, must be at least 1.
             */
            Processor(const ulong order) noexcept :
            m_order_of_decomposition(order),
            m_number_of_harmonics((order + 1) * (order + 1))
            {
                for(ulong i = 0; i < m_number_of_harmonics; i++)
                {
                    m_harmonics.push_back(Harmonic<Hoa3d, T>(i));
                }
            }
            
            //! The ambisonic destructor.
            /** The ambisonic destructor.
             */
            ~Processor()
            {
                m_harmonics.clear();
            }
            
            //! Retrieve the decomposition order.
            /** Retrieve the decomposition order.
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
            /** The degrees of the harmonics are in the range 0 to decomposition order.
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
                return degree * (degree + 1) + order;
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
    };    
}

#endif


