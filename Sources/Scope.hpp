/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_SCOPE_LIGHT
#define DEF_HOA_SCOPE_LIGHT

#include "Encoder.hpp"
#include "Planewaves.hpp"

namespace hoa
{
    
    template <Dimension D, typename T> class Scope;
    
    //! The ambisonic scope.
    /** The scope discretize a circle by a set of point and uses a decoder to project the circular harmonics on it. This class should be used for graphical interfaces outside the digital signal processing if the number of points to discretize the circle is very large. Then you should prefer to record snapshot of the circular harmonics and to call the process method at an interval adapted to a graphical rendering.
     */
    template <typename T> class Scope<Hoa2d, T> : public Encoder<Hoa2d, T>, protected Planewave<Hoa2d, T>::Processor
    {
    private:
        T*  m_matrix;
        T*  m_vector;
        T   m_maximum;
    public:
        
        //! The scope constructor.
        /**	The scope constructor allocates and initialize the member values to computes circular harmonics projection on a circle depending on a decomposition order and a circle discretization. The circle is discretized by the number of points. The order must be at least 1. The number of points and column should be at least 3 (but it's very low).
         
         @param     order            The order.
         @param     numberOfPoints   The number of points.
         */
        Scope(ulong order, ulong numberOfPoints) :
        Encoder<Hoa2d, T>(order),
        Planewave<Hoa2d, T>::Processor(numberOfPoints)
        {
            m_matrix = new T[Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves() * Encoder<Hoa2d, T>::getNumberOfHarmonics()];
            m_vector = new T[Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves()];
            
            const T factor = 1. / (T)(Encoder<Hoa2d, T>::getDecompositionOrder() + 1.);
            for(ulong i = 0; i < Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves(); i++)
            {
                Encoder<Hoa2d, T>::setAzimuth(Planewave<Hoa2d, T>::Processor::getPlanewaveAzimuth(i));
                Encoder<Hoa2d, T>::process(&factor, m_matrix + i * Encoder<Hoa2d, T>::getNumberOfHarmonics());
                m_matrix[i * Encoder<Hoa2d, T>::getNumberOfHarmonics()] = factor * 0.5;
            }
            for(ulong i = 0; i < Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves(); i++)
            {
                m_vector[i] = 0.;
            }
            m_maximum = 0;
        }
        //! The Scope destructor.
        /**	The Scope destructor free the memory.
         */
        ~Scope()
        {
            delete [] m_matrix;
            delete [] m_vector;
        }
        
        //! Retrieve the number of points.
        /**	Retrieve the number of points used to discretize the ambisonic circle.
         @return     This method returns the number of points used to discretize the circle.
         */
        inline ulong getNumberOfPoints() const noexcept
        {
            return Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves();
        }
        
        //! Retrieve the value of a point of the circular harmonics projection.
        /**	Retrieve the result value of the circular harmonics projection for a given point defined by an index. The absolute of the value can be used as the radius of the point for a 2 dimentionnal representation. For the index, 0 is the 0 azimtuh of the circle. The maximum index must be the number of points - 1.
         @param     index   The point index of the point.
         @return    This method returns the value of a point of the ambisonic circle.
         */
        inline T getPointValue(const ulong index) const noexcept
        {
            return m_vector[index];
        }
        
        //! Retrieve the radius of a point of the circular harmonics projection.
        /**	Retrieve the radius of the circular harmonics projection for a given point defined by an index. This the absolute of the result of the projection. For the index, 0 is the 0 azimtuh of the circle. The maximum index must be the number of points - 1.
         @param     pointIndex   The point index of the point.
         @return    This method returns the radius of a point of the ambisonic circle.
         */
        inline T getPointRadius(const ulong index) const noexcept
        {
            return fabs(m_vector[index]);
        }
        
        //! Retrieve the azimuth of a point of the circular harmonics projection.
        /**	Retrieve the azimuth of the circular harmonics projection for a given point defined by an index.The maximum index must be the number of points - 1.
         @param     pointIndex   The point index of the point.
         @return    This method returns the azimuth of a point of the ambisonic circle.
         */
        inline T getPointAzimuth(const ulong index) const noexcept
        {
            return Planewave<Hoa2d, T>::Processor::getPlanewaveAzimuth(index);
        }
        
        //! Retrieve the abscissa of a point of the circular harmonics projection.
        /**	Retrieve the abscissa of the circular harmonics projection for a given point defined by an index.The maximum index must be the number of points - 1.
         
         @param     pointIndex   The point index of the point.
         @return    This method returns the abscissa of a point of the ambisonic circle.
         
         @see       getOrdinate
         */
        inline double getPointAbscissa(const ulong index) const noexcept
        {
            return fabs(m_vector[index]) * Planewave<Hoa2d, T>::Processor::getPlanewaveAbscissa(index);
        }
        
        //! Retrieve the ordinate of a point of the circular harmonics projection.
        /**	Retrieve the ordinate of the circular harmonics projection for a given point defined by an index.The maximum index must be the number of points - 1.
         
         @param     pointIndex   The point index of the point.
         @return    This method returns the ordinate of a point of the ambisonic circle.
         
         @see       getAbscissa
         */
        inline double getPointOrdinate(const ulong index) const noexcept
        {
            return fabs(m_vector[index]) * Planewave<Hoa2d, T>::Processor::getPlanewaveOrdinate(index);
        }
        
        //! This method performs the circular harmonics projection.
        /**	You should use this method to compute the projection of the circular harmonics over an ambisonics circle. The inputs array contains the circular harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         */
        inline void process(const T* inputs) noexcept
        {
            Signal<T>::matrix_vector_mul(Encoder<Hoa2d, T>::getNumberOfHarmonics(), Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves(), inputs, m_matrix, m_vector);
            m_maximum = fabs(Signal<T>::vector_max(Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves(), m_vector));
            if(m_maximum > 1.)
            {
                Signal<T>::vector_scale(Planewave<Hoa2d, T>::Processor::getNumberOfPlanewaves(), (1. / m_maximum), m_vector);
            }
        }
    };
}

#endif



