/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_SCOPE
#define DEF_HOA_2D_SCOPE

#include "Ambisonic.h"
#include "Decoder.h"

namespace Hoa2D
{
    //! The ambisonic scope.
    /** The scope discretize a circle by a set of point and uses a decoder to project the circular harmonics on it. This class should be used for graphical interfaces outside the digital signal processing if the number of points to discretize the circle is very large. Then you should prefer to record snapshot of the circular harmonics and to call the process method at an interval adapted to a graphical rendering.
     */
    class Scope : public Ambisonic
    {
    private:
        unsigned int    m_number_of_points;
        double*         m_harmonics;
        double*         m_matrix;
        DecoderRegular* m_decoder;
    public:
        
        //! The scope constructor.
        /**	The scope constructor allocates and initialize the member values to computes circular harmonics projection on a circle depending on a decomposition order and a circle discretization. The circle is discretized by the number of points. The order must be at least 1. The number of points and column should be at least 3 (but it's very low).
         
            @param     order            The order.
            @param     numberOfPoints   The number of points.
         */
        Scope(unsigned int order, unsigned int numberOfPoints);
        
        //! The Scope destructor.
        /**	The Scope destructor free the memory.
         */
        ~Scope();
        
        //! Retrieve the number of points.
        /**	Retrieve the number of points used to discretize the ambisonic circle.
         
            @return     This method returns the number of points used to discretize the circle.
         */
        inline unsigned int getNumberOfPoints() const
        {
            return m_number_of_points;
        }
        
        //! Retrieve the value of a point of the circular harmonics projection.
        /**	Retrieve the result value of the circular harmonics projection for a given point defined by an index. The absolute of the value can be used as the radius of the point for a 2 dimentionnal representation. For the index, 0 is the 0 azimtuh of the circle. The maximum index must be the number of points - 1.
         
         @param     pointIndex   The point index of the point.
         @return    This method returns the value of a point of the ambisonic circle.
         
         @see       getradius
         @see       getAzimuth
         */
        inline double getValue(unsigned int pointIndex) const
        {
            assert(pointIndex < m_number_of_points);
            return m_matrix[pointIndex];
        }
        
        //! Retrieve the radius of a point of the circular harmonics projection.
        /**	Retrieve the radius of the circular harmonics projection for a given point defined by an index. This the absolute of the result of the projection. For the index, 0 is the 0 azimtuh of the circle. The maximum index must be the number of points - 1.
         
            @param     pointIndex   The point index of the point.
            @return    This method returns the radius of a point of the ambisonic circle.
            
            @see       getAzimuth
            @see       getValue
         */
        inline double getRadius(unsigned int pointIndex) const
        {
            assert(pointIndex < m_number_of_points);
            return fabs(m_matrix[pointIndex]);
        }
		
        //! Retrieve the azimuth of a point of the circular harmonics projection.
        /**	Retrieve the azimuth of the circular harmonics projection for a given point defined by an index.The maximum index must be the number of points - 1.
         
            @param     pointIndex   The point index of the point.
            @return    This method returns the azimuth of a point of the ambisonic circle.
         
            @see       getValue
            @see       getRadius
         */
        inline double getAzimuth(unsigned int pointIndex) const
        {
            assert(pointIndex < m_number_of_points);
            return (double)pointIndex * HOA_2PI / (double)m_number_of_points;
        }
        
        //! Retrieve the abscissa of a point of the circular harmonics projection.
        /**	Retrieve the abscissa of the circular harmonics projection for a given point defined by an index.The maximum index must be the number of points - 1.
         
            @param     pointIndex   The point index of the point.
            @return    This method returns the abscissa of a point of the ambisonic circle.

            @see       getOrdinate
         */
        inline double getAbscissa(unsigned int pointIndex) const
        {
            assert(pointIndex < m_number_of_points);
            return abscissa(fabs(m_matrix[pointIndex]), getAzimuth(pointIndex));
        }
		
        //! Retrieve the ordinate of a point of the circular harmonics projection.
        /**	Retrieve the ordinate of the circular harmonics projection for a given point defined by an index.The maximum index must be the number of points - 1.
         
            @param     pointIndex   The point index of the point.
            @return    This method returns the ordinate of a point of the ambisonic circle.
         
            @see       getAbscissa
         */
        inline double getOrdinate(unsigned int pointIndex) const
        {
            assert(pointIndex < m_number_of_points);
            return ordinate(fabs(m_matrix[pointIndex]), getAzimuth(pointIndex));
        }
        
        //! This method performs the circular harmonics projection with single precision.
        /**	You should use this method to compute the projection of the circular harmonics over an ambisonics circle. The inputs array contains the circular harmonics samples and the minimum size must be the number of harmonics.
         
            @param     inputs   The inputs array.
         */
        void process(const float* inputs);
        
        //! This method performs the circular harmonics projection with double precision.
        /**	You should use this method to compute the projection of the circular harmonics over an ambisonics circle. The inputs array contains the circular harmonics samples and the minimum size must be the number of harmonics.
         
            @param     inputs   The inputs array.
         */
        void process(const double* inputs);
    };
	
} // end of namespace Hoa3D

#endif


