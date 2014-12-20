/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_3D_ENCODER__
#define __DEF_HOA_3D_ENCODER__

#include "Ambisonic_3D.h"

namespace Hoa3D
{
    //! The ambisonic encoder.
    /** The encoder should be used to encode a signal in the spherical harmonics domain depending of an order of decomposition. It allows to control the azimuth and the elevation of the encoding. If you want to spatialize with distance compensation, you should use the Map class.
        
        @see Map
     */
    class Encoder : public Ambisonic
    {
        
    private:
        long            m_elevation;
        long            m_azimuth;
        double**        m_azimuth_matrix;
        double**        m_elevation_matrix;
        double*         m_normalization;
        
    public:
        
        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes spherical harmonics coefficients depending of a decomposition order. The order must be at least 1.
         
            @param     order	The order.
         */
        Encoder(unsigned int order);
        
        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
        ~Encoder();
        
        //! This method set the angle of azimuth.
        /**	The angle of azimuth in radian and you should prefer to use it between 0 and 2 Pi to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is Pi/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield.
         
            @param     azimuth	The azimuth.
            @see       setElevation()
         */
        void setAzimuth(const double azimuth);
        
        //! This method set the angle of elevation.
        /**	The angle of elevation in radian and you should prefer to use it between 0 and 2 Pi to avoid recursive wrapping of the value. The direction of rotation is from bottom to the top. The 0 radian is centered at the "front" of the soundfield, then Pi/2 is at the top, -Pi/2 is at the bottom and Pi is behind. Note that if the angle of elevation is between Pi/2 and 3*Pi/2, the azimuth is reversed.
         
            @param     elevation The elevation.
            @see       setElevation()
         */
        void setElevation(const double elevation);
        
        /**	Retreive the normalization of an harmonics
         
            @param     index The index of the harmonics.
         */
        double getNormalization(const unsigned int index) const
        {
            assert(index < m_number_of_harmonics);
            return m_normalization[index];
        };
        
        /**	Retreive the azimuth of a source
         */
        double getAzimuth() const
        {
            return (double)m_azimuth * HOA_2PI / (double)(NUMBEROFCIRCLEPOINTS - 1);
        };
        
        /**	Retreive the azimuth of a source
         */
        double getElevation() const
        {
            return (double)m_elevation * HOA_2PI / (double)(NUMBEROFCIRCLEPOINTS - 1);
        };
        
        //! This method performs the encoding with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     input    The input sample.
            @param     outputs  The outputs array.
         */
        void process(const float input, float* outputs);
        
        //! This method performs the encoding with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     input	The input sample.
            @param     outputs  The outputs array.
         */
        void process(const double input, double* outputs);
    };
}

#endif



