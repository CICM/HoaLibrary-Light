/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_MAP
#define DEF_HOA_2D_MAP

#include "Ambisonic.h"
#include "Encoder.h"
#include "Wider.h"

namespace Hoa2D
{
    //! The ambisonic multi-encoder with distance compensation.
    /** The map is a multi Encoder with distance compensation. It uses intances of the Wider class to decrease the directionnality of sources by simulating fractionnal orders when the sources are inside the ambisonic circle and a simple diminution of the gain when the sources get away from the ambisonic circle.
     
        @see Encoder
     */
    class Map : public Ambisonic
    {
        
    private:
        
        unsigned int			m_number_of_sources;
        double*					m_gains;
		bool*					m_muted;
        int						m_first_source;
        
        double*                 m_azimuth;
        double*                 m_cosx;
        double*                 m_sinx;
        long*                   m_wide;
        double*                 m_wide_matrix;
        
    public:
        
        //! The map constructor.
        /**	The map constructor allocates and initialize the member values and classes depending of a decomposition order and the number of sources. The order and the number of sources must be at least 1.
         
            @param     order            The order.
            @param     numberOfSources	The number of sources.
         */
        Map(unsigned int order, unsigned int numberOfSources);
        
        //! The map destructor.
        /**	The map destructor free the memory and deallocate the member classes.
         */
        ~Map();
        
        //! This method retrieve the number of sources.
        /** Retrieve the number of sources.
         
            @return The number of sources.
         */
        unsigned int getNumberOfSources() const
        {
            return m_number_of_sources;
        };
        
        //! This method set the angle of azimuth of a source.
        /**	The angle of azimuth in radian and you should prefer to use it between 0 and 2 Pi to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is Pi/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield. The index must be between 0 and the number of sources - 1.
         
            @param     index	The index of the source.
            @param     azimuth	The azimuth.
            @see       setRadius()
         */
        void setAzimuth(const unsigned int index, const double azimuth);
        
        //! This method set the radius of a source.
        /**	The radius is between 0 and infinity. At 0, the source is in the center of the ambisonic circle and at 1, the source is at the limit of the ambisonic circle. Over 1, the source get away the ambisonic circle. The index must be between 0 and the number of sources - 1.
         
            @param     index	The index of the source.
            @param     radius   The radius.
            @see       setAzimuth()
         */
        void setRadius(const unsigned int index, const double radius);
		
		//! This method mute or unmute a source.
        /**	Mute or unmute a source with a boolean value. The index must be between 0 and the number of sources - 1.
         
            @param     index	The index of the source.
            @param     muted	The mute state.
         */
        void setMute(const unsigned int index, const bool muted);
        
        //! This method retrieve the azimuth of a source.
        /** Retrieve the azimuth of a source.
         
            @param     index	The index of the source.
            @return The azimuth of the source if the source exists, otherwise the function generates an error.
         */
        double getAzimuth(const unsigned int index) const
        {
            assert(index < m_number_of_sources);
            return m_azimuth[index];
        }
		
        //! This method retrieve the radius of a source.
        /** Retrieve the radius of a source.
         
            @param     index	The index of the source.
            @return The radius of the source if the source exists, otherwise the function generates an error.
         */
        double getRadius(const unsigned int index) const
        {
            assert(index < m_number_of_sources);
            if(m_wide[index] / ((double)(NUMBEROFLINEARPOINTS - 1) * m_number_of_harmonics) < 1)
                return m_wide[index] / ((double)(NUMBEROFLINEARPOINTS - 1) * m_number_of_harmonics);
            else
                return 1. / sqrt(m_gains[index]);
        }
        
		//! This method retrieve the mute or unmute state of a source.
        /**	Get the Mute state of a source.
         
            @param     index	The index of the source.
            @return    The mute state of the source if the source exists, otherwise the function generates an error.
            @see       setMute()
         */
        bool getMute(const unsigned int index) const
        {
            assert(index < m_number_of_sources);
            return m_muted[index];
        }
		
        
        //! This method performs the encoding with distance compensation with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding with distance compensation sample by sample. The inputs array contains the samples of the sources and the minimum size sould be the number of sources. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
         @param     inputs  The inputs array.
         @param     outputs The outputs array.
         */
        void process(const float* inputs, float* outputs);
        
        //! This method performs the encoding with distance compensation with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding with distance compensation sample by sample. The inputs array contains the samples of the sources and the minimum size sould be the number of sources. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     inputs  The inputs array.
            @param     outputs The outputs array.
         */
        void process(const double* inputs, double* outputs);
    };
}

#endif
