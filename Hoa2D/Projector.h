/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_PROJECTOR
#define DEF_HOA_2D_PROJECTOR

#include "Ambisonic.h"
#include "Planewaves.h"
#include "Encoder.h"

namespace Hoa2D
{
    //! The ambisonic projector.
    /** The projector should be used to projects the circular harmonics on a circle and to give access to the planewaves domain. The projection is similar to the decoding exept that the circle discretization cannot be defined be the user. Then the number of channels (or planewaves) must be a least the number of harmonics, the first angle is 0 radian and the angular distances between the channels are equals.
     */
    class Projector : public Ambisonic, public Planewaves
    {
        
    private:
        double*         m_projector_matrix_double;
        float*          m_projector_matrix_float;        
    public:
        
        //! The projector constructor.
        /**	The projector constructor allocates and initialize the member values to project the circular harmonics on a circle. The order must be at least 1 and the number of channels must be a least the number of harmonics.
         
            @param     order				The order.
            @param     numberOfChannels     The number of channels.
         */
		Projector(unsigned int order, unsigned int numberOfChannels);
		
        //! The projector destructor.
        /**	The projector destructor free the memory.
         */
		~Projector();
        
		//! This method performs the projection with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the projection sample by sample. The inputs array contains the circular harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels (or planewaves) samples and the minimum size must be a least the number of channels.
         
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
		void process(const float* inputs, float* outputs);
		
		//! This method performs the projection with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the projection sample by sample. The inputs array contains the circular harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels (or planewaves) samples and the minimum size must be a least the number of channels.
         
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
		void process(const double* inputs, double* outputs);
    };
}


#endif


