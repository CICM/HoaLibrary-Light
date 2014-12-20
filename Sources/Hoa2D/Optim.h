/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_OPTIM
#define DEF_HOA_2D_OPTIM

#include "Ambisonic.h"

namespace Hoa2D
{
    //! The ambisonic optimization.
    /** The optimization should be used to optimize the ambisonic sound field. There are 3 optimization modes, Basic (no optimization), MaxRe (energy vector optimization) and InPhase (energy and velocity vector optimization). Basic has no effect, it should be used with a perfect ambisonic channels, arrengement where all the channels are to equal distance on a circle, and for a listener placed at the perfect center of the circle. MaxRe should be used for auditory confined to the center of the circle. InPhase should be used when the auditory covers the entire channels area and when the channels arragement is not a perfect circle or when the channels are not to equal distance. Note that the optimizations decrease the precision of the sound field restitution thus it can be compared to particular cases of the fractional orders.
     */    class Optim : public Ambisonic
    {
    public:
        enum Mode
        {
            Basic   = 0,	/**< Basic Optimization     */
            MaxRe   = 1,	/**< Max Re Optimization    */
            InPhase = 2     /**< In Phase Optimization  */
        };
        
    private:
        
        Mode            m_mode;
        double*         m_harmonics;
        
    public:
        
        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         
            @param     order	The order.
         */
        Optim(unsigned int order, Mode mode = Basic);
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~Optim();
        
        //! Retrieve the optimization mode.
        /** Retrieve the optimization mode : Basic, MaxRe or InPhase.
         
            @return The method returns the optimization mode.
         */
        Mode getMode() const
        {
            return m_mode;
        }
        
        //! This method set the optimization mode.
        /**	The mode should be one of the 3 optimization modes, Basic, MaxRe or InPhase.
         
            @param     mode The optimization mode.
         */
        void setMode(Mode mode);
        
        //! This method performs the optimization with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the optimization sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
        void process(const float* inputs, float* outputs);
        
        //! This method performs the optimization with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the optimization sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
        void process(const double* inputs, double* outputs);
    };
}

#endif



