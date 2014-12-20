/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_CONVERTER__
#define __DEF_HOA_CONVERTER__

#include "../Hoa2D/Hoa2d.h"
#include "../Hoa3D/Hoa3D.h"

namespace Hoa
{
    //! The converter class.
    /** The converter classe can be used to from/to Hoa and another format.
     */
    class Converter
    {
    public:

        enum Mode
        {
            HOA2D = 0,   /**< The Hoa 2D format */
            HOA3D = 1,   /**< The Hoa 3D format */
            FUMA_2D = 2, /**< The 2D Furse-Malham component ordering and normalization */
            FUMA_3D = 3, /**< The 3D Furse-Malham component ordering and normalization */
            ACN_N2D = 4, /**< The 2D ACN component ordering and N2D normalization */
            ACN_SN2D = 5, /**< The 2D ACN component ordering and SN2D normalization */
            ACN_N3D = 6, /**< The 3D ACN component ordering and N3D normalization */
            ACN_SN3D = 7, /**< The 3D ACN component ordering and SN3D normalization */
        };
    protected:

        unsigned int    m_order;
        Mode            m_input_mode;
        Mode            m_output_mode;

        unsigned int    m_number_of_inputs_harmonics;
        unsigned int    m_number_of_outputs_harmonics;
        
    public:
        
        //! The converter constructor.
        /** The converter constructor computes the weights and the ordering of the harmonics to convert from/to Hoa and another format. The order must be at least 1.
         
            @param     order	The order.
         */
        Converter(unsigned int order, Mode inputMode, Mode outputMode);
        
        //! The converter destructor.
        /**	The converter destructor free the memory.
         */
        ~Converter();
        
        //! Retrieve the decomposition order.
        /** Retrieve the decomposition order of an ambisonic class.
         */
        inline unsigned int getDecompositionOrder() const {return m_order;};
        
        //! Retrieve the number of inputs harmonics.
        /** Retrieve the number of inputs harmonics.
         */
        inline unsigned int getNumberOfInputHarmonics() const {return m_number_of_inputs_harmonics;};
        
        //! Retrieve the number of outputs harmonics.
        /** Retrieve the number of outputs harmonics.
         */
        inline unsigned int getNumberOfOutputHarmonics() const {return m_number_of_outputs_harmonics;};
        
        //! Retrieve the argument of an harmonic.
        /** The argument of an harmonic is in the range -band to band. The harmonics are sorted by their bands, from 0 to the decomposition order. In each band contains 2 * band + 1 harmonics, sorted by their arguments in the range -band to band. The harmonic input and output arrays in process method of ambisonic classes must have this configuration. 
			For the first bands, the harmonics arrangement is h[0, 0] h[1, 0] h[1, -1] h[1, 1] h[2, 0] h[2, -1] h[2, 1] h[2, -2] h[2, 2] etc.
			with h[band, argument].
         
            @param     index	The global index of an harmonic.
            @return    The method returns the argument of the harmonic if the harmonic exists, otherwise the function generates an error.
            @see       getHarmonicDegree()
            @see       getHarmonicName()
         */
        inline int getHarmonicOrder(const unsigned int index) const
        {
            assert(index < m_number_of_harmonics);
            return m_harmonics_arguments[index];
        };
        
        //! Retrieve the band of an harmonic.
        /** The bands of the harmonics are in the range 0 to the decomposition order. Each band contains 2 * band + 1 harmonics in the range -band to band. The harmonic input and output arrays in process method of ambisonic classes must have this configuration.
			For the first bands, the harmonics arrangement is h[0, 0] h[1, 0] h[1, -1] h[1, 1] h[2, 0] h[2, -1] h[2, 1] h[2, -2] h[2, 2] etc.
			with h[band, argument].
         
            @param     index	The global index of an harmonic.
            @return    The method returns the band of the harmonic if the harmonic exists, otherwise the function generates an error.
            @see       getHarmonicOrder()
            @see       getHarmonicName()
         */
        inline unsigned int getHarmonicDegree(const unsigned int index) const
        {
            assert(index < m_number_of_harmonics);
            return m_harmonics_bands[index];
        };
        
        //! Retrieve a name for an harmonic.
        /** Retrieve a name for an harmonic in a std::string format that will be "harmonic band argument".
         
            @param     index	The global index of an harmonic.
            @return    The method returns a name for the harmonic that contains its band and its argument if the harmonic exists, otherwise the function generates an error.
            @see       getHarmonicDegree()
            @see       getHarmonicOrder()
         */
        inline std::string getHarmonicName(const unsigned int index) const
        {
            assert(index < m_number_of_harmonics);
            return "Harmonic " + int_to_string(getHarmonicDegree(index)) + " " + int_to_string(getHarmonicOrder(index));
        };
        
        //! This method performs the convertion with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the convertion sample by sample. The inputs array contains the harmonics samples and the minimum size must be the number of inputs harmonics. The outputs array contains the harmonics samples and the minimum size must be the number of outputs harmonics.
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void process(const float* inputs, float* outputs);
        
        //! This method performs the encoding with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The inputs array contains the harmonics samples and the minimum size must be the number of inputs harmonics. The outputs array contains the harmonics samples and the minimum size must be the number of outputs harmonics.
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void process(const double* inputs, double* outputs);
    };
}

#endif


