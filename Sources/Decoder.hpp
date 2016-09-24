/*
// Copyright (c) 2012-2015 Eliott Paris & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_DECODER_LIGHT
#define DEF_HOA_DECODER_LIGHT

#include "Encoder.hpp"
#include "Hrir.hpp"

namespace hoa
{
    //! The decoder class decodes a sound field in the harmonics domain through the planewaves domain.
    /** The decoder should be used to decode a set the harmonics domain to a set of planewaves for loudspeakers. There are three types of decoder. Regular for a perfect circle or sphere of loudspeakers. Irregular when the loudspeakers are not equally spaced on the circle or the sphere. Binaural for headphone restitution.
     */
    template <Dimension D, typename T> class Decoder : public Processor<D, T>::Harmonics, public Processor<D, T>::Planewaves
    {
    public:
        //! The decoder constructor.
        /**	The decoder constructor allocates and initialize the base classes.
         @param     order                   The order
         @param     numberOfPlanewaves      The number of channels.
         */
        Decoder(const size_t order, const size_t numberOfPlanewaves) hoa_noexcept;

        //! The destructor.
        /** The destructor free the memory.
         */
		virtual ~Decoder() = 0;

        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimum size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        virtual void process(const T* inputs, T* outputs) hoa_noexcept;

        //! This method computes the decoding matrix.
        /**	You should use this method after changing the position of the loudspeakers and/or calling the process method.
         @param vectorsize The vector size for binaural decoding.
         */
        virtual void computeRendering(const size_t vectorsize = 64);
    };
    
    
    //! The regular decoder class decodes a sound field in the harmonics domain through the planewaves domain for a perfect circle or sphere of loudspeakers.
    /** The regular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if more or equal to the number of harmonics plus one and when the loudspeakers are equally spaced on the circle or the sphere.
     */
    template <Dimension D, typename T> class DecoderRegular : public Decoder<D, T>
    {
    public:
        
        //! The regular constructor.
        /**	The regular constructor allocates and initialize the decoding matrix depending on a order of decomposition and a number of channels. The order must be at least 1 and the number of channels must be at least the number of harmonics.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        DecoderRegular(const size_t order, const size_t numberOfPlanewaves) hoa_noexcept;
        
        //! The destructor.
        /** The destructor free the memory.
         */
        virtual ~DecoderRegular() = 0;
        
        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimum size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to channels.
         */
        virtual void process(const T* inputs, T* outputs) hoa_noexcept hoa_override;
        
        //! This method computes the decoding matrix.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        virtual void computeRendering(const size_t vectorsize = 64) hoa_override;
    };
    
    //! The irregular decoder class decodes a sound field in the harmonics domain through the planewaves domain for a irregular circle or sphere (2d only).
    /** The irregular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if less than the number of harmonics plus one or when the loudspeakers are not equally spaced on the circle or the sphere.
     */
    template <Dimension D, typename T> class DecoderIrregular : public Decoder<D, T>
    {
    public:
        //! The irregular constructor.
        /**	The irregular constructor allocates and initialize the decoding matrix depending on a order of decomposition and a number of channels. The order must be at least 1.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        DecoderIrregular(const size_t order, const size_t numberOfPlanewaves) hoa_noexcept;
        
        //! The destructor.
        /** The destructor free the memory.
         */
        virtual ~DecoderIrregular() = 0;
        
        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimum size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        virtual void process(const T* inputs, T* outputs) hoa_noexcept hoa_override;
        
        //! This method computes the decoding matrix.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        virtual void computeRendering(const size_t vectorsize = 64) hoa_override;
        
    };
    
    //! The binaural decoder class decodes a sound field in the harmonics domain for headphones.
    /** The binaural decoder should be used to decode an ambisonic sound field for headphones. It decodes the sound field through the planewaves domain an convolves the results with HRTF from the IRCAM database.
     */
    template <Dimension D, typename T> class DecoderBinaural : public Decoder<Hoa2d, T>
    {
    public:
        //! The binaural decoder constructor.
        /**	The binaural decoder constructor allocates and initialize the member values to the decoding matrix depending on a order of decomposition and a number of channels. The order and the number of channels must be at least 1.
         @param     order				The order
         */
        DecoderBinaural(const size_t order);
        
        
        //! The binaural decoder destructor.
        /**	The binaural decoder destructor free the memory.
         */
        virtual ~DecoderBinaural() = 0;
        
        //! This method computes the decoding matrix.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        virtual void computeRendering(const size_t vectorsize = 64)  hoa_override;
        
        //! This method performs the binaural decoding and the convolution.
        virtual  void processBlock() hoa_noexcept;
        
        //! This method performs the binaural decoding.
        /**	You should use this method for not-in-place processing and performs the binaural decoding sample by sample. The inputs array contains the spherical harmonics samples : inputs[number of harmonics] and the outputs array contains the headphones samples : outputs[2].
         
         @param     inputs	The input samples.
         @param     outputs  The output array that contains samples destinated to channels.
         */
        virtual void process(const T* inputs, T* outputs) hoa_noexcept hoa_override;
    };

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    //! The decoder class decodes a sound field in the harmonics domain through the planewaves domain.
    /** The decoder should be used to decode a set the harmonics domain to a set of planewaves for loudspeakers. There are three types of decoder. Regular for a perfect circle or sphere of loudspeakers. Irregular when the loudspeakers are not equally spaced on the circle or the sphere. Binaural for headphone restitution.
     */
    template <typename T> class Decoder<Hoa2d, T> : public Processor<Hoa2d, T>::Harmonics, public Processor<Hoa2d, T>::Planewaves
    {
    public:
        enum Mode
        {
            RegularMode = 0,
            IrregularMode = 1,
            BinauralMode = 2
        };

        //! The decoder constructor.
        /**	The decoder constructor allocates and initialize the base classes.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        Decoder(const size_t order, const size_t numberOfPlanewaves) hoa_noexcept :
        Processor<Hoa2d, T>::Harmonics(order),
        Processor<Hoa2d, T>::Planewaves(numberOfPlanewaves)
        {
            ;
        }

        //! The destructor.
        /** The destructor free the memory.
         */
        virtual ~Decoder() hoa_default_f

        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimum size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        inline virtual void process(const T* inputs, T* outputs) hoa_noexcept = 0;

        //! This method computes the decoding matrix.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        virtual void computeRendering(const size_t vectorsize = 64) = 0;

        //! This method retrieves the mode of the decoder.
        /**	This method retrieves the mode of the decoder.
         @retun The mode of the decoder.
         */
        inline virtual Mode getMode() const hoa_noexcept {return RegularMode;}
    };

    //! The ambisonic regular decoder.
    /** The regular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if more or equal to the number of harmonics plus one and when the loudspeakers are equally spaced.
     */
    template <typename T> class DecoderRegular<Hoa2d, T> : public Decoder<Hoa2d, T>
    {
    private:
        T*  m_matrix;
    public:

        //! The regular constructor.
        /**	The regular constructor allocates and initialize the decoding matrix depending on a order of decomposition and a number of channels. The order must be at least 1 and the number of channels must be at least the number of harmonics.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        DecoderRegular(const size_t order, const size_t numberOfPlanewaves) hoa_noexcept : Decoder<Hoa2d, T>(order, numberOfPlanewaves)
        {
            m_matrix = Signal<T>::alloc(Decoder<Hoa2d, T>::getNumberOfPlanewaves() * Decoder<Hoa2d, T>::getNumberOfHarmonics());
            computeRendering();
        }

        //! The destructor.
        /** The destructor free the memory.
         */
        ~DecoderRegular()
        {
            Signal<T>::free(m_matrix);
        }

        //! This method retrieves the mode of the decoder.
        /**	This method retrieves the mode of the decoder.
         @retun The mode of the decoder.
         */
        inline typename Decoder<Hoa2d, T>::Mode getMode() const hoa_noexcept hoa_override {return Decoder<Hoa2d, T>::RegularMode;};

        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimum size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to channels.
         */
        inline void process(const T* inputs, T* outputs) hoa_noexcept hoa_override
        {
            Signal<T>::mul(Decoder<Hoa2d, T>::getNumberOfHarmonics(), Decoder<Hoa2d, T>::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }

        //! This method computes the decoding matrix.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        void computeRendering(const size_t vectorsize = 64) hoa_override
        {
            EncoderBasic<Hoa2d, T> encoder(Decoder<Hoa2d, T>::getDecompositionOrder());
            const T factor = 1. / (T)(Decoder<Hoa2d, T>::getDecompositionOrder() + 1.);
            for(size_t i = 0; i < Decoder<Hoa2d, T>::getNumberOfPlanewaves(); i++)
            {
                encoder.setAzimuth(Decoder<Hoa2d, T>::getPlanewaveAzimuth(i));
                encoder.process(&factor, m_matrix + i * Decoder<Hoa2d, T>::getNumberOfHarmonics());
                m_matrix[i * encoder.getNumberOfHarmonics()] = factor * 0.5;
            }
        }
    };

    
    //! The ambisonic irregular decoder.
    /** The irregular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if less than the number of harmonics plus one or when the loudspeakers are not equally spaced.
     */
    template <typename T> class DecoderIrregular<Hoa2d, T> : public Decoder<Hoa2d, T>
    {
    private:
        T*  m_matrix;
    public:

        //! The irregular constructor.
        /**	The irregular constructor allocates and initialize the decoding matrix depending on a order of decomposition and a number of channels. The order must be at least 1.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        DecoderIrregular(const size_t order, const size_t numberOfPlanewaves) hoa_noexcept : Decoder<Hoa2d, T>(order, numberOfPlanewaves)
        {
            m_matrix = Signal<T>::alloc(Decoder<Hoa2d, T>::getNumberOfPlanewaves() * Decoder<Hoa2d, T>::getNumberOfHarmonics());
            computeRendering();
        }

        //! The destructor.
        /** The destructor free the memory.
         */
        ~DecoderIrregular()
        {
            Signal<T>::free(m_matrix);
        }

        //! This method retrieves the mode of the decoder.
        /**	This method retrieves the mode of the decoder.
         @retun The mode of the decoder.
         */
        inline typename Decoder<Hoa2d, T>::Mode getMode() const hoa_noexcept hoa_override {return Decoder<Hoa2d, T>::IrregularMode;};

        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimum size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        inline void process(const T* inputs, T* outputs) hoa_noexcept hoa_override
        {
            Signal<T>::mul(Decoder<Hoa2d, T>::getNumberOfHarmonics(), Decoder<Hoa2d, T>::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }

        //! This method computes the decoding matrix.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        void computeRendering(const size_t vectorsize = 64)  hoa_override
        {
            EncoderBasic<Hoa2d, T> encoder(Decoder<Hoa2d, T>::getDecompositionOrder());
            Signal<T>::clear(Decoder<Hoa2d, T>::getNumberOfPlanewaves() * Decoder<Hoa2d, T>::getNumberOfHarmonics(), m_matrix);
            T* vector_harmonics = Signal<T>::alloc(Decoder<Hoa2d, T>::getNumberOfHarmonics());

            if(Decoder<Hoa2d, T>::getNumberOfPlanewaves() == 1)
            {
                const size_t nls = size_t(Decoder<Hoa2d, T>::getDecompositionOrder() + 1.);
                const T factor = 1. / (T)(nls);
                for(size_t i = 0; i <nls; i++)
                {
                    encoder.setAzimuth(T(i) * HOA_2PI / T(nls));
                    encoder.process(&factor, vector_harmonics);
                    vector_harmonics[0] = factor * 0.5;
                    Signal<T>::add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix);
                }
            }
            else
            {
                T smallest_distance = (T)HOA_2PI;
                std::vector<Planewave<Hoa2d, T> > channels;
                for(size_t i = 0; i < Decoder<Hoa2d, T>::getNumberOfPlanewaves(); i++)
                {
                    channels.push_back(Planewave<Hoa2d, T>(i, Math<T>::wrap_twopi(Decoder<Hoa2d, T>::getPlanewaveAzimuth(i)), 0.));
                }

                sort(channels.begin(), channels.end(), Planewave<Hoa2d, T>::sort_azimuth);

                {
                    const T current_angle   = channels[0].getAzimuth(0., 0., 0.);
                    const T previous_angle  = channels[channels.size() - 1].getAzimuth(0., 0., 0.);
                    const T previous_portion= (HOA_2PI - previous_angle) + current_angle;
                    if(smallest_distance > previous_portion)
                    {
                        smallest_distance = previous_portion;
                    }
                    //post("channel %i : %f", (int)channels[0].getIndex(), (float)(channels[0].getAzimuth() / HOA_2PI * 360.f));
                }
                for(size_t i = 1; i < channels.size(); i++)
                {
                    const T current_angle   = channels[i].getAzimuth(0., 0., 0.);
                    const T previous_angle  = channels[i-1].getAzimuth(0., 0., 0.);
                    const T previous_portion= current_angle - previous_angle;
                    if(smallest_distance > previous_portion)
                    {
                        smallest_distance = previous_portion;
                    }
                    //post("channel %i : %f", (int)channels[i].getIndex(), (float)(channels[i].getAzimuth() / HOA_2PI * 360.f));
                }
                //post("");

                if(smallest_distance > HOA_2PI / T(Decoder<Hoa2d, T>::getNumberOfHarmonics() + 1.))
                {
                    smallest_distance = HOA_2PI / T(Decoder<Hoa2d, T>::getNumberOfHarmonics() + 1.);
                }
                const size_t nvirtual = (size_t)ceil(HOA_2PI / smallest_distance);
                const T factor = 1. / (T)(nvirtual);

                //post("number of virtual %i", nvirtual);
                for(size_t i = 0; i < nvirtual; i++)
                {
                    const T angle = T(i) / T(nvirtual) * HOA_2PI;
                    //post("virtual %i :  %f", (int)i , (float)(angle / HOA_2PI * 360.f));
                    if(angle < channels[0].getAzimuth(0., 0., 0.))
                    {
                        const T portion = (HOA_2PI - channels[channels.size()-1].getAzimuth(0., 0., 0.)) + channels[0].getAzimuth(0., 0., 0.);

                        const T factor1 = (1. - ((channels[0].getAzimuth(0., 0., 0.) - angle) / portion)) * factor;
                        encoder.setAzimuth(angle);
                        encoder.process(&factor1, vector_harmonics);
                        vector_harmonics[0] = factor1 * 0.5;
                        Signal<T>::add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[0].getIndex() * Decoder<Hoa2d, T>::getNumberOfHarmonics());

                        const T factor2 = ((channels[0].getAzimuth(0., 0., 0.) - angle) / portion) * factor;
                        encoder.process(&factor2, vector_harmonics);
                        vector_harmonics[0] = factor2 * 0.5;
                        Signal<T>::add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[channels.size() - 1].getIndex() * Decoder<Hoa2d, T>::getNumberOfHarmonics());

                        //post("portion : %f", (float)portion / HOA_2PI * 360.f);
                        //post("channel %i (%f) : %f", (int)channels[channels.size()-1].getIndex(),
                        //     (float)(channels[channels.size()-1].getAzimuth() / HOA_2PI * 360.f),
                        //     (float)(factor2 / factor) * 100.f);
                        //post("channel %i (%f) : %f", (int)channels[0].getIndex(),
                        //     (float)(channels[0].getAzimuth() / HOA_2PI * 360.f),
                        //     (float)(factor1 / factor) * 100.f);
                    }
                    else if(angle >= channels[channels.size() - 1].getAzimuth(0., 0., 0.))
                    {
                        const T portion = (HOA_2PI - channels[channels.size()-1].getAzimuth(0., 0., 0.)) + channels[0].getAzimuth(0., 0., 0.);

                        const T factor1 = (1. - ((angle - channels[channels.size()-1].getAzimuth(0., 0., 0.)) / portion)) * factor;
                        encoder.setAzimuth(angle);
                        encoder.process(&factor1, vector_harmonics);
                        vector_harmonics[0] = factor1 * 0.5;
                        Signal<T>::add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[channels.size()-1].getIndex() * Decoder<Hoa2d, T>::getNumberOfHarmonics());

                        const T factor2 = ((angle - channels[channels.size()-1].getAzimuth(0., 0., 0.)) / portion) * factor;
                        encoder.process(&factor2, vector_harmonics);
                        vector_harmonics[0] = factor2 * 0.5;
                        Signal<T>::add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[0].getIndex() * Decoder<Hoa2d, T>::getNumberOfHarmonics());

                        //post("portion : %f", (float)portion / HOA_2PI * 360.f);
                        //post("channel %i (%f) : %f", (int)channels[channels.size()-1].getIndex(),
                        //     (float)(channels[channels.size()-1].getAzimuth() / HOA_2PI * 360.f),
                        //     (float)(factor1 / factor) * 100.f);
                        //post("channel %i (%f) : %f", (int)channels[0].getIndex(),
                        //     (float)(channels[0].getAzimuth() / HOA_2PI * 360.f),
                        //     (float)(factor2 / factor) * 100.f);
                    }
                    else
                    {
                        for(size_t j = 1; j < channels.size(); j++)
                        {
                            if(angle < channels[j].getAzimuth(0., 0., 0.) && angle >= channels[j-1].getAzimuth(0., 0., 0.))
                            {
                                const T portion = (channels[j].getAzimuth(0., 0., 0.) - channels[j-1].getAzimuth(0., 0., 0.));

                                const T factor1 = (1. - ((channels[j].getAzimuth(0., 0., 0.) - angle) / portion)) * factor;
                                encoder.setAzimuth(angle);
                                encoder.process(&factor1, vector_harmonics);
                                vector_harmonics[0] = factor1 * 0.5;
                                Signal<T>::add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[j].getIndex() * Decoder<Hoa2d, T>::getNumberOfHarmonics());

                                const T factor2 = ((channels[j].getAzimuth(0., 0., 0.) - angle) / portion) * factor;
                                encoder.process(&factor2, vector_harmonics);
                                vector_harmonics[0] = factor2 * 0.5;
                                Signal<T>::add(Decoder<Hoa2d, T>::getNumberOfHarmonics(), vector_harmonics, m_matrix + channels[j-1].getIndex() * Decoder<Hoa2d, T>::getNumberOfHarmonics());

                                //post("portion : %f", (float)portion / HOA_2PI * 360.f);
                                //post("channel %i (%f) : %f", (int)channels[j-1].getIndex(),
                                //     (float)(channels[j-1].getAzimuth() / HOA_2PI * 360.f),
                                //     (float)(factor2 / factor) * 100.f);
                                //post("channel %i (%f) : %f", (int)channels[j].getIndex(),
                                //     (float)(channels[j].getAzimuth() / HOA_2PI * 360.f),
                                //     (float)(factor1 / factor) * 100.f);

                                break;
                            }
                        }
                    }
                    //post("");
                }
                channels.clear();
            }
            Signal<T>::free(vector_harmonics);
        }
    };

    //! The ambisonic binaural decoder.
    /** The binaural decoder should be used to decode an ambisonic sound field for headphones.
     */
    template <typename T> class DecoderBinaural<Hoa2d, T> : public Decoder<Hoa2d, T>
    {
    private:
        size_t       m_vector_size;
        size_t       m_crop_size;
        T*          m_input;
        T*          m_result;
        T*          m_left;
        T*          m_right;

        void clear()
        {
            m_input  = Signal<T>::free(m_input);
            m_result = Signal<T>::free(m_result);
            m_left   = Signal<T>::free(m_left);
            m_right  = Signal<T>::free(m_right);
        }
    public:

        //! The binaural decoder constructor.
        /**	The binaural decoder constructor allocates and initialize the member values to the decoding matrix depending on a order of decomposition and a number of channels. The order and the number of channels must be at least 1.
         @param     order				The order
         */
        DecoderBinaural(const size_t order) hoa_noexcept : Decoder<Hoa2d, T>(order, 2),
        m_vector_size(0ul),
        m_input(hoa_nullptr),
        m_result(hoa_nullptr),
        m_left(hoa_nullptr),
        m_right(hoa_nullptr)
        {
            Decoder<Hoa2d, T>::setPlanewaveAzimuth(0, (T)(HOA_PI2*3.));
            Decoder<Hoa2d, T>::setPlanewaveAzimuth(1, (T)(HOA_PI2));
            setCropSize(0ul);
        }

        //! This method retrieves the mode of the decoder.
        /**	This method retrieves the mode of the decoder.
         @retun The mode of the decoder.
         */
        inline typename Decoder<Hoa2d, T>::Mode getMode() const hoa_noexcept hoa_override {return Decoder<Hoa2d, T>::BinauralMode;};

        //! The binaural decoder destructor.
        /**	The binaural decoder destructor free the memory.
         */
        ~DecoderBinaural() hoa_noexcept
        {
            clear();
        }

        //! This method sets the crop size of the responses.
        /**	This method sets the crop size of the responses.
         @param size The crop size.
         */
        inline void setCropSize(const size_t size) hoa_noexcept
        {
            if(!size || size > Hrir<Hoa2d, T>::getNumberOfRows())
                m_crop_size = Hrir<Hoa2d, T>::getNumberOfRows();
            else
                m_crop_size = size;
        }

        //! This method gets the crop size of the responses.
        /**	This method gets the crop size of the responses.
         @return The crop size.
         */
        inline size_t getCropSize() const hoa_noexcept
        {
            if(m_crop_size == Hrir<Hoa2d, T>::getNumberOfRows())
            {
                return 0;
            }
            else
            {
                return m_crop_size;
            }
        }

        //! This method computes the decoding matrix.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        void computeRendering(const size_t vectorsize = 64)  hoa_override
        {
            clear();
            m_vector_size  = vectorsize;
            m_input  = Signal<T>::alloc(Hrir<Hoa2d, T>::getNumberOfColumns() * m_vector_size);
            m_result = Signal<T>::alloc(Hrir<Hoa2d, T>::getNumberOfRows() * m_vector_size);
            m_left   = Signal<T>::alloc(Hrir<Hoa2d, T>::getNumberOfRows() + m_vector_size);
            m_right  = Signal<T>::alloc(Hrir<Hoa2d, T>::getNumberOfRows() + m_vector_size);
        }

    private:
        inline void processChannel(const T* harmonics, const T* response, T* vector, T* output) hoa_noexcept
        {
            const size_t l = Hrir<Hoa2d, T>::getNumberOfColumns();   // Harmonics size aka 11
            const size_t m = m_crop_size;      // Impulses size
            const size_t n = m_vector_size;    // Vector size
            Signal<T>::mul(m, n, l, response, harmonics, m_result);
            for(size_t i = 0; i < n; i ++)
            {
                Signal<T>::add(m, m_result + i, n, vector + i, 1ul);
            }
            Signal<T>::copy(m_vector_size, vector, output);
            Signal<T>::copy(m, vector + m_vector_size, vector);
            Signal<T>::clear(m_vector_size, vector + m);
        }
    public:

        //! This method performs the binaural decoding and the convolution.
        inline void processBlock(const T** inputs, T** outputs) hoa_noexcept
        {
            T* input = m_input;
            for(size_t i = 0; i < Hrir<Hoa2d, T>::getNumberOfColumns() && i < Decoder<Hoa2d, T>::getNumberOfHarmonics(); i++)
            {
                Signal<T>::copy(m_vector_size, inputs[i], input);
                input += m_vector_size;
            }
            processChannel(m_input, Hrir<Hoa2d, T>::getLeftMatrix(), m_left, outputs[0]);
            processChannel(m_input, Hrir<Hoa2d, T>::getRightMatrix(), m_right, outputs[1]);
        }

        inline void process(const T* inputs, T* outputs) hoa_noexcept hoa_override {}
    };


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    //! The decoder class decodes a sound field in the harmonics domain through the planewaves domain.
    /** The decoder should be used to decode a set the harmonics domain to a set of planewaves for loudspeakers. There are three types of decoder. Regular for a perfect circle or sphere of loudspeakers. Irregular when the loudspeakers are not equally spaced on the circle or the sphere. Binaural for headphone restitution.
     */
    template <typename T> class Decoder<Hoa3d, T> : public Processor<Hoa3d, T>::Harmonics, public Processor<Hoa3d, T>::Planewaves
    {
    public:

        enum Mode
        {
            RegularMode = 0,
            BinauralMode = 2
        };

        //! The decoder constructor.
        /**	The decoder constructor allocates and initialize the base classes.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        Decoder(const size_t order, const size_t numberOfPlanewaves) hoa_noexcept :
        Processor<Hoa3d, T>::Harmonics(order),
        Processor<Hoa3d, T>::Planewaves(numberOfPlanewaves)
        {
            ;
        }

        //! The destructor.
        /** The destructor free the memory.
         */
        virtual ~Decoder() {}

        //! This method retrieves the mode of the decoder.
        /**	This method retrieves the mode of the decoder.
         @retun The mode of the decoder.
         */
        inline virtual Mode getMode() const hoa_noexcept = 0;

        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimum size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        virtual void process(const T* inputs, T* outputs) hoa_noexcept = 0;

        //! This method computes the decoding matrix.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        virtual void computeRendering(const size_t vectorsize = 64) = 0;
    };

    //! The ambisonic regular decoder.
    /** The regular decoder should be used to decode an ambisonic sound field when the number of loudspeakers if more or equal to the number of harmonics plus one and when the loudspeakers are equally spaced.
     */
    template <typename T> class DecoderRegular<Hoa3d, T> : public Decoder<Hoa3d, T>
    {
    private:
        T*  m_matrix;
    public:

        //! The regular constructor.
        /**	The regular constructor allocates and initialize the decoding matrix depending on a order of decomposition and a number of channels. The order must be at least 1 and the number of channels must be at least the number of harmonics.
         @param     order				The order
         @param     numberOfPlanewaves     The number of channels.
         */
        DecoderRegular(const size_t order, const size_t numberOfPlanewaves) hoa_noexcept : Decoder<Hoa3d, T>(order, numberOfPlanewaves)
        {
            m_matrix = Signal<T>::alloc(Decoder<Hoa3d, T>::getNumberOfPlanewaves() * Decoder<Hoa3d, T>::getNumberOfHarmonics());
            computeRendering();
        }

        //! The regular destructor.
        /**	The regular destructor free the memory.
         */
        ~DecoderRegular()
        {
            Signal<T>::free(m_matrix);
        }
        
        //! This method retrieves the mode of the decoder.
        /**	This method retrieves the mode of the decoder.
         @retun The mode of the decoder.
         */
        inline typename Decoder<Hoa3d, T>::Mode getMode() const hoa_noexcept hoa_override {return Decoder<Hoa3d, T>::RegularMode;}

        //! This method performs the decoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimum size must be the number of channels.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to channels.
         */
        inline void process(const T* inputs, T* outputs) hoa_noexcept hoa_override
        {
            Signal<T>::mul(Decoder<Hoa3d, T>::getNumberOfHarmonics(), Decoder<Hoa3d, T>::getNumberOfPlanewaves(), inputs, m_matrix, outputs);
        }

        //! This method computes the decoding matrix.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        void computeRendering(const size_t vectorsize = 64)  hoa_override
        {
            EncoderBasic<Hoa3d, T> encoder(Decoder<Hoa3d, T>::getDecompositionOrder());
            const T factor = 1. / (T)(Decoder<Hoa3d, T>::getNumberOfPlanewaves());
            for(size_t i = 0; i < Decoder<Hoa3d, T>::getNumberOfPlanewaves(); i++)
            {
                encoder.setAzimuth(Decoder<Hoa3d, T>::getPlanewaveAzimuth(i));
                encoder.setElevation(Decoder<Hoa3d, T>::getPlanewaveElevation(i));
                encoder.process(&factor, m_matrix + i * Decoder<Hoa3d, T>::getNumberOfHarmonics());
                for(size_t j = 0; j < Decoder<Hoa3d, T>::getNumberOfHarmonics(); j++)
                {
                    const size_t l = Decoder<Hoa3d, T>::getHarmonicDegree(j);
                    if(encoder.getHarmonicOrder(j) == 0)
                    {
                        m_matrix[i * encoder.getNumberOfHarmonics() + j] *= (2. * l + 1.);
                    }
                    else
                    {
                        m_matrix[i * encoder.getNumberOfHarmonics() + j] *= T(2. * l + 1.) * 4. * HOA_PI;
                    }
                }
            }
        }
    };

    //! The ambisonic binaural decoder.
    /** The binaural decoder should be used to decode an ambisonic sound field for headphones.
     */
    template <typename T> class DecoderBinaural<Hoa3d, T> : public Decoder<Hoa3d, T>
    {
        size_t       m_vector_size;
        size_t       m_crop_size;
        T*          m_input;
        T*          m_result;
        T*          m_left;
        T*          m_right;

        void clear()
        {
            m_input  = Signal<T>::free(m_input);
            m_result = Signal<T>::free(m_result);
            m_left   = Signal<T>::free(m_left);
            m_right  = Signal<T>::free(m_right);
        }

    public:

        //! The binaural decoder constructor.
        /**	The binaural decoder constructor allocates and initialize the member values to the decoding matrix depending on a order of decomposition and a number of channels. The order and the number of channels must be at least 1.
         @param     order				The order
         */
        DecoderBinaural(const size_t order) : Decoder<Hoa3d, T>(order, 2),
        m_vector_size(0ul),
        m_input(hoa_nullptr),
        m_result(hoa_nullptr),
        m_left(hoa_nullptr),
        m_right(hoa_nullptr)
        {
            Decoder<Hoa3d, T>::setPlanewaveAzimuth(0, (T)(HOA_PI2*3.));
            Decoder<Hoa3d, T>::setPlanewaveAzimuth(1, (T)HOA_PI2);
            setCropSize(0ul);
        }

        //! This method retrieves the mode of the decoder.
        /**	This method retrieves the mode of the decoder.
         @retun The mode of the decoder.
         */
        inline typename Decoder<Hoa3d, T>::Mode getMode() const hoa_noexcept hoa_override {return Decoder<Hoa3d, T>::BinauralMode;};

        //! The binaural decoder destructor.
        /**	The binaural decoder destructor free the memory.
         */
        ~DecoderBinaural()
        {
            clear();
        }

        //! This method sets the crop size of the responses.
        /**	This method sets the crop size of the responses.
         @param size The crop size.
         */
        inline void setCropSize(const size_t size) hoa_noexcept
        {
            if(!size || size > Hrir<Hoa3d, T>::getNumberOfRows())
            {
                m_crop_size = Hrir<Hoa3d, T>::getNumberOfRows();
            }
            else
            {
                m_crop_size = size;
            }
        }

        //! This method gets the crop size of the responses.
        /**	This method gets the crop size of the responses.
         @return The crop size.
         */
        inline size_t getCropSize() const hoa_noexcept
        {
            if(m_crop_size == Hrir<Hoa3d, T>::getNumberOfRows())
            {
                return 0;
            }
            else
            {
                return m_crop_size;
            }
        }

        //! This method computes the decoding matrix.
        /**	You should use this method after changing the position of the loudspeakers.
         @param vectorsize The vector size for binaural decoding.
         */
        void computeRendering(const size_t vectorsize = 64) hoa_override
        {
            clear();
            m_vector_size  = vectorsize;
            m_input  = Signal<T>::alloc(Hrir<Hoa3d, T>::getNumberOfColumns() * m_vector_size);
            m_result = Signal<T>::alloc(Hrir<Hoa3d, T>::getNumberOfRows() * m_vector_size);
            m_left   = Signal<T>::alloc(Hrir<Hoa3d, T>::getNumberOfRows() + m_vector_size);
            m_right  = Signal<T>::alloc(Hrir<Hoa3d, T>::getNumberOfRows() + m_vector_size);
        }

    private:
        inline void processChannel(const T* harmonics, const T* response, T* vector, T* output) hoa_noexcept
        {
            const size_t l = Hrir<Hoa3d, T>::getNumberOfColumns();   // Harmonics size aka 11
            const size_t m = m_crop_size;      // Impulses size
            const size_t n = m_vector_size;    // Vector size
            Signal<T>::mul(m, n, l, response, harmonics, m_result);
            for(size_t i = 0; i < n; i ++)
            {
                Signal<T>::add(m, m_result + i, n, vector + i, 1ul);
            }

            Signal<T>::copy(m_vector_size, vector, output);
            Signal<T>::copy(m, vector + m_vector_size, vector);
            Signal<T>::clear(m_vector_size, vector + m);
        }
    public:

        //! This method performs the binaural decoding and the convolution.
        inline void processBlock(const T** inputs, T** outputs) hoa_noexcept
        {
            T* input = m_input;
            for(size_t i = 0; i < Hrir<Hoa3d, T>::getNumberOfColumns() && i < Decoder<Hoa3d, T>::getNumberOfHarmonics(); i++)
            {
                Signal<T>::copy(m_vector_size, inputs[i], input);
                input += m_vector_size;
            }
            processChannel(m_input, Hrir<Hoa3d, T>::getLeftMatrix(), m_left, outputs[0]);
            processChannel(m_input, Hrir<Hoa3d, T>::getRightMatrix(), m_right, outputs[1]);
        }

        inline void process(const T* inputs, T* outputs) hoa_noexcept hoa_override {}

    };

}

#endif
