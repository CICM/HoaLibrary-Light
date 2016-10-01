/*
// Copyright (c) 2012-2016 CICM - Universite Paris 8 - Labex Arts H2H.
// Authors :
// 2012: Pierre Guillot, Eliott Paris & Julien Colafrancesco.
// 2012-2015: Pierre Guillot & Eliott Paris.
// 2015: Pierre Guillot & Eliott Paris & Thomas Le Meur (Light version)
// 2016: Pierre Guillot & Eliott Paris.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_DECODER_LIGHT
#define DEF_HOA_DECODER_LIGHT

#include "Hoa_Encoder.hpp"
#include "Hoa_Hrir.hpp"

namespace hoa
{
    //! The decoder class decodes a sound field in the harmonics domain through the planewaves domain.
    /** The decoder should be used to decode a set the harmonics domain to a set of planewaves for loudspeakers. There are three types of decoder. Regular for a perfect circle or sphere of loudspeakers. Irregular when the loudspeakers are not equally spaced on the circle or the sphere. Binaural for headphone restitution.
     */
    template <Dimension D, typename T> class Decoder : public ProcessorHarmonics<D, T>, public ProcessorPlanewaves<D, T>
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
        virtual void prepare(const size_t vectorsize = 64);
    };


    //! @brief The class decodes a sound field from the harmonics domain to a set of planewaves.
    //! @details The regular decoder should be used to decode from a set of signals associated
    //! to harmonics to a set signals associated to plane waves. The plane waves can match to
    //! loudspeakers or microphones. The decoding amounts to compute for each plane wave the
    //! dot product of the signal associated to the harmonics with the harmonic coefficients
    //! of a unitary signal in the direction of the plane wave.
    //! \f[S(\theta, \varphi)=\sum_{l=0}^{N}\sum_{m=-l}^{l}S_{l,m}Y_{l,m}^{\frac{1}{N+1}}(\theta,\varphi) \f]
    //! with \f$S_{l,m}\f$, the signals associated to the harmonics to decode, \f$Y_{l,m}^{\frac{1}{N+1}}\f$
    //! the harmonic coefficients of a unitary signal \f$\frac{1}{N+1}\f$ in the direction of
    //! the plane wave, \f$N\f$ the order of decomposition, \f$l\f$ the degree, \f$m\f$ the
    //! azimuthal order, \f$\theta\f$ the azimuth and \f$\varphi\f$ the elevation in radian of
    //! the plane wave.\n
    //! The plane waves must be equally spaced on a circle or a sphere and the number of plane
    //! waves must be at least the number of harmonics.
    template <Dimension D, typename T> class DecoderRegular : public Decoder<D, T>
    {
    public:

        //! @brief The constructor.
        //! @param order The order
        //! @param nplws The number of channels.
        DecoderRegular(size_t order, size_t nplws) : Decoder<D, T>(order, nplws)
        {
            m_matrix = new T[Decoder<D, T>::getNumberOfPlanewaves() * Decoder<D, T>::getNumberOfHarmonics()];
            prepare();
        }

        //! @brief The destructor.
        ~DecoderRegular() {
            delete [] m_matrix;
        }

        //! @brief The method performs the decoding of the harmonics signal.
        //! @details The input pointer must be the harmonics signal to decoder and the outputs
        //! array contains the plane waves samples thus the minimum size of the array must
        //! be the number of harmonics and the number of plane waves.
        //! @param inputs  The inputs array.
        //! @param outputs The outputs array.
        void process(const T* inputs, T* outputs) hoa_noexcept hoa_override
        {
            const size_t nharm = Decoder<D, T>::getNumberOfHarmonics();
            const size_t nplws = Decoder<D, T>::getNumberOfPlanewaves();
            const T* matrix = m_matrix;
            for(size_t i = 0ul; i < nplws; i++)
            {
                T result = 0;
                const T* in1 = inputs;
                for(size_t j = nharm>>3; j; --j, in1 += 8, matrix += 8)
                {
                    result += in1[0] * matrix[0]; result += in1[1] * matrix[1];
                    result += in1[2] * matrix[2]; result += in1[3] * matrix[3];
                    result += in1[4] * matrix[4]; result += in1[5] * matrix[5];
                    result += in1[6] * matrix[6]; result += in1[7] * matrix[7];
                }
                for(size_t j = nharm&7; j; --j, in1++, matrix++)
                {
                    result += in1[0] * matrix[0];
                }
                outputs[i] = result;
            }
        }

        //! @brief Prepare the decoder for processing.
        void prepare(const size_t vectorsize = 64) hoa_override
        {
            const size_t order = Decoder<D, T>::getDecompositionOrder();
            const size_t nharm = Decoder<D, T>::getNumberOfHarmonics();
            const size_t nplws = Decoder<D, T>::getNumberOfPlanewaves();
            
            Encoder<D, T> encoder(order);
            const T factor = T(1.) / T(order * 2 + 1);// / ((D == Hoa2d) ? T(order + 1) : T(1.));
            for(size_t i = 0; i < nplws; i++)
            {
                encoder.setAzimuth(Decoder<D, T>::getPlanewaveAzimuth(i));
                encoder.setElevation(Decoder<D, T>::getPlanewaveElevation(i));
                encoder.process(&factor, m_matrix + i * nharm);
                for(size_t j = 0; j < nharm; ++j)
                {
                    const size_t degree = Decoder<D, T>::getHarmonicDegree(j);
                     m_matrix[i * nharm + j] *= (T(2) * T(degree) + T(1));
                }
            }
        }
        
        //! @brief Return the type of the decoder.
        inline typename Decoder<D, T>::Mode getMode() const hoa_noexcept {return Decoder<D, T>::RegularMode;}
        
    private:
        T*  m_matrix;
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
        virtual void prepare(const size_t vectorsize = 64) hoa_override;

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
        virtual void prepare(const size_t vectorsize = 64)  hoa_override;

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
    template <typename T> class Decoder<Hoa2d, T> : public ProcessorHarmonics<Hoa2d, T>, public ProcessorPlanewaves<Hoa2d, T>
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
        ProcessorHarmonics<Hoa2d, T>(order),
        ProcessorPlanewaves<Hoa2d, T>(numberOfPlanewaves)
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
        virtual void prepare(const size_t vectorsize = 64) = 0;

        //! This method retrieves the mode of the decoder.
        /**	This method retrieves the mode of the decoder.
         @retun The mode of the decoder.
         */
        inline virtual Mode getMode() const hoa_noexcept {return RegularMode;}
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
            prepare();
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
        void prepare(const size_t vectorsize = 64)  hoa_override
        {
            Encoder<Hoa2d, T> encoder(Decoder<Hoa2d, T>::getDecompositionOrder());
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
                    T nazimtuh = Decoder<Hoa2d, T>::getPlanewaveAzimuth(i);
                    while(nazimtuh < T(0.)) {
                        nazimtuh += T(HOA_2PI);
                    }
                    while(nazimtuh >= T(HOA_2PI)) {
                        nazimtuh -= T(HOA_2PI);
                    }
                    channels.push_back(Planewave<Hoa2d, T>(i, nazimtuh, 0.));
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
        void prepare(const size_t vectorsize = 64)  hoa_override
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
    template <typename T> class Decoder<Hoa3d, T> : public ProcessorHarmonics<Hoa3d, T>, public ProcessorPlanewaves<Hoa3d, T>
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
        ProcessorHarmonics<Hoa3d, T>(order),
        ProcessorPlanewaves<Hoa3d, T>(numberOfPlanewaves)
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
        inline virtual Mode getMode() const hoa_noexcept { return RegularMode; }

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
        virtual void prepare(const size_t vectorsize = 64) = 0;
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
        void prepare(const size_t vectorsize = 64) hoa_override
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
