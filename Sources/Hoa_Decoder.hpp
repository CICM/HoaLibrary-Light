/*
 // Copyright (c) 2012-2017 CICM - Universite Paris 8 - Labex Arts H2H.
 // Authors :
 // 2012: Pierre Guillot, Eliott Paris & Julien Colafrancesco.
 // 2012-2015: Pierre Guillot & Eliott Paris.
 // 2015: Pierre Guillot & Eliott Paris & Thomas Le Meur (Light version)
 // 2016-2017: Pierre Guillot.
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#pragma once

#include "Hoa_Encoder.hpp"
#include "Hoa_Hrir.hpp"

#include <Eigen/Dense>

namespace hoa
{
    // ================================================================================ //
    // DECODER //
    // ================================================================================ //
    
    //! @brief The decoder class decodes a sound field in the harmonics domain through the planewaves domain.
    //! @details The decoder should be used to decode a set of harmonics to a set of planewaves for loudspeakers.
    template <Dimension D, typename T>
    class Decoder
    : public ProcessorHarmonics<D, T>
    , public ProcessorPlanewaves<D, T>
    {
    public:
        
        //! @brief Decoding mode :
        //! @details There are three types of decoders :
        //! - Regular is for a perfect circle or sphere of loudspeakers.
        //! - Irregular is suitable when the loudspeakers are not equally distibuted on the circle or the sphere.
        //! - Binaural is for headphone restitution.
        enum Mode { RegularMode = 0, IrregularMode, BinauralMode };
        
        //! @brief The decoder constructor.
        //! @details Allocates and initialize the base classes.
        //! @param order                   The order
        //! @param numberOfPlanewaves      The number of channels.
        Decoder(const size_t order, const size_t numberOfPlanewaves) noexcept;
        
        //! @brief The destructor.
        virtual ~Decoder() = 0;
        
        //! @brief This method performs the decoding.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! The inputs array contains the spherical harmonics samples (minimum is the number of harmonics)
        //! The outputs array contains the channels samples (minimum is the number of channels).
        //! @param inputs  The input array that contains the samples of the harmonics.
        //! @param outputs The output array that contains samples destinated to the channels.
        virtual void process(const T* inputs, T* outputs) noexcept;
        
        //! @brief This method computes the decoding matrix.
        //! @details You should use this method after changing the position of the loudspeakers
        //! and/or calling the process method.
        //! @param vectorsize The vector size for binaural decoding.
        virtual void prepare(const size_t vectorsize = 64);
    };
    
    // ================================================================================ //
    // DECODER REGULAR //
    // ================================================================================ //
    
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
    //! the plane wave.<br>
    //! The plane waves must be equally spaced on a circle or a sphere and the number of plane
    //! waves must be at least the number of harmonics.
    template <Dimension D, typename T>
    class DecoderRegular
    : public Decoder<D, T>
    {
    public:
        
        //! @brief The constructor.
        //! @param order The order
        //! @param nplws The number of channels.
        DecoderRegular(size_t order, size_t nplws)
        : Decoder<D, T>(order, nplws)
        , m_matrix(Decoder<D, T>::getNumberOfPlanewaves()
                   * Decoder<D, T>::getNumberOfHarmonics())
        {
            prepare();
        }
        
        //! @brief The destructor.
        ~DecoderRegular() = default;
        
        //! @brief The method performs the decoding of the harmonics signal.
        //! @details The input pointer must be the harmonics signal to decoder and the outputs
        //! array contains the plane waves samples thus the minimum size of the array must
        //! be the number of harmonics and the number of plane waves.
        //! @param inputs  The inputs array.
        //! @param outputs The outputs array.
        void process(const T* inputs, T* outputs) noexcept override
        {
            const size_t nharm = Decoder<D, T>::getNumberOfHarmonics();
            const size_t nplws = Decoder<D, T>::getNumberOfPlanewaves();
            const T* matrix = m_matrix.data();
            
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
        void prepare(const size_t vectorsize = 64) override
        {
            hoa_unused(vectorsize);
            
            const size_t order = Decoder<D, T>::getDecompositionOrder();
            const size_t nharm = Decoder<D, T>::getNumberOfHarmonics();
            const size_t nplws = Decoder<D, T>::getNumberOfPlanewaves();
            
            Encoder<D, T> encoder(order);
            const T factor = T(1.) / T(order * 2 + 1);// / ((D == Hoa2d) ? T(order + 1) : T(1.));
            for(size_t i = 0; i < nplws; i++)
            {
                encoder.setAzimuth(Decoder<D, T>::getPlanewaveAzimuth(i));
                encoder.setElevation(Decoder<D, T>::getPlanewaveElevation(i));
                encoder.process(&factor, m_matrix.data() + i * nharm);
                if(D == Hoa2d)
                {
                    m_matrix[i * nharm] *= T(0.5);
                }
                else
                {
                    for(size_t j = 0; j < nharm; ++j)
                    {
                        const size_t degree = Decoder<D, T>::getHarmonicDegree(j);
                        m_matrix[i * nharm + j] *= (T(2) * T(degree) + T(1));
                    }
                }
            }
        }
        
        //! @brief Returns the type of the decoder.
        inline typename Decoder<D, T>::Mode getMode() const noexcept override
        {
            return Decoder<D, T>::RegularMode;
        }
        
    private:
        
        std::vector<T> m_matrix;
    };
    
    // ================================================================================ //
    // DECODER IRREGULAR //
    // ================================================================================ //
    
    //! @brief The irregular decoder class decodes a sound field in the harmonics domain
    //! through the planewaves domain for an irregular circle or sphere (2d only).
    //! @details The irregular decoder should be used to decode an ambisonic sound field
    //! when the number of loudspeakers is less than the number of harmonics plus one
    //! or when the loudspeakers are not equally distributed on the circle or the sphere.
    template <Dimension D, typename T>
    class DecoderIrregular {};
    
    // ================================================================================ //
    // DECODER BINAURAL //
    // ================================================================================ //
    
    //! @brief The binaural decoder decodes a sound field in the harmonics domain for headphones.
    //! @details It decodes the sound field through the planewaves domain
    //! and convolves the results with HRIRs.
    template <Dimension D, typename T, typename HrirType>
    class DecoderBinaural {};
    
    // ================================================================================ //
    // DECODER 2D //
    // ================================================================================ //
    
    //! @brief Decodes a sound field in the harmonics domain through the planewaves domain.
    //! @details There are three types of decoder.
    //! - Regular for a perfect circle of loudspeakers.
    //! - Irregular when the loudspeakers are not equally spaced on the circle.
    //! - Binaural for headphone restitution.
    template <typename T>
    class Decoder<Hoa2d, T>
    : public ProcessorHarmonics<Hoa2d, T>
    , public ProcessorPlanewaves<Hoa2d, T>
    {
    public:
        
        enum Mode { RegularMode = 0, IrregularMode, BinauralMode };
        
        //! @brief The decoder constructor.
        //! @param order The order
        //! @param channels The number of channels.
        Decoder(const size_t order, const size_t channels)
        : ProcessorHarmonics<Hoa2d, T>(order)
        , ProcessorPlanewaves<Hoa2d, T>(channels)
        {}
        
        //! @brief Destructor.
        virtual ~Decoder() = default;
        
        //! @brief This method performs the decoding.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! @param inputs  The input array that contains the samples of the harmonics.
        //! @param outputs The output array that contains samples destinated to the channels.
        inline virtual void process(const T* inputs, T* outputs) noexcept = 0;
        
        //! @brief This method computes the decoding matrix.
        //! @details You should use this method after changing the position of the loudspeakers.
        //! @param vectorsize The vector size for binaural decoding.
        virtual void prepare(const size_t vectorsize = 64) = 0;
        
        //! @brief This method retrieves the mode of the decoder.
        //! @details    This method retrieves the mode of the decoder.
        //! @return The mode of the decoder.
        inline virtual Mode getMode() const noexcept {return RegularMode;}
    };
    
    // ================================================================================ //
    // DECODER 2D IRREGULAR //
    // ================================================================================ //
    
    //! @brief The ambisonic irregular decoder.
    //! @details The irregular decoder should be used to decode an ambisonic sound field
    //! when the number of loudspeakers if less than the number of harmonics plus one
    //! or when the loudspeakers are not equally spaced.
    template <typename T>
    class DecoderIrregular<Hoa2d, T>
    : public Decoder<Hoa2d, T>
    {
    private:
        T*  m_matrix;
    public:
        
        //! @brief Constructor.
        //! @param order The order
        //! @param channels The number of channels.
        DecoderIrregular(const size_t order, const size_t channels)
        : Decoder<Hoa2d, T>(order, channels)
        {
            m_matrix = Signal<T>::alloc(Decoder<Hoa2d, T>::getNumberOfPlanewaves()
                                        * Decoder<Hoa2d, T>::getNumberOfHarmonics());
            prepare();
        }
        
        //! @brief Destructor.
        ~DecoderIrregular()
        {
            Signal<T>::free(m_matrix);
        }
        
        //! @brief Returns the mode of the decoder.
        //! @return The mode of the decoder.
        inline typename Decoder<Hoa2d, T>::Mode getMode() const noexcept override
        {
            return Decoder<Hoa2d, T>::IrregularMode;
        }
        
        //! @brief This method performs the decoding.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! @param inputs The input array that contains the samples of the harmonics.
        //! @param outputs The output array that contains samples destinated to the channels.
        inline void process(const T* inputs, T* outputs) noexcept override
        {
            Signal<T>::mul(Decoder<Hoa2d, T>::getNumberOfHarmonics(),
                           Decoder<Hoa2d, T>::getNumberOfPlanewaves(),
                           inputs, m_matrix, outputs);
        }
        
        //! @brief This method computes the decoding matrix.
        //! @details You should use this method after changing the position of the loudspeakers.
        //! @param vectorsize The vector size for binaural decoding.
        void prepare(const size_t vectorsize = 64)  override
        {
            hoa_unused(vectorsize);
            
            const auto order = Decoder<Hoa2d, T>::getDecompositionOrder();
            const auto harmonics = Decoder<Hoa2d, T>::getNumberOfHarmonics();
            const auto planewaves = Decoder<Hoa2d, T>::getNumberOfPlanewaves();
            
            Encoder<Hoa2d, T> encoder(order);
            Signal<T>::clear(planewaves * harmonics, m_matrix);
            
            T* vector_harmonics = Signal<T>::alloc(harmonics);
            
            if(planewaves == 1)
            {
                const size_t nls = size_t(order + 1);
                const T factor = static_cast<T>(1.) / static_cast<T>(nls);
                for(size_t i = 0; i <nls; i++)
                {
                    encoder.setAzimuth(T(i) * HOA_2PI / T(nls));
                    encoder.process(&factor, vector_harmonics);
                    vector_harmonics[0] = factor * 0.5;
                    Signal<T>::add(harmonics, vector_harmonics, m_matrix);
                }
            }
            else
            {
                T smallest_distance = static_cast<T>(HOA_2PI);
                std::vector<Planewave<Hoa2d, T> > channels;
                for(size_t i = 0; i < planewaves; i++)
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
                
                sort(channels.begin(), channels.end(), Planewave<Hoa2d, T>::compare_azimuth);
                
                {
                    const T current_angle   = channels[0].getAzimuth();
                    const T previous_angle  = channels[channels.size() - 1].getAzimuth();
                    const T previous_portion= (HOA_2PI - previous_angle) + current_angle;
                    if(smallest_distance > previous_portion)
                    {
                        smallest_distance = previous_portion;
                    }
                    //post("channel %i : %f", (int)channels[0].getIndex(), (float)(channels[0].getAzimuth() / HOA_2PI * 360.f));
                }
                for(size_t i = 1; i < channels.size(); i++)
                {
                    const T current_angle   = channels[i].getAzimuth();
                    const T previous_angle  = channels[i-1].getAzimuth();
                    const T previous_portion= current_angle - previous_angle;
                    if(smallest_distance > previous_portion)
                    {
                        smallest_distance = previous_portion;
                    }
                    //post("channel %i : %f", (int)channels[i].getIndex(), (float)(channels[i].getAzimuth() / HOA_2PI * 360.f));
                }
                //post("");
                
                if(smallest_distance > HOA_2PI / T(harmonics + 1.))
                {
                    smallest_distance = HOA_2PI / T(harmonics + 1.);
                }
                const size_t nvirtual = static_cast<size_t>(std::ceil(static_cast<T>(HOA_2PI) / smallest_distance));
                const T factor = static_cast<T>(1) / static_cast<T>(nvirtual);
                
                //post("number of virtual %i", nvirtual);
                for(size_t i = 0; i < nvirtual; i++)
                {
                    const T angle = T(i) / T(nvirtual) * HOA_2PI;
                    //post("virtual %i :  %f", (int)i , (float)(angle / HOA_2PI * 360.f));
                    if(angle < channels[0].getAzimuth())
                    {
                        const T portion = (HOA_2PI - channels[channels.size()-1].getAzimuth()) + channels[0].getAzimuth();
                        
                        const T factor1 = (1. - ((channels[0].getAzimuth() - angle) / portion)) * factor;
                        encoder.setAzimuth(angle);
                        encoder.process(&factor1, vector_harmonics);
                        vector_harmonics[0] = factor1 * 0.5;
                        Signal<T>::add(harmonics, vector_harmonics, m_matrix + channels[0].getIndex() * harmonics);
                        
                        const T factor2 = ((channels[0].getAzimuth() - angle) / portion) * factor;
                        encoder.process(&factor2, vector_harmonics);
                        vector_harmonics[0] = factor2 * 0.5;
                        Signal<T>::add(harmonics, vector_harmonics, m_matrix + channels[channels.size() - 1].getIndex() * harmonics);
                        
                        //post("portion : %f", (float)portion / HOA_2PI * 360.f);
                        //post("channel %i (%f) : %f", (int)channels[channels.size()-1].getIndex(),
                        //     (float)(channels[channels.size()-1].getAzimuth() / HOA_2PI * 360.f),
                        //     (float)(factor2 / factor) * 100.f);
                        //post("channel %i (%f) : %f", (int)channels[0].getIndex(),
                        //     (float)(channels[0].getAzimuth() / HOA_2PI * 360.f),
                        //     (float)(factor1 / factor) * 100.f);
                    }
                    else if(angle >= channels[channels.size() - 1].getAzimuth())
                    {
                        const T portion = (HOA_2PI - channels[channels.size()-1].getAzimuth()) + channels[0].getAzimuth();
                        
                        const T factor1 = (1. - ((angle - channels[channels.size()-1].getAzimuth()) / portion)) * factor;
                        encoder.setAzimuth(angle);
                        encoder.process(&factor1, vector_harmonics);
                        vector_harmonics[0] = factor1 * 0.5;
                        Signal<T>::add(harmonics, vector_harmonics, m_matrix + channels[channels.size()-1].getIndex() * harmonics);
                        
                        const T factor2 = ((angle - channels[channels.size()-1].getAzimuth()) / portion) * factor;
                        encoder.process(&factor2, vector_harmonics);
                        vector_harmonics[0] = factor2 * 0.5;
                        Signal<T>::add(harmonics, vector_harmonics, m_matrix + channels[0].getIndex() * harmonics);
                        
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
                            if(angle < channels[j].getAzimuth() && angle >= channels[j-1].getAzimuth())
                            {
                                const T portion = (channels[j].getAzimuth() - channels[j-1].getAzimuth());
                                
                                const T factor1 = (1. - ((channels[j].getAzimuth() - angle) / portion)) * factor;
                                encoder.setAzimuth(angle);
                                encoder.process(&factor1, vector_harmonics);
                                vector_harmonics[0] = factor1 * 0.5;
                                Signal<T>::add(harmonics, vector_harmonics, m_matrix + channels[j].getIndex() * harmonics);
                                
                                const T factor2 = ((channels[j].getAzimuth() - angle) / portion) * factor;
                                encoder.process(&factor2, vector_harmonics);
                                vector_harmonics[0] = factor2 * 0.5;
                                Signal<T>::add(harmonics, vector_harmonics, m_matrix + channels[j-1].getIndex() * harmonics);
                                
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
    
    // ================================================================================ //
    // DECODER 2D BINAURAL //
    // ================================================================================ //
    
    //! @brief The ambisonic binaural decoder.
    //! @details The binaural decoder should be used to decode an ambisonic sound field for headphones.
    template <typename T, typename HrirType>
    class DecoderBinaural<Hoa2d, T, HrirType>
    : public Decoder<Hoa2d, T>
    {
    public:
        
        static_assert(HrirType::dimension == Hoa2d, "not a valid 2D response");
        
        using hrir_t = Hrir<Hoa2d, HrirType>;
        
        //! @brief Constructor.
        //! @param order The order
        DecoderBinaural(const size_t order)
        : Decoder<Hoa2d, T>(order, 2)
        {
            Decoder<Hoa2d, T>::setPlanewaveAzimuth(0, static_cast<T>(HOA_PI2*3.));
            Decoder<Hoa2d, T>::setPlanewaveAzimuth(1, static_cast<T>(HOA_PI2));
            setCropSize(0ul);
        }
        
        //! @brief Destructor.
        ~DecoderBinaural() noexcept
        {
            clear();
        }
        
        //! @brief This method retrieves the mode of the decoder.
        //! @return The mode of the decoder.
        inline typename Decoder<Hoa2d, T>::Mode getMode() const noexcept override
        {
            return Decoder<Hoa2d, T>::BinauralMode;
        }
        
        //! @brief This method sets the crop size of the responses.
        //! @param size The crop size.
        inline void setCropSize(const size_t size) noexcept
        {
            const auto num_rows = Hrir<Hoa2d, T>::getNumberOfRows();
            m_crop_size = (size == 0ul || size > num_rows) ? num_rows : size;
        }
        
        //! @brief This method gets the crop size of the responses.
        //! @return The crop size.
        inline size_t getCropSize() const noexcept
        {
            const auto num_rows = Hrir<Hoa2d, T>::getNumberOfRows();
            return (m_crop_size == num_rows) ? 0ul : m_crop_size;
        }
        
        //! @brief This method computes the decoding matrix.
        //! @param vectorsize The vector size for binaural decoding.
        void prepare(const size_t vectorsize = 64) override
        {
            clear();
            m_vector_size  = vectorsize;
            m_input  = Signal<T>::alloc(hrir_t::getNumberOfColumns() * m_vector_size);
            m_result = Signal<T>::alloc(hrir_t::getNumberOfRows() * m_vector_size);
            m_left   = Signal<T>::alloc(hrir_t::getNumberOfRows() + m_vector_size);
            m_right  = Signal<T>::alloc(hrir_t::getNumberOfRows() + m_vector_size);
        }
        
    public:
        
        //! @brief This method performs the binaural decoding and the convolution.
        inline void processBlock(const T** inputs, T** outputs) noexcept
        {
            T* input = m_input;
            for(size_t i = 0; i < Hrir<Hoa2d, T>::getNumberOfColumns() && i < Decoder<Hoa2d, T>::getNumberOfHarmonics(); i++)
            {
                Signal<T>::copy(m_vector_size, inputs[i], input);
                input += m_vector_size;
            }
            processChannel(m_input, hrir_t::template getLeftMatrix<T>(), m_left, outputs[0]);
            processChannel(m_input, hrir_t::template getRightMatrix<T>(), m_right, outputs[1]);
        }
        
        inline void process(const T* inputs, T* outputs) noexcept override
        {
            hoa_unused(inputs);
            hoa_unused(outputs);
            
            assert(true && "use the processBlock method instead");
        }
        
    private:
        
        inline void processChannel(const T* harmonics, const T* response, T* vector, T* output) noexcept
        {
            const size_t l = hrir_t::getNumberOfColumns();   // Harmonics size aka 11
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
        
        void clear()
        {
            m_input  = Signal<T>::free(m_input);
            m_result = Signal<T>::free(m_result);
            m_left   = Signal<T>::free(m_left);
            m_right  = Signal<T>::free(m_right);
        }
        
    private:
        size_t m_vector_size = 0ul;
        size_t m_crop_size = 0ul;
        T* m_input = nullptr;
        T* m_result = nullptr;
        T* m_left = nullptr;
        T* m_right = nullptr;
    };
    
    // ================================================================================ //
    // DECODER 3D //
    // ================================================================================ //
    
    //! @brief The decoder 3D decodes a sound field in the harmonics domain through the planewaves domain.
    //! @details The decoder should be used to decode a set of harmonics
    //! to a set of planewaves for loudspeakers.
    //! There are two types of decoders.
    //! - Regular for a perfect circle or sphere of loudspeakers.
    //! - Binaural for headphone restitution.
    template <typename T>
    class Decoder<Hoa3d, T>
    : public ProcessorHarmonics<Hoa3d, T>
    , public ProcessorPlanewaves<Hoa3d, T>
    {
    public:
        
        enum Mode { RegularMode = 0, BinauralMode = 2 };
        
        //! @brief Constructor.
        //! @param order The order
        //! @param channels The number of channels.
        Decoder(const size_t order, const size_t channels)
        : ProcessorHarmonics<Hoa3d, T>(order)
        , ProcessorPlanewaves<Hoa3d, T>(channels)
        {}
        
        //! @brief Destructor.
        virtual ~Decoder() {}
        
        //! @brief This method retrieves the mode of the decoder.
        //! @return The mode of the decoder.
        inline virtual Mode getMode() const noexcept { return RegularMode; }
        
        //! @brief This method performs the decoding.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! @param inputs  The input array that contains the samples of the harmonics.
        //! @param outputs The output array that contains samples destinated to the channels.
        virtual void process(const T* inputs, T* outputs) noexcept = 0;
        
        //! @brief This method computes the decoding matrix.
        //! @details You should use this method after changing the position of the loudspeakers.
        //! @param vectorsize The vector size for binaural decoding.
        virtual void prepare(const size_t vectorsize = 64) = 0;
    };
    
    // ================================================================================ //
    // DECODER 3D BINAURAL //
    // ================================================================================ //
    
    //! @brief The ambisonic binaural decoder.
    //! @details The binaural decoder should be used to decode an ambisonic sound field for headphones.
    template <typename T, typename HrirType>
    class DecoderBinaural<Hoa3d, T, HrirType>
    : public Decoder<Hoa3d, T>
    {
    public:
        
        static_assert(HrirType::dimension == Hoa3d, "not a valid 3D response");
        
        using hrir_t = Hrir<Hoa3d, HrirType>;
        
    public:
        
        //! @brief Constructor.
        //! @param order The order.
        DecoderBinaural(const size_t order)
        : Decoder<Hoa3d, T>(order, 2)
        {
            Decoder<Hoa3d, T>::setPlanewaveAzimuth(0, static_cast<T>(HOA_PI2*3.));
            Decoder<Hoa3d, T>::setPlanewaveAzimuth(1, static_cast<T>(HOA_PI2));
        }
        
        //! @brief This method retrieves the mode of the decoder.
        //! @return The mode of the decoder.
        typename Decoder<Hoa3d, T>::Mode getMode() const noexcept override
        {
            return Decoder<Hoa3d, T>::BinauralMode;
        }
        
        //! @brief Destructor.
        ~DecoderBinaural() = default;
        
        //! @brief This method sets the crop size of the responses.
        //! @param size The crop size.
        void setCropSize(const size_t size) noexcept
        {
            const auto num_rows = hrir_t::getNumberOfRows();
            m_crop_size = (size == 0ul || size > num_rows) ? num_rows : size;
        }
        
        //! @brief This method gets the crop size of the responses.
        //! @return The crop size.
        size_t getCropSize() const noexcept
        {
            const auto num_rows = hrir_t::getNumberOfRows();
            return (m_crop_size == num_rows) ? 0ul : m_crop_size;
        }
        
        //! @brief This method computes the decoding matrix.
        //! @param vectorsize The vector size for binaural decoding.
        void prepare(const size_t vectorsize = 64) override
        {
            m_vector_size = vectorsize;
            
            const auto response_size = hrir_t::getNumberOfRows();
            const auto number_of_harmonics = hrir_t::getNumberOfColumns();
            
            // the setZero method resize and set matrice and vector coefficients to 0.
            
            m_input.setZero(number_of_harmonics, m_vector_size);
            m_result.setZero(response_size, m_vector_size);
            
            const auto side_vecsize = response_size + m_vector_size;
            m_left.setZero(side_vecsize);
            m_right.setZero(side_vecsize);
        }
        
    public:
        
        //! @brief This method performs the binaural decoding and the convolution.
        inline void processBlock(const T** inputs, T** outputs) noexcept
        {
            const auto ins = std::min(Decoder<Hoa3d, T>::getNumberOfHarmonics(),
                                      hrir_t::getNumberOfColumns());
            
            for(auto i = 0; i < ins; ++i)
            {
                m_input.row(i).noalias() = vector_t::Map(inputs[i], m_vector_size);
            }
            
            processChannel(m_input, m_left_matrix, m_left, outputs[0]);
            processChannel(m_input, m_right_matrix, m_right, outputs[1]);
        }
        
        inline void process(const T* inputs, T* outputs) noexcept override
        {
            hoa_unused(inputs);
            hoa_unused(outputs);
            
            assert(true && "use the processBlock method instead");
        }
        
    private:
        
        using vector_t = Eigen::RowVectorX<T>;
        using matrix_t = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using hrir_matrix_t = Eigen::Map<const Eigen::Matrix<T, hrir_t::getNumberOfRows(), hrir_t::getNumberOfColumns(), Eigen::RowMajor>>;
        
        void processChannel(matrix_t const& harmonics_matrix,
                            hrir_matrix_t const& response_matrix,
                            vector_t& buffer,
                            T* output)
        {
            const size_t m = m_crop_size;    // crop <= response size
            const size_t n = m_vector_size;

            m_result.topRows(m).noalias() = response_matrix.topRows(m) * harmonics_matrix;
            
            for(size_t i = 0; i < n; ++i)
            {
                buffer.segment(i, m).noalias() += m_result.col(i).head(m);
            }
            
            vector_t::Map(output, n) = buffer.head(n);
            buffer.head(m) = buffer.segment(n, m);
            buffer.segment(m, n).setZero();
        }
        
    private:
        
        size_t  m_vector_size = 0ul;
        size_t  m_crop_size = hrir_t::getNumberOfRows();
        matrix_t m_input = {};
        matrix_t m_result = {};
        vector_t m_left = {};
        vector_t m_right = {};
        
        const hrir_matrix_t m_left_matrix = {
            hrir_t::template getLeftMatrix<T>(),
            hrir_t::getNumberOfRows(),
            hrir_t::getNumberOfColumns()
        };
        
        const hrir_matrix_t m_right_matrix = {
            hrir_t::template getRightMatrix<T>(),
            hrir_t::getNumberOfRows(),
            hrir_t::getNumberOfColumns()
        };
    };
    
    // syntactic sugar
    
    template<typename T>
    using decoder_binaural_2d_t = DecoderBinaural<Hoa2d, T, hrir::Sadie_D2_2D>;
    
    template<typename T>
    using decoder_binaural_3d_t = DecoderBinaural<Hoa3d, T, hrir::Sadie_D2_3D>;
}
