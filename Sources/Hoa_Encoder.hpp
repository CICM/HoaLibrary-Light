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

#ifndef DEF_HOA_ENCODER_LIGHT
#define DEF_HOA_ENCODER_LIGHT

#include "Hoa_Processor.hpp"

namespace hoa
{

    //! The encoder class generates the harmonics for one or several signal according to an azimuth, an elevation and a radius.
    /** The encoder should be used to encode a signal in the harmonics domain depending on an order of decomposition. It allows to control the azimuth, the elevation and potentially the radius of the signal.
     */
    template <Dimension D, typename T> class Encoder : public ProcessorHarmonics<D, T>
    {
    public:

        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes harmonics coefficients for the encoding. The order must be at least 1.
         @param     order	The order.
         */
        Encoder(const size_t order) hoa_noexcept;

        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
		virtual ~Encoder() hoa_noexcept = 0;
    };
    
    //! The basic encoder class generates the harmonics for one signal according to an azimuth and an elevation.
    /** The basic encoder should be used to encode a signal in the harmonics domain depending on an order of decomposition. It allows to control the azimuth and the elevation of the signal.
     */
    template <Dimension D, typename T> class EncoderBasic : public Encoder<D, T>
    {
    public:
        //! The basic constructor.
        /**	The basic constructor allocates and initialize the member values to computes harmonics coefficients for the encoding. The order must be at least 1.
         @param     order	The order.
         */
        EncoderBasic(const size_t order) hoa_noexcept;
        
        //! The basic destructor.
        /**	The basic destructor free the memory.
         */
        virtual ~EncoderBasic() hoa_noexcept = 0;
        
        //! Mute or unmute the process.
        /**	This method mutes or unmutes the process.
         @param     muted	The mute state.
         */
        virtual void setMute(const bool muted) hoa_noexcept;
        
        //! Get the mute or unmute state of the process.
        /**	This method gets mute state of the process.
         @return    The mute state of the process.
         */
        virtual bool getMute() const hoa_noexcept;
        
        //! Set the azimuth.
        /**	This method  sets the azimuth \f$\theta\f$ in radian and you should prefer to use it between \f$0\f$ and \f$2\pi\f$ to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The \f$0\f$ radian is \f$\frac{\pi}{2}\f$ phase shifted relative to a mathematical representation of a circle, then the \f$0\f$ radian is at the "front" of the soundfield.
         @param     azimuth	The azimuth.
         */
        virtual void setAzimuth(const T azimuth) hoa_noexcept;
        
        //! Get the azimuth
        /** The method returns the azimuth \f$\theta\f$ between \f$0\f$ and \f$2\pi\f$.
         @return     The azimuth.
         */
        virtual T getAzimuth() const hoa_noexcept;
        
        //! Set the elevation.
        /**	This method  sets the elevation \f$\varphi\f$ in radian and you should prefer to use it between \f$-\pi\f$ and \f$\pi\f$ to avoid recursive wrapping of the value. The direction of rotation is from bottom to the top. The \f$0\f$ radian is centered at the "front" of the soundfield, then \f$\frac{\pi}{2}\f$ is at the top, \f$-\frac{\pi}{2}\f$ is at the bottom and \f$\pi\f$ is behind. Note that if the angle of elevation is between \f$\frac{\pi}{2}\f$ and \f$\frac{3\pi}{2}\f$, the azimuth is reversed.
         @param     elevation The elevation.
         */
        virtual void setElevation(const T elevation) hoa_noexcept;
        
        //!	Get the elevation.
        /** The method returns the elevation \f$\varphi\f$ between \f$-\pi\f$ and \f$\pi\f$.
         @return     The elevation.
         */
        virtual T getElevation()  const hoa_noexcept;
        
        //! This method performs the encoding.
        /**	You should use this method for not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         \f[Y_{l,m}(\theta, \varphi) = k_{l, m} P_{l, \left|m\right|}(\cos{(\varphi)}) e^{+im\theta} \f]
         with \f$e^{+im\theta}\f$ the azimuth part of the equation with \f$i\f$ the imaginary, \f$P_{l, \left|m\right|}(\cos{(\varphi)})\f$ the elevation part of the equation with \f$P_{l, \left|m\right|}(x)\f$ the associated Legendre polynomials, \f$k_{l, m}\f$ the normalization, \f$l\f$ the degree, \f$m\f$ the order, \f$\theta\f$ the azimuth in radian and \f$\varphi\f$ the elevation in radian.\n
         
         The azimuth part in the imaginary form \f$e^{+im\theta}\f$ can be expressed with the real form :\n
         if \f$m \geq 0\f$ then
         \f[e^{+im\theta} = \cos{(\left|m\right|\theta)}\f]
         else
         \f[e^{+im\theta} = sin{(\left|m\right|\theta)}\f]
         The elevation part \f$P_{l, \left|m\right|}(x)\f$ can be expressed with recursives formulas :
         \f[P_{l+1,l+1}(x) = -(2l+1)\sqrt{(1-x^2)}P_{(l,l)}(x) \f]
         \f[P_{l+1,l}(x) = x(2l+1)P_{(l,l)}(x) \f]
         \f[P_{l+1,m}(x) = \frac{x(2l+1)P_{(l,m)}(x) - (l+m)P_{(l-1,m)}(x)}{l-m+1} \f]
         and with \f[P_{0, 0}(x) = 1\f]
         The normalization part \f$k_{l, m}\f$ is equivalent to :\n
         if \f$m = 0\f$ then
         \f[k_{l, m} = 1\f]
         else
         \f[k_{l, m} = \sqrt{\frac{(l - \left|m\right|)!}{(l + \left|m\right|)!}}\sqrt{2} \f]
         
         @param     input    The pointer to the input sample.
         @param     outputs  The outputs array.
         */
        virtual void process(const T* input, T* outputs) hoa_noexcept hoa_override;
        
        //! This method performs the encoding but add the result to the outputs.
        /**	You should use this method for not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @see process
         @param     input    The pointer to the input sample.
         @param     outputs  The outputs array.
         */
        virtual void processAdd(const T* input, T* outputs) hoa_noexcept;
    };
    
    //! The dc encoder class generates the harmonics for one signal according to an azimuth, an elevation and a radius.
    /** The dc encoder should be used to encode a signal in the harmonics domain depending on an order of decomposition. It allows to control the azimuth, the elevation and the radius of the signal. The distance compensation is performed with the simulation of fractional orders when the signal is inside the ambisonic circle or sphere and with gain attenuation when the signal is outside the ambisonics circle or sphere.
     */
    template <Dimension D, typename T> class EncoderDC : public Encoder<D, T>
    {
    public:
        
        //! The dc constructor.
        /**	The dc constructor allocates and initialize the member values to computes harmonics coefficients for the encoding. The order must be at least 1.
         @param     order	The order.
         */
        EncoderDC(const size_t order) hoa_noexcept;
        
        //! The dc destructor.
        /**	The dc destructor free the memory.
         */
        virtual ~EncoderDC() hoa_noexcept = 0;
        
        //! Mute or unmute the process.
        /**	This method mutes or unmutes the process.
         @param     muted	The mute state.
         */
        virtual void setMute(const bool muted) hoa_noexcept;
        
        //! Get the mute or unmute state of the process.
        /**	This method gets mute state of the process.
         @return    The mute state of the process.
         */
        virtual bool getMute() const hoa_noexcept;
        
        //! Set the azimuth.
        /**	This method  sets the azimuth \f$\theta\f$ in radian and you should prefer to use it between \f$0\f$ and \f$2\pi\f$ to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The \f$0\f$ radian is \f$\frac{\pi}{2}\f$ phase shifted relative to a mathematical representation of a circle, then the \f$0\f$ radian is at the "front" of the soundfield.
         @param     azimuth	The azimuth.
         */
        virtual void setAzimuth(const T azimuth) hoa_noexcept;
        
        //! Get the azimuth.
        /** The method returns the azimuth \f$\theta\f$ between \f$0\f$ and \f$2\pi\f$.
         @return     The azimuth.
         */
        virtual T getAzimuth() const hoa_noexcept;
        
        //! Set the elevation.
        /**	This method  sets the elevation \f$\varphi\f$ in radian and you should prefer to use it between \f$-\pi\f$ and \f$\pi\f$ to avoid recursive wrapping of the value. The direction of rotation is from bottom to the top. The \f$0\f$ radian is centered at the "front" of the soundfield, then \f$\frac{\pi}{2}\f$ is at the top, \f$-\frac{\pi}{2}\f$ is at the bottom and \f$\pi\f$ is behind. Note that if the angle of elevation is between \f$\frac{\pi}{2}\f$ and \f$\frac{3\pi}{2}\f$, the azimuth is reversed.
         @param     elevation The elevation.
         */
        virtual void setElevation(const T elevation) hoa_noexcept;
        
        //!	Get the elevation.
        /** The method returns the elevation \f$\varphi\f$ between \f$-\pi\f$ and \f$\pi\f$.
         @return     The elevation.
         */
        virtual T getElevation()  const hoa_noexcept;
        
        //! Set the radius.
        /**	This method  sets the radius \f$\rho\f$ between \f$0\f$ and \f$+\infty\f$. \f$0\f$ is the center of the soundfield, \f$1\f$ is the radius of the ambisonics circle or sphere, beyond this limit the gain decreases and before the sound field is widened.
         @param     radius The radius.
         */
        virtual void setRadius(const T radius) hoa_noexcept;
        
        //! Get the radius.
        /** The method returns the radius \f$\rho\f$ between \f$0\f$ and \f$+\infty\f$.
         @return     The radius.
         */
        virtual T getRadius() const hoa_noexcept;
        
        //! This method performs the encoding.
        /**	You should use this method for not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         \f[Y^{dc}_{l,m}(\theta, \varphi, \rho) = (\frac{1}{\max{(\rho, 1)}})Y^{widened}_{l,m}(\rho) \leftarrow Y_{l,m}(\theta, \varphi) \f]
         with \f$Y_{l,m}\f$ the basic encoding, \f$Y^{widened}_{l,m}\f$ the widening operation, \f$l\f$ the degree, \f$m\f$ the order, \f$\theta\f$ the azimuth in radian, \f$\varphi\f$ the elevation in radian and \f$\rho\f$ the radius.\n
         @see Basic
         @see Wider
         @param     input    The pointer to the input sample.
         @param     outputs  The outputs array.
         */
        virtual void process(const T* input, T* outputs) hoa_noexcept;
        
        //! This method performs the encoding but add the result to the outputs.
        /**	You should use this method for not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @see process
         @param     input    The pointer to the input sample.
         @param     outputs  The outputs array.
         */
        inline void processAdd(const T* input, T* outputs) hoa_noexcept;
        
    };
    
    //! The multi encoder class generates the harmonics for several signals according to an azimuth, an elevation and a radius for each one.
    /** The multi encoder should be used to encode several signals in the harmonics domain depending on an order of decomposition. It allows to control the azimuth, the elevation and the radius of each signal. The class uses a set of dc encoders.
     */
    template <Dimension D, typename T> class EncoderMulti : public Encoder<D, T>
    {
    public:
        
        //! The multi encoder constructor.
        /**	The multi encoder constructor allocates and initialize the member values and classes depending on a order of decomposition and the number of sources. The order and the number of sources must be at least 1.
         
         @param     order            The order.
         @param     numberOfSources	The number of sources.
         */
        EncoderMulti(const size_t order, size_t numberOfSources) hoa_noexcept;
        
        //! The multi encoder destructor.
        /**	The multi encoder destructor free the memory and deallocate the member classes.
         */
        virtual ~EncoderMulti() hoa_noexcept = 0;
        
        //! This method retrieve the number of sources.
        /** Retrieve the number of sources.
         @return The number of sources.
         */
        virtual size_t getNumberOfSources() const hoa_noexcept;
        
        //! Set the azimuth of a signal.
        /**	This method  sets the azimuth \f$\theta_{index}\f$ of a signal in radian and you should prefer to use it between \f$0\f$ and \f$2\pi\f$ to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The \f$0\f$ radian is \f$\frac{\pi}{2}\f$ phase shifted relative to a mathematical representation of a circle, then the \f$0\f$ radian is at the "front" of the soundfield.
         @param index   The index of the signal.
         @param azimuth	The azimuth.
         */
        virtual void setAzimuth(const size_t index, const T azimuth) hoa_noexcept;
        
        //! Set the elevation  of a signal.
        /**	This method  sets the elevation \f$\varphi_{index}\f$  of a signal in radian and you should prefer to use it between \f$-\pi\f$ and \f$\pi\f$ to avoid recursive wrapping of the value. The direction of rotation is from bottom to the top. The \f$0\f$ radian is centered at the "front" of the soundfield, then \f$\frac{\pi}{2}\f$ is at the top, \f$-\frac{\pi}{2}\f$ is at the bottom and \f$\pi\f$ is behind. Note that if the angle of elevation is between \f$\frac{\pi}{2}\f$ and \f$\frac{3\pi}{2}\f$, the azimuth is reversed.
         @param index       The index of the signal.
         @param elevation   The elevation.
         */
        virtual void setElevation(const size_t index, const T elevation) hoa_noexcept;
        
        //! Set the radius.
        /**	This method  sets the radius \f$\rho_{index}\f$ between \f$0\f$ and \f$+\infty\f$. \f$0\f$ is the center of the soundfield, \f$1\f$ is the radius of the ambisonic circle or sphere, beyond this limit the gain decreases and before the sound field is widened.
         @param index   The index of the signal.
         @param radius  The radius.
         */
        virtual void setRadius(const size_t index, const T radius) hoa_noexcept;
        
        //! This method mute or unmute a signal.
        /**	Mute or unmute a signal with a boolean value.
         @param     index	The index of the signal.
         @param     muted	The mute state.
         */
        virtual void setMute(const size_t index, const bool muted) hoa_noexcept;
        
        //! Get the azimuth of a signal.
        /** The method returns the azimuth \f$\theta_{index}\f$ between \f$0\f$ and \f$2\pi\f$.
         @param index The index of the signal.
         @return    The azimuth.
         */
        virtual T getAzimuth(const size_t index) const hoa_noexcept;
        
        //!	Get the elevation of a signal.
        /** The method returns the elevation \f$\varphi_{index}\f$ between \f$-\pi\f$ and \f$\pi\f$.
         @param index The index of the signal.
         @return     The elevation.
         */
        virtual T getElevation(const size_t index) const hoa_noexcept;
        
        //! Get the radius of a signal.
        /** The method returns the radius \f$\rho_{index}\f$ between \f$0\f$ and \f$+\infty\f$.
         @param index The index of the signal.
         @return     The radius.
         */
        virtual T getRadius(const size_t index) const hoa_noexcept;
        
        //! Get the mute or unmute state of a signal.
        /**	This method gets mute state of a signal.
         @param index The index of the signal.
         @return    The mute state of the signal.
         */
        virtual bool getMute(const size_t index) const hoa_noexcept;
        
        
        //! This method performs the encoding with distance compensation.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The input array contains the samples of the sources and the minimum size should be the number of sources. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         \f[Y^{multi}_{l,m}(\theta_0^n, \varphi_0^n, \rho_0^n) = \sum_{i=0}^n Y^{dc}_{l,m}(\theta_i, \varphi_i, \rho_i) \f]
         @param     input  The input array.
         @param     outputs The outputs array.
         */
        virtual void process(const T* input, T* outputs) hoa_noexcept;
    };


    
    
    
    
    
    
    
    
    
    
    //! The encoder class generates the harmonics for one or several signal according to an azimuth, an elevation and a radius.
    /** The encoder should be used to encode a signal in the harmonics domain depending on an order of decomposition. It allows to control the azimuth, the elevation and potentially the radius of the signal.
     */
    template <typename T> class Encoder<Hoa2d, T> : public ProcessorHarmonics<Hoa2d, T>
    {
    public:

        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes harmonics coefficients for the encoding. The order must be at least 1.
         @param     order	The order.
         */
        Encoder(const size_t order) hoa_noexcept : ProcessorHarmonics<Hoa2d, T>(order) {}

        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
        virtual ~Encoder() hoa_noexcept hoa_default_f
        
        //! This method performs the encoding.
        /**	You should use this method for not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The pointer to the input sample of the inputs samples.
         @param     outputs  The outputs array.
         */
        virtual void process(const T* input, T* outputs) hoa_noexcept = 0;
    };

    //! The basic encoder class generates the harmonics for one signal according to an azimuth and an elevation.
    /** The basic encoder should be used to encode a signal in the harmonics domain depending on an order of decomposition. It allows to control the azimuth and the elevation of the signal.
     */
    template <typename T> class EncoderBasic<Hoa2d, T> : public Encoder<Hoa2d, T>
    {
    private:
        T    m_azimuth;
        T    m_cosx;
        T    m_sinx;
        bool m_muted;
    public:

        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes circular harmonics coefficients for the encoding. The order must be at least 1.
            @param     order	The order.
         */
        EncoderBasic(const size_t order) hoa_noexcept: Encoder<Hoa2d, T>(order)
        {
            setMute(false);
            setAzimuth(0.);
        }

        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
        ~EncoderBasic() hoa_noexcept hoa_default_f

        //! This method set the azimuth.
        /**	The azimuth in radian and you should prefer to use it between 0 and 2π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield.
            @param     azimuth	The azimuth.
         */
        inline void setAzimuth(const T azimuth) hoa_noexcept
        {
            m_azimuth = azimuth;
            m_cosx    = std::cos(m_azimuth);
            m_sinx    = std::sin(m_azimuth);
        }

        //! Get the azimuth.
        /** The method returns the azimuth between 0 and 2π.
            @return     The azimuth.
         */
        inline T getAzimuth() const hoa_noexcept { return Math<T>::wrap_twopi(m_azimuth);}

        //! This method mute or unmute.
        /**	Mute or unmute.
         @param     muted	The mute state.
         */
        inline void setMute(const bool muted) hoa_noexcept { m_muted = muted; }

        //! This method retrieve the mute or unmute state of a source.
        /**	Get the Mute state of a source.
         @return    The mute state of the source.
         */
        inline bool getMute() const hoa_noexcept { return m_muted; }

        //! This method performs the encoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
        // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
         @param     input	The input sample.
         @param     outputs The output array.
         */
        inline void process(const T* input, T* outputs) hoa_noexcept hoa_override
        {
            if(!m_muted)
            {
                T cos_x = m_cosx;
                T sin_x = m_sinx;
                T tcos_x = cos_x;
                (*outputs++)    = (*input);                         // Hamonic [0,0]
                (*outputs++)    = (*input) * sin_x;                 // Hamonic [1,-1]
                (*outputs++)    = (*input) * cos_x;                 // Hamonic [1,1]
                for(size_t i = 2; i <= ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder(); i++)
                {
                    cos_x   = tcos_x * m_cosx - sin_x * m_sinx;
                    sin_x   = tcos_x * m_sinx + sin_x * m_cosx;
                    tcos_x  = cos_x;
                    (*outputs++)    = (*input) * sin_x;            // Hamonic [l,-l]
                    (*outputs++)    = (*input) * cos_x;            // Hamonic [l,l]
                }
            }
            else
            {
                for(size_t i = 0; i < ProcessorHarmonics<Hoa2d, T>::getNumberOfHarmonics(); i++)
                {
                    (*outputs++) = 0.;
                }
            }
        }

        //! This method performs the encoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
         // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
         @param     input	The input sample.
         @param     outputs The output array.
         */
        inline void processAdd(const T* input, T* outputs) hoa_noexcept
        {
            if(!m_muted)
            {
                T cos_x = m_cosx;
                T sin_x = m_sinx;
                T tcos_x = cos_x;
                (*outputs++)    += (*input);                         // Hamonic [0,0]
                (*outputs++)    += (*input) * sin_x;                 // Hamonic [1,-1]
                (*outputs++)    += (*input) * cos_x;                 // Hamonic [1,1]
                for(size_t i = 2; i <= ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder(); i++)
                {
                    cos_x   = tcos_x * m_cosx - sin_x * m_sinx;
                    sin_x   = tcos_x * m_sinx + sin_x * m_cosx;
                    tcos_x  = cos_x;

                    (*outputs++)    += (*input) * sin_x;            // Hamonic [i,-i]
                    (*outputs++)    += (*input) * cos_x;            // Hamonic [i,i]
                }
            }
        }
    };

    //! The dc encoder class generates the harmonics for one signal according to an azimuth, an elevation and a radius.
    /** The dc encoder should be used to encode a signal in the harmonics domain depending on an order of decomposition. It allows to control the azimuth, the elevation and the radius of the signal. The distance compensation is performed with the simulation of fractional orders when the signal is inside the ambisonic circle or sphere and with gain attenuation when the signal is outside the ambisonics circle or sphere.
     */
    template <typename T> class EncoderDC<Hoa2d, T> : public Encoder<Hoa2d, T>
    {
    private:
        T   m_azimuth;
        T   m_cosx;
        T   m_sinx;
        T   m_factor;
        T   m_gain;
        T   m_radius;
        T   m_distance;
        bool m_muted;
    public:

        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes circular harmonics coefficients for the encoding. The order must be at least 1.
         @param     order	The order.
         */
        EncoderDC(const size_t order) hoa_noexcept: Encoder<Hoa2d, T>(order)
        {
            setAzimuth(0.);
            setRadius(1.);
            setMute(false);
        }

        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
        ~EncoderDC() hoa_noexcept hoa_default_f

        //! This method set the azimuth.
        /**	The azimuth in radian and you should prefer to use it between 0 and 2π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield.
         @param     azimuth	The azimuth.
         */
        inline void setAzimuth(const T azimuth) hoa_noexcept
        {
            m_azimuth = azimuth;
            m_cosx    = std::cos(m_azimuth);
            m_sinx    = std::sin(m_azimuth);
        }

        //! Get the azimuth.
        /** The method returns the azimuth between 0 and 2π.
         @return     The azimuth.
         */
        inline T getAzimuth() const hoa_noexcept { return Math<T>::wrap_twopi(m_azimuth); }

        //! This method set the radius.
        /**	The radius is between 0 and infinity. At 0, the source is in the center of the ambisonic circle and at 1, the source is at the limit of the ambisonic circle. Over 1, the source get away the ambisonic circle.
         @param     radius   The radius.
         @see       setAzimuth()
         */
        inline void setRadius(const T radius) hoa_noexcept
        {
            m_radius = std::max(radius, (T)0.);
            if(m_radius < 1.)
            {
                m_factor    = T((1. - m_radius) * HOA_PI);
                m_gain      = T((sin(m_factor - HOA_PI2) + 1.) * 0.5);
                m_distance  = 1.;
            }
            else
            {
                m_factor    = 0;
                m_gain      = 0;
                m_distance  = T(1. / radius);
            }
        }

        //! Get the azimuth.
        /** The method returns the azimuth between 0 and 2π.
         @return     The azimuth.
         */
        inline T getRadius() const hoa_noexcept { return m_radius; }

        //! This method mute or unmute.
        /**	Mute or unmute.
         @param     muted	The mute state.
         */
        inline void setMute(const bool muted) hoa_noexcept { m_muted = muted; }

        //! This method retrieve the mute or unmute state of a source.
        /**	Get the Mute state of a source.
         @return    The mute state of the source.
         */
        inline bool getMute() const hoa_noexcept { return m_muted; }

        //! This method performs the encoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
         // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
         @param     input	The input sample.
         @param     outputs The output array.
         */
        inline void process(const T* input, T* outputs) hoa_noexcept hoa_override
        {
            if(!m_muted)
            {
                T cos_x = m_cosx;
                T sin_x = m_sinx;
                T tcos_x = cos_x;
                const T gain1   = (m_gain * ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder());
                const T factor1 = (cos(Math<T>::clip(m_factor, 0., T(HOA_PI))) + 1.) * T(0.5) * T((gain1 - m_gain) + m_distance);

                (*outputs++) = (*input) * (gain1 + m_distance);            // Hamonic [0,0]
                (*outputs++) = (*input) * sin_x * factor1;                 // Hamonic [1,-1]
                (*outputs++) = (*input) * cos_x * factor1;                 // Hamonic [1,1]
                for(size_t i = 2; i <= ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder(); i++)
                {
                    const T gain    = (m_gain * (ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder() - i) + m_distance);
                    const T factor  = (cos(Math<T>::clip(m_factor * i, 0., T(HOA_PI))) + 1.) * T(0.5);

                    cos_x   = tcos_x * m_cosx - sin_x * m_sinx;
                    sin_x   = tcos_x * m_sinx + sin_x * m_cosx;
                    tcos_x  = cos_x;

                    (*outputs++)    = (*input) * sin_x * factor * gain;    // Hamonic [i,-i]
                    (*outputs++)    = (*input) * cos_x * factor * gain;    // Hamonic [i,i]
                }
            }
            else
            {
                for(size_t i = 0; i < ProcessorHarmonics<Hoa2d, T>::getNumberOfHarmonics(); i++)
                {
                    (*outputs++) = 0.;
                }
            }
        }

        //! This method performs the encoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
         // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
         @param     input	The input sample.
         @param     outputs The output array.
         */
        inline void processAdd(const T* input, T* outputs) hoa_noexcept
        {
            if(!m_muted)
            {
                T cos_x = m_cosx;
                T sin_x = m_sinx;
                T tcos_x = cos_x;
                const T gain1   = (m_gain * ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder());
                const T factor1 = (cos(Math<T>::clip(m_factor, 0., HOA_PI)) + 1.) * 0.5 * ((gain1 - m_gain) + m_distance);

                (*outputs++) += (*input) * (gain1 + m_distance);            // Hamonic [0,0]
                (*outputs++) += (*input) * sin_x * factor1;                 // Hamonic [1,-1]
                (*outputs++) += (*input) * cos_x * factor1;                 // Hamonic [1,1]
                for(size_t i = 2; i <= ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder(); i++)
                {
                    const T gain    = (m_gain * (ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder() - i) + m_distance);
                    const T factor  = (cos(Math<T>::clip(m_factor * i, 0., HOA_PI)) + 1.) * 0.5 ;

                    cos_x   = tcos_x * m_cosx - sin_x * m_sinx;
                    sin_x   = tcos_x * m_sinx + sin_x * m_cosx;
                    tcos_x  = cos_x;

                    (*outputs++)    += (*input) * sin_x * factor * gain;    // Hamonic [i,-i]
                    (*outputs++)    += (*input) * cos_x * factor * gain;    // Hamonic [i,i]
                }
            }
        }
    };


    //! The multi encoder class generates the harmonics for several signals according to an azimuth, an elevation and a radius for each one.
    /** The multi encoder should be used to encode several signals in the harmonics domain depending on an order of decomposition. It allows to control the azimuth, the elevation and the radius of each signal. The class uses a set of dc encoders.
     */
    template <typename T> class EncoderMulti<Hoa2d, T> : public Encoder<Hoa2d, T>
    {
    private:
        const size_t                     m_number_of_sources;
        std::vector< EncoderDC<Hoa2d, T>* >  m_encoders;
    public:

        //! The map constructor.
        /**	The map constructor allocates and initialize the member values and classes depending on a order of decomposition and the number of sources. The order and the number of sources must be at least 1.

         @param     order            The order.
         @param     numberOfSources	The number of sources.
         */
        EncoderMulti(const size_t order, size_t numberOfSources) hoa_noexcept : Encoder<Hoa2d, T>(order),
        m_number_of_sources(numberOfSources)
        {
            for(size_t i = 0; i < m_number_of_sources; i++)
            {
                m_encoders.push_back(new EncoderDC<Hoa2d, T>(order));
            }
        }

        //! The map destructor.
        /**	The map destructor free the memory and deallocate the member classes.
         */
        ~EncoderMulti() hoa_noexcept
        {
            for(size_t i = 0; i < m_number_of_sources; i++)
            {
                delete m_encoders[i];
            }
            m_encoders.clear();
        }

        //! This method retrieve the number of sources.
        /** Retrieve the number of sources.
         @return The number of sources.
         */
        inline size_t getNumberOfSources() const hoa_noexcept { return m_number_of_sources; }

        //! This method set the angle of azimuth of a source.
        /**	The angle of azimuth in radian and you should prefer to use it between 0 and 2 π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield. The index must be between 0 and the number of sources - 1.

         @param     index	The index of the source.
         @param     azimuth	The azimuth.
         @see       setRadius()
         */
        inline void setAzimuth(const size_t index, const T azimuth) hoa_noexcept { m_encoders[index]->setAzimuth(azimuth); }

        //! This method set the radius of a source.
        /**	The radius is between 0 and infinity. At 0, the source is in the center of the ambisonic circle and at 1, the source is at the limit of the ambisonic circle. Over 1, the source get away the ambisonic circle. The index must be between 0 and the number of sources - 1.
         @param     index	The index of the source.
         @param     radius   The radius.
         @see       setAzimuth()
         */
        inline void setRadius(const size_t index, const T radius) hoa_noexcept { m_encoders[index]->setRadius(radius); }

        //! This method mute or unmute a source.
        /**	Mute or unmute a source with a boolean value. The index must be between 0 and the number of sources - 1.
         @param     index	The index of the source.
         @param     muted	The mute state.
         */
        inline void setMute(const size_t index, const bool muted) hoa_noexcept { m_encoders[index]->setMute(muted); }

        //! This method retrieve the azimuth of a source.
        /** Retrieve the azimuth of a source.
         @param     index	The index of the source.
         @return The azimuth of the source if the source exists, otherwise the function generates an error.
         */
        inline T getAzimuth(const size_t index) const hoa_noexcept { return m_encoders[index]->getAzimuth(); }

        //! This method retrieve the radius of a source.
        /** Retrieve the radius of a source.

         @param     index	The index of the source.
         @return The radius of the source if the source exists, otherwise the function generates an error.
         */
        inline T getRadius(const size_t index) const hoa_noexcept { return m_encoders[index]->getRadius(); }

        //! This method retrieve the mute or unmute state of a source.
        /**	Get the Mute state of a source.

         @param     index	The index of the source.
         @return    The mute state of the source if the source exists, otherwise the function generates an error.
         @see       setMute()
         */
        inline bool getMute(const size_t index) const hoa_noexcept { return m_encoders[index]->getMute(); }


        //! This method performs the encoding with distance compensation.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The input array contains the samples of the sources and the minimum size should be the number of sources. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     input  The input array.
         @param     outputs The outputs array.
         */
        inline void process(const T* input, T* outputs) hoa_noexcept hoa_override
        {
            m_encoders[0]->process(input, outputs);
            for(size_t i = 1; i < m_number_of_sources; i++)
            {
                m_encoders[i]->processAdd(++input, outputs);
            }
        }
    };

    
    
    
    
    
    
    
    
    template <typename T> class Encoder<Hoa3d, T> : public ProcessorHarmonics<Hoa3d, T>
    {
    public:

        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes harmonics coefficients for the encoding. The order must be at least 1.
         @param     order	The order.
         */
        Encoder(const size_t order) hoa_noexcept : ProcessorHarmonics<Hoa3d, T>(order) {}

        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
        virtual ~Encoder() hoa_noexcept hoa_default_f

        //! This method performs the encoding.
        /**	You should use this method for not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The pointer to the input sample of the inputs samples.
         @param     outputs  The outputs array.
         */
        virtual void process(const T* input, T* outputs) hoa_noexcept = 0;
    };

    //! The basic encoder class generates the harmonics for one signal according to an azimuth and an elevation.
    /** The basic encoder should be used to encode a signal in the harmonics domain depending on an order of decomposition. It allows to control the azimuth and the elevation of the signal.
     */
    template <typename T> class EncoderBasic<Hoa3d, T> : public Encoder<Hoa3d, T>
    {
    private:
        T  m_azimuth;
        T  m_elevation;
        T  m_cos_phi;
        T  m_sin_phi;
        T  m_cos_theta;
        T  m_sqrt_rmin;
        T* m_normalization;
        bool m_muted;
    public:

        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes circular harmonics coefficients for the encoding. The order must be at least 1.
         @param     order	The order.
         */
        EncoderBasic(const size_t order) hoa_noexcept : Encoder<Hoa3d, T>(order)
        {
            m_normalization = Signal<T>::alloc(ProcessorHarmonics<Hoa3d, T>::getNumberOfHarmonics());
            for(size_t i = 0; i < ProcessorHarmonics<Hoa3d, T>::getNumberOfHarmonics(); i++)
            {
                m_normalization[i] = ProcessorHarmonics<Hoa3d, T>::getHarmonicSemiNormalization(i);
            }
            setMute(false);
            setAzimuth(0.);
            setElevation(0.);
        }

        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
        ~EncoderBasic() hoa_noexcept
        {
            Signal<T>::free(m_normalization);
        }

        //! This method mute or unmute.
        /**	Mute or unmute.
         @param     muted	The mute state.
         */
        inline void setMute(const bool muted) hoa_noexcept { m_muted = muted; }

        //! This method retrieve the mute or unmute state of a source.
        /**	Get the Mute state of a source.
         @return    The mute state of the source.
         */
        inline bool getMute() const hoa_noexcept { return m_muted; }

        //! This method set the angle of azimuth.
        /**	The angle of azimuth in radian and you should prefer to use it between 0 and 2 π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield.
         @param     azimuth	The azimuth.
         @see       setElevation()
         */
        inline void setAzimuth(const T azimuth) hoa_noexcept
        {
            m_azimuth = azimuth;
            m_cos_phi = std::cos(m_azimuth);
            m_sin_phi = std::sin(m_azimuth);
        }

        //! Get the azimuth angle
        /** The method returns the azimuth between 0 and 2π.
         @return     The azimuth.
         */
        inline T getAzimuth() const hoa_noexcept { return Math<T>::wrap_twopi(m_azimuth); }

        //! This method set the angle of elevation.
        /**	The angle of elevation in radian and you should prefer to use it between -π and π to avoid recursive wrapping of the value. The direction of rotation is from bottom to the top. The 0 radian is centered at the "front" of the soundfield, then π/2 is at the top, -π/2 is at the bottom and π is behind. Note that if the angle of elevation is between π/2 and 3*π/2, the azimuth is reversed.
         @param     elevation The elevation.

         */
        inline void setElevation(const T elevation) hoa_noexcept
        {
            m_elevation = Math<T>::wrap_pi(elevation);
            m_cos_theta = T(std::cos(HOA_PI2 + m_elevation));
            m_sqrt_rmin = T(std::sqrt(1 - m_cos_theta * m_cos_theta));
        }

        //!	Get the elevation angle
        /** The method returns the elevation between -π and π.
         @return     The elevation.
         */
        inline T getElevation()  const hoa_noexcept { return m_elevation; }

        //! This method performs the encoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics. For the elevation, the function uses three recurrence formulas :
         \f[P(l, l)(x) = (-1)^l \times (2l - 1)!! \times (1 - x^2)^{0.5l}\f]
         \f[P(l + 1, l)(x) = x \times (2l + 1) \times P(l, l)\f]
         \f[P(l + 1, m)(x) = \frac{(2l + 1) \times x \times P(m, l) - (l + m) \times P(m, l - 1)}{(l - m + 1)}\f]
         with \f$0 \leq l\f$ and \f$-l \leq m \leq +l\f$ and \f[P(0, 0] = 1\f]
         // P(l+1,l+1)   = -(2l+1)√(1-x²)P(l,l)
         // P(l,l+1)     = x(2l+1)P(l,l)
         // P(l+1,m) = (x(2l+1)P(l,m) - (l+m)P(l-1,m))/(l-m+1)
         // P(l+1,0) = (x(2l+1)P(l,0) - l * P(l-1,0))/(l+1)
         For the azimuth, the function uses formulas : sin(lϕ) and cos(lϕ)
         @param     input    The input sample.
         @param     outputs  The outputs array.
         */
        void process(const T* input, T* outputs) hoa_noexcept hoa_override
        {
            if(!m_muted)
            {
                const size_t order = ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder();
                const T cos_theta = m_cos_theta;
                const T sqr_theta = -m_sqrt_rmin;
                const T cos_phi   = (m_elevation >= -HOA_PI2 && m_elevation <= HOA_PI2) ? m_cos_phi : -m_cos_phi;
                const T sin_phi   = (m_elevation >= -HOA_PI2 && m_elevation <= HOA_PI2) ? m_sin_phi : -m_sin_phi;
                const T* norm     = m_normalization;

                // Elevation
                T leg_l2 = 1.;
                T leg_l1 = cos_theta;
                T tleg_l;
                T pleg_l = leg_l2;

                // Azimtuh
                T cos_x = cos_phi;
                T sin_x = sin_phi;
                T tcos_x = cos_x;

                // For m[0] and l{0...N}
                *(outputs)      = (*input);                 // Hamonic [0, 0]
                *(outputs+2)    = (*input) * leg_l1;      // Hamonic [1, 0]
                size_t index = 6;
                for(size_t i = 2; i <= order; i++, index += 2 * i)
                {
                    tleg_l  = (cos_theta * leg_l1 * (T)(2 * (i - 1) + 1) - (T)(i - 1) * leg_l2) / (T)(i);
                    leg_l2  = leg_l1;
                    leg_l1   = tleg_l;
                    *(outputs+index) = (*input) * leg_l1;   // Hamonic [i, 0]
                }

                // For m{1...N-1} and l{m...N}
                index = 1;
                for(size_t i = 1; i < order; i++, index = i * i)
                {
                    size_t inc = 2;
                    leg_l2 = sqr_theta * pleg_l * (T)(2 * (i - 1) + 1);
                    leg_l1 = cos_theta  * leg_l2 * (T)(2 * i + 1);
                    pleg_l = leg_l2;

                    *(outputs+index) = (*input) * leg_l2 * sin_x * *(norm+index);   // Hamonic [i,-i]
                    index += 2*i;
                    *(outputs+index) = (*input) * leg_l2 * cos_x * *(norm+index);   // Hamonic [i, i]
                    index += inc;

                    *(outputs+index) = (*input) * leg_l1 * sin_x * *(norm+index);   // Hamonic [i+1,-i]
                    index += 2*i;
                    *(outputs+index) = (*input) * leg_l1 * cos_x * *(norm+index);   //Hamonic [i+1, i]
                    inc += 2;
                    index += inc;

                    for(size_t j = i + 2; j <= order; j++)
                    {
                        tleg_l  = (cos_theta * leg_l1 * (T)(2 * (j - 1) + 1) - (T)(j - 1 + i) * leg_l2) / (T)(j - i);
                        leg_l2  = leg_l1;
                        leg_l1   = tleg_l;

                        *(outputs+index) = (*input) * leg_l1 * sin_x * *(norm+index);   //Hamonic [j,-i]
                        index += 2*i;
                        *(outputs+index) = (*input) * leg_l1 * cos_x * *(norm+index);   //Hamonic [j, i]
                        inc += 2;
                        index += inc;
                    }
                    cos_x   = tcos_x * cos_phi - sin_x * sin_phi;
                    sin_x   = tcos_x * sin_phi + sin_x * cos_phi;
                    tcos_x  = cos_x;
                }

                // For m[N] and l[N]
                index = order * order;
                leg_l2 = sqr_theta * pleg_l * (T)(2 * (order - 1) + 1);
                *(outputs+index) = (*input) * leg_l2 * sin_x * *(norm+index);
                index += 2 * order;
                *(outputs+index) = (*input) * leg_l2 * cos_x * *(norm+index);
            }
            else
            {
                for(size_t i = 0; i < ProcessorHarmonics<Hoa3d, T>::getNumberOfHarmonics(); i++)
                {
                    (*outputs++) = 0.;
                }
            }
        }

        //! This method performs the encoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics. For the elevation, the function uses three recurrence formulas :
         \f[P(l, l)(x) = (-1)^l \times (2l - 1)!! \times (1 - x^2)^{0.5l}\f]
         \f[P(l + 1, l)(x) = x \times (2l + 1) \times P(l, l)\f]
         \f[P(l + 1, m)(x) = \frac{(2l + 1) \times x \times P(m, l) - (l + m) \times P(m, l - 1)}{(l - m + 1)}\f]
         with \f$0 \leq l\f$ and \f$-l \leq m \leq +l\f$ and \f[P(0, 0] = 1\f]
         // P(l+1,l+1)   = -(2l+1)√(1-x²)P(l,l)
         // P(l,l+1)     = x(2l+1)P(l,l)
         // P(l+1,m) = (x(2l+1)P(l,m) - (l+m)P(l-1,m))/(l-m+1)
         // P(l+1,0) = (x(2l+1)P(l,0) - l * P(l-1,0))/(l+1)
         For the azimuth, the function uses formulas : sin(lϕ) and cos(lϕ)
         @param     input    The input sample.
         @param     outputs  The outputs array.
         */
        void processAdd(const T* input, T* outputs) hoa_noexcept
        {
            if(!m_muted)
            {
                const size_t order = ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder();
                const T cos_theta = m_cos_theta;
                const T sqr_theta = -m_sqrt_rmin;
                const T cos_phi   = (m_elevation >= -HOA_PI2 && m_elevation <= HOA_PI2) ? m_cos_phi : -m_cos_phi;
                const T sin_phi   = (m_elevation >= -HOA_PI2 && m_elevation <= HOA_PI2) ? m_sin_phi : -m_sin_phi;
                const T* norm     = m_normalization;

                // Elevation
                T leg_l2 = 1.;
                T leg_l1 = cos_theta;
                T tleg_l;
                T pleg_l = leg_l2;

                // Azimtuh
                T cos_x = cos_phi;
                T sin_x = sin_phi;
                T tcos_x = cos_x;

                // For m[0] and l{0...N}
                *(outputs)      += (*input);                 // Hamonic [0, 0]
                *(outputs+2)    += (*input) * leg_l1;      // Hamonic [1, 0]
                size_t index = 6;
                for(size_t i = 2; i <= order; i++, index += 2 * i)
                {
                    tleg_l      = (cos_theta * leg_l1 * (T)(2 * (i - 1) + 1) - (T)(i - 1) * leg_l2) / (T)(i);
                    leg_l2      = leg_l1;
                    leg_l1      = tleg_l;
                    *(outputs+index) += (*input) * leg_l1;   // Hamonic [i, 0]
                }

                // For m{1...N-1} and l{m...N}
                index = 1;
                for(size_t i = 1; i < order; i++, index = i * i)
                {
                    size_t inc = 2;
                    leg_l2 = sqr_theta * pleg_l * (T)(2 * (i - 1) + 1);
                    leg_l1 = cos_theta  * leg_l2 * (T)(2 * i + 1);
                    pleg_l = leg_l2;

                    *(outputs+index) += (*input) * leg_l2 * sin_x * *(norm+index);   // Hamonic [i,-i]
                    index += 2*i;
                    *(outputs+index) += (*input) * leg_l2 * cos_x * *(norm+index);   // Hamonic [i, i]
                    index += inc;

                    *(outputs+index) += (*input) * leg_l1 * sin_x * *(norm+index);   // Hamonic [i+1,-i]
                    index += 2*i;
                    *(outputs+index) += (*input) * leg_l1 * cos_x * *(norm+index);   //Hamonic [i+1, i]
                    inc += 2;
                    index += inc;

                    for(size_t j = i + 2; j <= order; j++)
                    {
                        tleg_l  = (cos_theta * leg_l1 * (T)(2 * (j - 1) + 1) - (T)(j - 1 + i) * leg_l2) / (T)(j - i);
                        leg_l2  = leg_l1;
                        leg_l1   = tleg_l;

                        *(outputs+index) += (*input) * leg_l1 * sin_x * *(norm+index);   //Hamonic [j,-i]
                        index += 2*i;
                        *(outputs+index) += (*input) * leg_l1 * cos_x * *(norm+index);   //Hamonic [j, i]
                        inc += 2;
                        index += inc;
                    }
                    cos_x   = tcos_x * cos_phi - sin_x * sin_phi;
                    sin_x   = tcos_x * sin_phi + sin_x * cos_phi;
                    tcos_x  = cos_x;
                }

                // For m[N] and l[N]
                index = order * order;
                leg_l2 = sqr_theta * pleg_l * (T)(2 * (order - 1) + 1);
                *(outputs+index) += (*input) * leg_l2 * sin_x * *(norm+index);
                index += 2 * order;
                *(outputs+index) += (*input) * leg_l2 * cos_x * *(norm+index);
            }
            else
            {
                for(size_t i = 0; i < ProcessorHarmonics<Hoa3d, T>::getNumberOfHarmonics(); i++)
                {
                    (*outputs++) = 0.;
                }
            }
        }
    };

    //! The dc encoder class generates the harmonics for one signal according to an azimuth, an elevation and a radius.
    /** The dc encoder should be used to encode a signal in the harmonics domain depending on an order of decomposition. It allows to control the azimuth, the elevation and the radius of the signal. The distance compensation is performed with the simulation of fractional orders when the signal is inside the ambisonic circle or sphere and with gain attenuation when the signal is outside the ambisonics circle or sphere.
     */
    template <typename T> class EncoderDC<Hoa3d, T> : public Encoder<Hoa3d, T>
    {
    private:
        T  m_azimuth;
        T  m_elevation;
        T  m_cos_phi;
        T  m_sin_phi;
        T  m_cos_theta;
        T  m_sqrt_rmin;
        T  m_radius;
        T* m_normalization;
        T* m_distance;
        bool m_muted;
    public:

        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes circular harmonics coefficients for the encoding. The order must be at least 1.
         @param     order	The order.
         */
        EncoderDC(const size_t order) hoa_noexcept : Encoder<Hoa3d, T>(order)
        {
            m_normalization = Signal<T>::alloc(ProcessorHarmonics<Hoa3d, T>::getNumberOfHarmonics());
            for(size_t i = 0; i < ProcessorHarmonics<Hoa3d, T>::getNumberOfHarmonics(); i++)
            {
                m_normalization[i] = ProcessorHarmonics<Hoa3d, T>::getHarmonicSemiNormalization(i);
            }
            m_distance = Signal<T>::alloc(ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder());
            setMute(false);
            setAzimuth(0.);
            setElevation(0.);
            setRadius(1.);
        }

        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
        ~EncoderDC() hoa_noexcept
        {
            Signal<T>::free(m_normalization);
            Signal<T>::free(m_distance);

        }

        //! This method set the azimuth.
        /**	The azimuth in radian and you should prefer to use it between 0 and 2π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield.
         @param     azimuth	The azimuth.
         */
        inline void setAzimuth(const T azimuth) hoa_noexcept
        {
            m_azimuth = azimuth;
            m_cos_phi = std::cos(m_azimuth);
            m_sin_phi = std::sin(m_azimuth);
        }

        //! Get the azimuth.
        /** The method returns the azimuth between 0 and 2π.
         @return     The azimuth.
         */
        inline T getAzimuth() const hoa_noexcept
        {
            return Math<T>::wrap_twopi(m_azimuth);
        }

        //! This method set the angle of elevation.
        /**	The angle of elevation in radian and you should prefer to use it between -π and π to avoid recursive wrapping of the value. The direction of rotation is from bottom to the top. The 0 radian is centered at the "front" of the soundfield, then π/2 is at the top, -π/2 is at the bottom and π is behind. Note that if the angle of elevation is between π/2 and 3*π/2, the azimuth is reversed.
         @param     elevation The elevation.

         */
        inline void setElevation(const T elevation) hoa_noexcept
        {
            m_elevation = elevation;
            m_cos_theta = std::cos(HOA_PI2 + m_elevation);
            m_sqrt_rmin = std::sqrt(1 - m_cos_theta * m_cos_theta);
        }

        //!	Get the elevation angle
        /** The method returns the elevation between -π and π.
         @return     The elevation.
         */
        inline T getElevation()  const hoa_noexcept
        {
            return Math<T>::wrap_pi(m_elevation);
        }

        //! This method set the radius.
        /**	The radius is between 0 and infinity. At 0, the source is in the center of the ambisonic circle and at 1, the source is at the limit of the ambisonic circle. Over 1, the source get away the ambisonic circle.
         @param     radius   The radius.
         @see       setAzimuth()
         */
        inline void setRadius(const T radius) hoa_noexcept
        {
            m_radius = max(radius, (T)0.);
            T factor, gain, dist;
            if(m_radius < 1.)
            {
                factor      = (1. - m_radius) * HOA_PI;
                gain        = (sin(factor - HOA_PI2) + 1.) * 0.5;
                dist        = 1.;
            }
            else
            {
                factor      = 0;
                gain        = 0;
                dist        = 1. / radius;
            }

            const T gain1   = (gain * ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder());
            m_distance[0] = (gain1 + dist);
            m_distance[1] = (cos(Math<T>::clip(factor, 0., HOA_PI)) + 1.) * 0.5 * ((gain1 - gain) + dist);

            for(size_t i = 2; i <= ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder(); i++)
            {
                const T gain2   = (gain * (ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder() - i) + dist);
                const T factor1 = (cos(Math<T>::clip(factor * i, 0., HOA_PI)) + 1.) * 0.5;
                m_distance[i]   = factor1 * gain2;
            }
        }

        //! Get the azimuth.
        /** The method returns the azimuth between 0 and 2π.
         @return     The azimuth.
         */
        inline T getRadius() const hoa_noexcept
        {
            return m_radius;
        }

        //! This method mute or unmute.
        /**	Mute or unmute.
         @param     muted	The mute state.
         */
        inline void setMute(const bool muted) hoa_noexcept
        {
            m_muted = muted;
        }

        //! This method retrieve the mute or unmute state of a source.
        /**	Get the Mute state of a source.
         @return    The mute state of the source.
         */
        inline bool getMute() const hoa_noexcept
        {
            return m_muted;
        }

        //! This method performs the encoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
         // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
         @param     input	The input sample.
         @param     outputs The output array.
         */
        void process(const T* input, T* outputs) hoa_noexcept hoa_override
        {
            if(!m_muted)
            {
                const size_t order = ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder();
                const T cos_theta = m_cos_theta;
                const T sqr_theta = -m_sqrt_rmin;
                const T cos_phi   = (m_elevation >= -HOA_PI2 && m_elevation <= HOA_PI2) ? m_cos_phi : -m_cos_phi;
                const T sin_phi   = (m_elevation >= -HOA_PI2 && m_elevation <= HOA_PI2) ? m_sin_phi : -m_sin_phi;
                const T* norm     = m_normalization;
                const T* dist     = m_distance;

                // Elevation
                T leg_l2 = 1.;
                T leg_l1 = cos_theta;
                T tleg_l;
                T pleg_l = leg_l2;

                // Azimtuh
                T cos_x = cos_phi;
                T sin_x = sin_phi;
                T tcos_x = cos_x;

                // For m[0] and l{0...N}
                *(outputs)      = (*input) * *(dist);                 // Hamonic [0, 0]
                *(outputs+2)    = (*input) * leg_l1 * *(dist+1);      // Hamonic [1, 0]
                size_t index = 6;
                for(size_t i = 2; i <= order; i++, index += 2 * i)
                {
                    tleg_l  = (cos_theta * leg_l1 * (T)(2 * (i - 1) + 1) - (T)(i - 1) * leg_l2) / (T)(i);
                    leg_l2  = leg_l1;
                    leg_l1   = tleg_l;
                    *(outputs+index) = (*input) * leg_l1 * *(dist+i);   // Hamonic [i, 0]
                }

                // For m{1...N-1} and l{m...N}
                index = 1;
                for(size_t i = 1; i < order; i++, index = i * i)
                {
                    size_t inc = 2;
                    leg_l2 = sqr_theta * pleg_l * (T)(2 * (i - 1) + 1);
                    leg_l1 = cos_theta  * leg_l2 * (T)(2 * i + 1);
                    pleg_l = leg_l2;

                    *(outputs+index) = (*input) * leg_l2 * sin_x * *(norm+index) * *(dist+i);   // Hamonic [i,-i]
                    index += 2*i;
                    *(outputs+index) = (*input) * leg_l2 * cos_x * *(norm+index) * *(dist+i);   // Hamonic [i, i]
                    index += inc;

                    *(outputs+index) = (*input) * leg_l1 * sin_x * *(norm+index) * *(dist+i+1);   // Hamonic [i+1,-i]
                    index += 2*i;
                    *(outputs+index) = (*input) * leg_l1 * cos_x * *(norm+index) * *(dist+i+1);   //Hamonic [i+1, i]
                    inc += 2;
                    index += inc;

                    for(size_t j = i + 2; j <= order; j++)
                    {
                        tleg_l  = (cos_theta * leg_l1 * (T)(2 * (j - 1) + 1) - (T)(j - 1 + i) * leg_l2) / (T)(j - i);
                        leg_l2  = leg_l1;
                        leg_l1   = tleg_l;

                        *(outputs+index) = (*input) * leg_l1 * sin_x * *(norm+index) * *(dist+j);   //Hamonic [j,-i]
                        index += 2*i;
                        *(outputs+index) = (*input) * leg_l1 * cos_x * *(norm+index) * *(dist+j);   //Hamonic [j, i]
                        inc += 2;
                        index += inc;
                    }
                    cos_x   = tcos_x * cos_phi - sin_x * sin_phi;
                    sin_x   = tcos_x * sin_phi + sin_x * cos_phi;
                    tcos_x  = cos_x;
                }

                // For m[N] and l[N]
                index = order * order;
                leg_l2 = sqr_theta * pleg_l * (T)(2 * (order - 1) + 1);
                *(outputs+index) = (*input) * leg_l2 * sin_x * *(norm+index) * *(dist+order);
                index += 2 * order;
                *(outputs+index) = (*input) * leg_l2 * cos_x * *(norm+index) * *(dist+order);
            }
            else
            {
                for(size_t i = 0; i < ProcessorHarmonics<Hoa3d, T>::getNumberOfHarmonics(); i++)
                {
                    (*outputs++) = 0.;
                }
            }
        }

        //! This method performs the encoding.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
         // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
         @param     input	The input sample.
         @param     outputs The output array.
         */
        void processAdd(const T* input, T* outputs) hoa_noexcept
        {
            if(!m_muted)
            {
                const size_t order = ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder();
                const T cos_theta = m_cos_theta;
                const T sqr_theta = -m_sqrt_rmin;
                const T cos_phi   = (m_elevation >= -HOA_PI2 && m_elevation <= HOA_PI2) ? m_cos_phi : -m_cos_phi;
                const T sin_phi   = (m_elevation >= -HOA_PI2 && m_elevation <= HOA_PI2) ? m_sin_phi : -m_sin_phi;
                const T* norm     = m_normalization;
                const T* dist     = m_distance;

                // Elevation
                T leg_l2 = 1.;
                T leg_l1 = cos_theta;
                T tleg_l;
                T pleg_l = leg_l2;

                // Azimtuh
                T cos_x = cos_phi;
                T sin_x = sin_phi;
                T tcos_x = cos_x;

                // For m[0] and l{0...N}
                *(outputs)      += (*input) * *(dist);                 // Hamonic [0, 0]
                *(outputs+2)    += (*input) * leg_l1 * *(dist+1);      // Hamonic [1, 0]
                size_t index = 6;
                for(size_t i = 2; i <= order; i++, index += 2 * i)
                {
                    tleg_l  = (cos_theta * leg_l1 * (T)(2 * (i - 1) + 1) - (T)(i - 1) * leg_l2) / (T)(i);
                    leg_l2  = leg_l1;
                    leg_l1   = tleg_l;
                    *(outputs+index) += (*input) * leg_l1 * *(dist+i);   // Hamonic [i, 0]
                }

                // For m{1...N-1} and l{m...N}
                index = 1;
                for(size_t i = 1; i < order; i++, index = i * i)
                {
                    size_t inc = 2;
                    leg_l2 = sqr_theta * pleg_l * (T)(2 * (i - 1) + 1);
                    leg_l1 = cos_theta  * leg_l2 * (T)(2 * i + 1);
                    pleg_l = leg_l2;

                    *(outputs+index) += (*input) * leg_l2 * sin_x * *(norm+index) * *(dist+i);   // Hamonic [i,-i]
                    index += 2*i;
                    *(outputs+index) += (*input) * leg_l2 * cos_x * *(norm+index) * *(dist+i);   // Hamonic [i, i]
                    index += inc;

                    *(outputs+index) += (*input) * leg_l1 * sin_x * *(norm+index) * *(dist+i+1);   // Hamonic [i+1,-i]
                    index += 2*i;
                    *(outputs+index) += (*input) * leg_l1 * cos_x * *(norm+index) * *(dist+i+1);   //Hamonic [i+1, i]
                    inc += 2;
                    index += inc;

                    for(size_t j = i + 2; j <= order; j++)
                    {
                        tleg_l  = (cos_theta * leg_l1 * (T)(2 * (j - 1) + 1) - (T)(j - 1 + i) * leg_l2) / (T)(j - i);
                        leg_l2  = leg_l1;
                        leg_l1   = tleg_l;

                        *(outputs+index) += (*input) * leg_l1 * sin_x * *(norm+index) * *(dist+j);   //Hamonic [j,-i]
                        index += 2*i;
                        *(outputs+index) += (*input) * leg_l1 * cos_x * *(norm+index) * *(dist+j);   //Hamonic [j, i]
                        inc += 2;
                        index += inc;
                    }
                    cos_x   = tcos_x * cos_phi - sin_x * sin_phi;
                    sin_x   = tcos_x * sin_phi + sin_x * cos_phi;
                    tcos_x  = cos_x;
                }

                // For m[N] and l[N]
                index = order * order;
                leg_l2 = sqr_theta * pleg_l * (T)(2 * (order - 1) + 1);
                *(outputs+index) += (*input) * leg_l2 * sin_x * *(norm+index) * *(dist+order);
                index += 2 * order;
                *(outputs+index) += (*input) * leg_l2 * cos_x * *(norm+index) * *(dist+order);
            }
        }
    };

    //! The multi encoder class generates the harmonics for several signals according to an azimuth, an elevation and a radius for each one.
    /** The multi encoder should be used to encode several signals in the harmonics domain depending on an order of decomposition. It allows to control the azimuth, the elevation and the radius of each signal. The class uses a set of dc encoders.
     */
    template <typename T> class EncoderMulti<Hoa3d, T> : public Encoder<Hoa3d, T>
    {
    private:
        const size_t                     m_number_of_sources;
        std::vector<EncoderDC<Hoa3d, T> *> m_encoders;
    public:

        //! The map constructor.
        /**	The map constructor allocates and initialize the member values and classes depending on a order of decomposition and the number of sources. The order and the number of sources must be at least 1.

         @param     order            The order.
         @param     numberOfSources	The number of sources.
         */
        EncoderMulti(const size_t order, size_t numberOfSources) hoa_noexcept : Encoder<Hoa3d, T>(order),
        m_number_of_sources(numberOfSources)
        {
            for(size_t i = 0; i < m_number_of_sources; i++)
            {
                m_encoders.push_back(new EncoderDC<Hoa3d, T>(order));
            }
        }

        //! The map destructor.
        /**	The map destructor free the memory and deallocate the member classes.
         */
        ~EncoderMulti() hoa_noexcept
        {
            for(size_t i = 0; i < m_number_of_sources; i++)
            {
                delete m_encoders[i];
            }
            m_encoders.clear();
        }

        //! This method retrieve the number of sources.
        /** Retrieve the number of sources.

         @return The number of sources.
         */
        inline size_t getNumberOfSources() const hoa_noexcept
        {
            return m_number_of_sources;
        }

        //! This method set the angle of azimuth of a source.
        /**	The angle of azimuth in radian and you should prefer to use it between 0 and 2 π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield. The index must be between 0 and the number of sources - 1.

         @param     index	The index of the source.
         @param     azimuth	The azimuth.
         @see       setRadius()
         */
        inline void setAzimuth(const size_t index, const T azimuth) hoa_noexcept
        {
            m_encoders[index]->setAzimuth(azimuth);
        }

        //! This method set the angle of azimuth of a source.
        /**	The angle of azimuth in radian and you should prefer to use it between 0 and 2 π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield. The index must be between 0 and the number of sources - 1.

         @param     index	The index of the source.
         @param     azimuth	The azimuth.
         @see       setRadius()
         */
        inline void setElevation(const size_t index, const T elevation) hoa_noexcept
        {
            m_encoders[index]->setElevation(elevation);
        }

        //! This method set the radius of a source.
        /**	The radius is between 0 and infinity. At 0, the source is in the center of the ambisonic circle and at 1, the source is at the limit of the ambisonic circle. Over 1, the source get away the ambisonic circle. The index must be between 0 and the number of sources - 1.

         @param     index	The index of the source.
         @param     radius   The radius.
         @see       setAzimuth()
         */
        inline void setRadius(const size_t index, const T radius) hoa_noexcept
        {
            m_encoders[index]->setRadius(radius);
        }

        //! This method mute or unmute a source.
        /**	Mute or unmute a source with a boolean value. The index must be between 0 and the number of sources - 1.

         @param     index	The index of the source.
         @param     muted	The mute state.
         */
        inline void setMute(const size_t index, const bool muted) hoa_noexcept
        {
            m_encoders[index]->setMute(muted);
        }

        //! This method retrieve the azimuth of a source.
        /** Retrieve the azimuth of a source.

         @param     index	The index of the source.
         @return The azimuth of the source if the source exists, otherwise the function generates an error.
         */
        inline T getAzimuth(const size_t index) const hoa_noexcept
        {
            return m_encoders[index]->getAzimuth();
        }

        //! This method retrieve the elevation of a source.
        /** Retrieve the elevation of a source.

         @param     index	The index of the source.
         @return The elevation of the source if the source exists, otherwise the function generates an error.
         */
        inline T getElevation(const size_t index) const hoa_noexcept
        {
            return m_encoders[index]->getElevation();
        }

        //! This method retrieve the radius of a source.
        /** Retrieve the radius of a source.

         @param     index	The index of the source.
         @return The radius of the source if the source exists, otherwise the function generates an error.
         */
        inline T getRadius(const size_t index) const hoa_noexcept
        {
            return m_encoders[index]->getRadius();
        }

        //! This method retrieve the mute or unmute state of a source.
        /**	Get the Mute state of a source.

         @param     index	The index of the source.
         @return    The mute state of the source if the source exists, otherwise the function generates an error.
         @see       setMute()
         */
        inline bool getMute(const size_t index) const hoa_noexcept
        {
            return m_encoders[index]->getMute();
        }


        //! This method performs the encoding with distance compensation.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The input array contains the samples of the sources and the minimum size should be the number of sources. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     input  The input array.
         @param     outputs The outputs array.
         */
        inline void process(const T* input, T* outputs) hoa_noexcept hoa_override
        {
            m_encoders[0]->process(input, outputs);
            for(size_t i = 1; i < m_number_of_sources; i++)
            {
                m_encoders[i]->processAdd(++input, outputs);
            }
        }
    };
}

#endif
