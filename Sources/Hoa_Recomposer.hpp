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

#ifndef DEF_HOA_RECOMPOSER_LIGHT
#define DEF_HOA_RECOMPOSER_LIGHT

#include "Encoder.hpp"
#include "Planewaves.hpp"

namespace hoa
{
    //! The ambisonic recomposer.
    /** The recomposer should be in the planewaves domain to come back the the circular harmonics domain. The recomposition is similar to the several encoding except that we consider planewaves (or virtual microphones) instead of sources. The number of channels (or planewaves) must be a least the number of harmonics, the first angle is 0 radian and the angular distances between the channels are equals.
     */
    template <Dimension D, typename T> class Recomposer : public EncoderBasic<D, T>, public Processor<D, T>::Planewaves
    {
        enum Recomposition
        {
            Fixe    = 0,
            Fisheye = 1,
            Free    = 2
        };
        
    public:
        //! The recomposer constructor.
        /**	The resomposer constructor allocates and initialize the base classes.
         @param     order                   The order
         @param     numberOfPlanewaves      The number of channels.
         */
        Recomposer(size_t order, size_t numberOfPlanewaves) hoa_noexcept :
        EncoderBasic<D, T>(order),
        Processor<D, T>::Planewaves(numberOfPlanewaves) {}
        
        //! The destructor.
        /** The destructor free the memory.
         */
        virtual ~Recomposer() = 0;
        
        //! This method performs the recomposition.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the planewaves samples and the minimum size must be the number of planewaves and the outputs array contains the harmonic samples and the minimum size must be the number of harmonics.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        virtual void process(const T* inputs, T* outputs) hoa_noexcept {}
    };
    
    //! The ambisonic recomposer.
    /** The recomposer should be in the planewaves domain to come back the the circular harmonics domain. The recomposition is similar to the several encoding except that we consider planewaves (or virtual microphones) instead of sources. The number of channels (or planewaves) must be a least the number of harmonics, the first angle is 0 radian and the angular distances between the channels are equals.
     */
    template <Dimension D, typename T> class RecomposerFixe : public Recomposer<D, T>
    {
        enum Recomposition
        {
            Fixe    = 0,
            Fisheye = 1,
            Free    = 2
        };
        
    public:
        //! The recomposer constructor.
        /**	The resomposer constructor allocates and initialize the base classes.
         @param     order                   The order
         @param     numberOfPlanewaves      The number of channels.
         */
        RecomposerFixe(size_t order, size_t numberOfPlanewaves) hoa_noexcept;
        
        //! The destructor.
        /** The destructor free the memory.
         */
        virtual ~RecomposerFixe() = 0;
        
        //! This method performs the recomposition.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the planewaves samples and the minimum size must be the number of planewaves and the outputs array contains the harmonic samples and the minimum size must be the number of harmonics.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        virtual void process(const T* inputs, T* outputs) hoa_noexcept hoa_override;
    };
    
    
    template <Dimension D, typename T> class RecomposerFisheye : public Recomposer<D, T>
    {
    public:
        //! The decoder constructor.
        /**	The decoder constructor allocates and initialize the base classes.
         @param     order                   The order
         @param     numberOfPlanewaves      The number of channels.
         */
        RecomposerFisheye(size_t order, size_t numberOfPlanewaves) hoa_noexcept;
        
        //! The destructor.
        /** The destructor free the memory.
         */
        virtual ~RecomposerFisheye() = 0;
        
        //! Set the fishEye value.
        /**	The fishEye value is between \f$0\f$ and \f$1\f$. At \f$0\f$, the sound field is intact and at \f$1\f$ the sound field is centered in front of the audience.
         @param     fisheye   The fisheye value.
         */
        virtual void setFisheye(const T fisheye) hoa_noexcept = 0;
        
        //! This method performs the recomposition.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the planewaves samples and the minimum size must be the number of planewaves and the outputs array contains the harmonic samples and the minimum size must be the number of harmonics.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        virtual void process(const T* inputs, T* outputs) hoa_noexcept hoa_override;
    };
    
    
    template <Dimension D, typename T> class RecomposerFree : public Recomposer<D, T>
    {
    public:
        //! The decoder constructor.
        /**	The decoder constructor allocates and initialize the base classes.
         @param     order                   The order
         @param     numberOfPlanewaves      The number of channels.
         */
        RecomposerFree(size_t order, size_t numberOfPlanewaves) hoa_noexcept;
        
        //! The destructor.
        /** The destructor free the memory.
         */
        virtual ~RecomposerFree() = 0;
        
        //! Set the azimuth.
        /**	The azimuth value is between \f$0\f$ and \f$2π\f$.
         @param     azim   The azimuth.
         */
        virtual void setAzimuth(const size_t index, const T azim) hoa_noexcept = 0;
        
        //! Set the widening value.
        /**	The the widening value is between \f$0\f$ and \f$1\f$.
         @param     radius   The radius.
         */
        virtual void setWidening(const size_t index, const T radius) hoa_noexcept = 0;
        
        //! Get the azimuth.
        /**	The azimuth value is between \f$0\f$ and \f$2π\f$.
         @param     index   The index of the planewave.
         @return The azimuth value.
         */
        virtual T getAzimuth(const size_t index) const hoa_noexcept;
        
        //! Get the widening value.
        /**	The the widening value is between \f$0\f$ and \f$1\f$.
         @param   index   The index of planewave.
         @return the widening value.
         */
        virtual T getWidening(const size_t index) const hoa_noexcept = 0;
        
        //! This method performs the recomposition.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the planewaves samples and the minimum size must be the number of planewaves and the outputs array contains the harmonic samples and the minimum size must be the number of harmonics.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        virtual void process(const T* inputs, T* outputs) hoa_noexcept hoa_override;
    };
    
    
    
    
    
    
    
    
    
    
    
    
    
    //! The ambisonic recomposer.
    /** The recomposer should be in the planewaves domain to come back the the circular harmonics domain. The recomposition is similar to the several encoding except that we consider planewaves (or virtual microphones) instead of sources. The number of channels (or planewaves) must be a least the number of harmonics, the first angle is 0 radian and the angular distances between the channels are equals.
     */
    template <typename T> class Recomposer<Hoa2d, T> : public EncoderBasic<Hoa2d, T>, public Processor<Hoa2d, T>::Planewaves
    {
    public:
        enum Recomposition
        {
            Fixe    = 0,
            Fisheye = 1,
            Free    = 2
        };
        
    
        //! The recomposer constructor.
        /**	The resomposer constructor allocates and initialize the base classes.
         @param     order                   The order
         @param     numberOfPlanewaves      The number of channels.
         */
        Recomposer(size_t order, size_t numberOfPlanewaves) hoa_noexcept :
        EncoderBasic<Hoa2d, T>(order),
        Processor<Hoa2d, T>::Planewaves(numberOfPlanewaves) {}
        
        //! The destructor.
        /** The destructor free the memory.
         */
        virtual ~Recomposer() hoa_noexcept hoa_default_f
        
        //! This method performs the recomposition.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the planewaves samples and the minimum size must be the number of planewaves and the outputs array contains the harmonic samples and the minimum size must be the number of harmonics.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        virtual void process(const T* inputs, T* outputs) hoa_noexcept hoa_noexcept = 0;
    };
    
    //! The ambisonic recomposer.
    /** The recomposer should be in the planewaves domain to come back the the circular harmonics domain. The recomposition is similar to the several encoding except that we consider planewaves (or virtual microphones) instead of sources. The number of channels (or planewaves) must be a least the number of harmonics, the first angle is 0 radian and the angular distances between the channels are equals.
     */
    template <typename T> class RecomposerFixe<Hoa2d, T> : public Recomposer<Hoa2d, T>
    {
        enum Recomposition
        {
            Fixe    = 0,
            Fisheye = 1,
            Free    = 2
        };
        
    private:
        T* m_matrix;
        
    public:
        //! The recomposer constructor.
        /**	The resomposer constructor allocates and initialize the base classes.
         @param     order                   The order
         @param     numberOfPlanewaves      The number of channels.
         */
        RecomposerFixe(size_t order, size_t numberOfPlanewaves) hoa_noexcept : Recomposer<Hoa2d, T>(order, numberOfPlanewaves)
        {
            const T factor = 1.;
            T* vector   = Signal<T>::alloc(Encoder<Hoa2d, T>::getNumberOfHarmonics());
            m_matrix    = Signal<T>::alloc(Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves() * Encoder<Hoa2d, T>::getNumberOfHarmonics());
            for(size_t i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                EncoderBasic<Hoa2d, T>::setAzimuth(Processor<Hoa2d, T>::Planewaves::getPlanewaveAzimuth(i));
                EncoderBasic<Hoa2d, T>::process(&factor, vector);
                for(size_t j = 0; j < Encoder<Hoa2d, T>::getNumberOfHarmonics(); j++)
                {
                    m_matrix[j * Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves() + i] = vector[j];
                }
            }
            Signal<T>::free(vector);
        }
        
        //! The destructor.
        /** The destructor free the memory.
         */
        ~RecomposerFixe()
        {
            Signal<T>::free(m_matrix);
        }
        
        //! This method performs the recomposition.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the planewaves samples and the minimum size must be the number of planewaves and the outputs array contains the harmonic samples and the minimum size must be the number of harmonics.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        void process(const T* inputs, T* outputs) hoa_noexcept hoa_override
        {
            Signal<T>::mul(Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(), Encoder<Hoa2d, T>::getNumberOfHarmonics(), inputs, m_matrix, outputs);
        }
    };


    template <typename T> class RecomposerFisheye<Hoa2d, T> : public Recomposer<Hoa2d, T>
    {
    private:
        std::vector< EncoderBasic<Hoa2d, T>*>m_encoders;
    public:
        //! The decoder constructor.
        /**	The decoder constructor allocates and initialize the base classes.
         @param     order                   The order
         @param     numberOfPlanewaves      The number of channels.
         */
        RecomposerFisheye(size_t order, size_t numberOfPlanewaves) hoa_noexcept : Recomposer<Hoa2d, T>(order, numberOfPlanewaves)
        {
            for(size_t i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_encoders.push_back(new EncoderBasic<Hoa2d, T>(order));
            }
        }

        //! The destructor.
        /** The destructor free the memory.
         */
        ~RecomposerFisheye()
        {
            for(size_t i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                delete m_encoders[i];
            }
            m_encoders.clear();
        }

        //! Set the fishEye value.
        /**	The fishEye value is between \f$0\f$ and \f$1\f$. At \f$0\f$, the sound field is intact and at \f$1\f$ the sound field is centered in front of the audience.
         @param     fisheye   The fisheye value.
         */
        inline void setFisheye(const T fisheye) hoa_noexcept hoa_override
        {
            T factor = 1. - Math<T>::clip(fisheye, (T)0., (T)1.);
            for(size_t i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                T azimuth = (T)i / (T)Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves() * HOA_2PI;
                if(azimuth < HOA_PI)
                {
                    azimuth *= factor;
                }
                else
                {
                    azimuth = HOA_2PI - ((HOA_2PI - azimuth) * factor);
                }
                m_encoders[i]->setAzimuth(azimuth);
            }
        }

        //! This method performs the recomposition.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the planewaves samples and the minimum size must be the number of planewaves and the outputs array contains the harmonic samples and the minimum size must be the number of harmonics.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        inline void process(const T* inputs, T* outputs) hoa_noexcept hoa_override
        {
            m_encoders[0]->process(inputs, outputs);
            for(size_t i = 1; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_encoders[i]->processAdd(++inputs, outputs);
            }
        }
    };

    template <typename T> class RecomposerFree<Hoa2d, T> : public Recomposer<Hoa2d, T>
    {
    private:
        std::vector< EncoderDC<Hoa2d, T>* >  m_encoders;
    public:
        //! The decoder constructor.
        /**	The decoder constructor allocates and initialize the base classes.
         @param     order                   The order
         @param     numberOfPlanewaves      The number of channels.
         */
        RecomposerFree(size_t order, size_t numberOfPlanewaves) hoa_noexcept : Recomposer<Hoa2d, T>(order, numberOfPlanewaves)
        {
            for(size_t i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_encoders.push_back(new EncoderDC<Hoa2d, T>(order));
                m_encoders[i]->setAzimuth(i * (HOA_2PI / numberOfPlanewaves));
            }
        }

        //! The destructor.
        /** The destructor free the memory.
         */
        ~RecomposerFree()
        {
            for(size_t i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                delete m_encoders[i];
            }
            m_encoders.clear();
        }

        //! Set the azimuth.
        /**	The azimuth value is between \f$0\f$ and \f$2π\f$.
        @param     azim   The azimuth.
         */
        inline void setAzimuth(const size_t index, const T azim) hoa_noexcept hoa_override
        {
            m_encoders[index]->setAzimuth(azim);
        }

        //! Set the widening value.
        /**	The the widening value is between \f$0\f$ and \f$1\f$.
         @param     radius   The radius.
         */
        inline void setWidening(const size_t index, const T radius) hoa_noexcept hoa_override
        {
            m_encoders[index]->setRadius(Math<T>::clip(radius, (T)0, (T)1));
        }

        //! Get the azimuth.
        /**	The azimuth value is between \f$0\f$ and \f$2π\f$.
         @param     index   The index of the planewave.
         @return The azimuth value.
         */
        inline T getAzimuth(const size_t index) const hoa_noexcept hoa_override
        {
            return m_encoders[index]->getAzimuth();
        }

        //! Get the widening value.
        /**	The the widening value is between \f$0\f$ and \f$1\f$.
         @param   index   The index of planewave.
         @return the widening value.
         */
        inline T getWidening(const size_t index) const hoa_noexcept hoa_override
        {
            return m_encoders[index]->getRadius();
        }

        //! This method performs the recomposition.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the planewaves samples and the minimum size must be the number of planewaves and the outputs array contains the harmonic samples and the minimum size must be the number of harmonics.
         @param     inputs  The input array that contains the samples of the harmonics.
         @param     outputs The output array that contains samples destinated to the channels.
         */
        inline void process(const T* inputs, T* outputs) hoa_noexcept hoa_override
        {
            m_encoders[0]->process(inputs, outputs);
            for(size_t i = 1; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_encoders[i]->processAdd(++inputs, outputs);
            }
        }
    };
}

#endif



