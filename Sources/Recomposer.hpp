/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_RECOMPOSER_LIGHT
#define DEF_HOA_RECOMPOSER_LIGHT

#include "Encoder.hpp"
#include "Planewaves.hpp"

namespace hoa
{
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    enum Recomposition
    {
        Fixe    = 0,
        Fisheye = 1,
        Free    = 2
    };
    
    //! The ambisonic recomposer.
    /** The recomposer should be in the planewaves domain to come back the the circular harmonics domain. The recomposition is similar to the several encoding exept that we consider planewaves (or virtual microphones) instead of sources. The number of channels (or planewaves) must be a least the number of harmonics, the first angle is 0 radian and the angular distances between the channels are equals.
     */
    template <Dimension D, typename T, Recomposition M> class Recomposer;
    
    
    template <typename T> class Recomposer<Hoa2d, T, Fixe> : public Encoder<Hoa2d, T>::Basic, public Processor<Hoa2d, T>::Planewaves
    {
    private:
        T* m_matrix;
        
    public:
        
        Recomposer(ulong order, ulong numberOfPlanewaves) noexcept :
        Encoder<Hoa2d, T>::Basic(order),
        Processor<Hoa2d, T>::Planewaves(numberOfPlanewaves)
        {
            const T factor = 1.;
            T* vector   = new T[Encoder<Hoa2d, T>::getNumberOfHarmonics()];
            m_matrix    = new T[Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves() * Encoder<Hoa2d, T>::getNumberOfHarmonics()];
            for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                Encoder<Hoa2d, T>::Basic::setAzimuth(Processor<Hoa2d, T>::Planewaves::getPlanewaveAzimuth(i));
                Encoder<Hoa2d, T>::Basic::process(&factor, vector);
                for(ulong j = 0; j < Encoder<Hoa2d, T>::getNumberOfHarmonics(); j++)
                {
                    m_matrix[j * Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves() + i] = vector[j];
                }
            }
            delete [] vector;
        }
        
        ~Recomposer()
        {
            delete [] m_matrix;
        }
        
        void process(const T* inputs, T* outputs) noexcept override
        {
            Signal<T>::matrix_vector_mul(Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(), Encoder<Hoa2d, T>::getNumberOfHarmonics(), inputs, m_matrix, outputs);
        }
    };

    
    template <typename T> class Recomposer<Hoa2d, T, Fisheye> : public Processor<Hoa2d, T>::Harmonics, public Processor<Hoa2d, T>::Planewaves
    {
    private:
        vector< typename Encoder<Hoa2d, T>::Basic*>m_encoders;
    public:
        
        Recomposer(ulong order, ulong numberOfPlanewaves) noexcept :
        Processor<Hoa2d, T>::Harmonics(order),
        Processor<Hoa2d, T>::Planewaves(numberOfPlanewaves)
        {
            for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_encoders.push_back(new typename Encoder<Hoa2d, T>::Basic(order));
            }
        }
        
        ~Recomposer()
        {
            for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                delete m_encoders[i];
            }
            m_encoders.clear();
        }
    
        inline void setFisheye(const T fisheye) noexcept
        {
            T factor = 1. - Math<T>::clip(fisheye, (T)0., (T)1.);
            for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                T azimuth = (T)i / (T)Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves() * HOA_2PI - HOA_PI;
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
        
        inline void process(const T* inputs, T* outputs) noexcept
        {
            m_encoders[0]->process(inputs, outputs);
            for(ulong i = 1; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_encoders[i]->processAdd(inputs++, outputs);
            }
        }
    };
    
    template <typename T> class Recomposer<Hoa2d, T, Free> : public Processor<Hoa2d, T>::Harmonics, public Processor<Hoa2d, T>::Planewaves
    {
    private:
        vector< typename Encoder<Hoa2d, T>::DC* >  m_encoders;
    public:
        
        Recomposer(ulong order, ulong numberOfPlanewaves) noexcept :
        Processor<Hoa2d, T>::Harmonics(order),
        Processor<Hoa2d, T>::Planewaves(numberOfPlanewaves)
        {
            for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_encoders.push_back(new typename Encoder<Hoa2d, T>::DC(order));
            }
        }
        
        ~Recomposer()
        {
            for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                delete m_encoders[i];
            }
            m_encoders.clear();
        }

        inline void setAzimuth(const ulong index, const T azimuth) noexcept
        {
            m_encoders[index]->setAzimuth(azimuth);
        }
        
        inline void setWidening(const ulong index, const T radius) noexcept
        {
            m_encoders[index]->setRadius(Math<T>::clip(radius, (T)0, (T)1));
        }
        
        inline T getAzimuth(const ulong index) const noexcept
        {
            return m_encoders[index]->getAzimuth();
        }

        inline T getWidening(const ulong index) const noexcept
        {
            return m_encoders[index]->getRadius();
        }

        inline void process(const T* inputs, T* outputs) noexcept
        {
            m_encoders[0]->process(inputs, outputs);
            for(ulong i = 1; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_encoders[i]->processAdd(inputs++, outputs);
            }
        }
    };
#endif
}

#endif



