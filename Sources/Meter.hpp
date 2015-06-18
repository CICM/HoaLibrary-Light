/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco, Thomas Le Meur & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_METER_LIGHT
#define DEF_HOA_METER_LIGHT

#include "Planewaves.hpp"
#include "Voronoi.hpp"

namespace hoa
{
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    //! The meter class.
    /** The meter should be used to draw an hoa meter.
     */
    template <Dimension D, typename T> class Meter;

    template <typename T> class Meter<Hoa2d, T> : public Processor<Hoa2d, T>::Planewaves
    {
    private:
        ulong   m_ramp;
        ulong   m_vector_size;
        T*      m_channels_peaks;
        T*      m_channels_azimuth_mapped;
        T*      m_channels_azimuth_width;
        ulong*  m_over_leds;

    public:
        //! The meter constructor.
        /**	The meter constructor allocates and initialize the base classes.
         @param     numberOfPlanewaves      The number of channels.
         */
        Meter(ulong numberOfPlanewaves) noexcept :
        Processor<Hoa2d, T>::Planewaves(numberOfPlanewaves)
        {
            m_ramp                      = 0;
            m_vector_size               = 0;
            m_channels_peaks            = Signal<T>::alloc(Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves());
            m_channels_azimuth_width    = Signal<T>::alloc(Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves());
            m_channels_azimuth_mapped   = Signal<T>::alloc(Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves());
            m_over_leds                 = new ulong[Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves()];
            for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_channels_peaks[i] = 0;
                m_over_leds[i]      = 0;
            }
        }

        //! The destructor.
        /** The destructor free the memory.
         */
        ~Meter()
        {
            Signal<T>::free(m_channels_peaks);
            Signal<T>::free(m_channels_azimuth_width);
            Signal<T>::free(m_channels_azimuth_mapped);
            delete [] m_over_leds;
        }

        //! Set the vector size.
        /** Set the vector size.
        @param vectorSize    The new vector size.
         */
        inline void setVectorSize(ulong vectorSize) noexcept
        {
            m_vector_size   = vectorSize;
            m_ramp          = 0;
        }

        //! Get the vector size.
        /** Get the vector size.
        @return The vector size.
         */
        inline ulong getVectorSize() const noexcept
        {
            return m_vector_size;
        }

        //! This method establish a model of a hoa meter.
        /** This method establish a model of a hoa meter.
         */
        void computeRendering()
        {
            if(Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves() == 1)
            {
                m_channels_azimuth_width[0] = HOA_2PI;
                m_channels_azimuth_mapped[0]= 0.;
            }
            else
            {
                vector<Planewave<Hoa2d, T> > channels;
                for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
                {
                    channels.push_back(Planewave<Hoa2d, T>(i, Math<T>::wrap_twopi(Processor<Hoa2d, T>::Planewaves::getPlanewaveAzimuth(i)), 0.));
                }
                std::sort(channels.begin(), channels.end(), Planewave<Hoa2d, T>::sort_azimuth);
                {
                    const T current_angle   = channels[0].getAzimuth(0., 0., 0.);
                    const T previous_angle  = channels[channels.size() - 1].getAzimuth(0., 0., 0.);
                    const T next_angle      = channels[1].getAzimuth(0., 0., 0.);
                    const T previous_portion= (HOA_2PI - previous_angle) + current_angle;
                    const T next_portion    = next_angle - current_angle;
                    m_channels_azimuth_width[channels[0].getIndex()] = (previous_portion + next_portion) * 0.5;
                    m_channels_azimuth_mapped[channels[0].getIndex()]= Math<T>::wrap_twopi((current_angle - previous_portion * 0.5) + m_channels_azimuth_width[channels[0].getIndex()] * 0.5);
                }
                for(ulong i = 1; i < channels.size() - 1; i++)
                {
                    const T current_angle   = channels[i].getAzimuth(0., 0., 0.);
                    const T previous_angle  = channels[i-1].getAzimuth(0., 0., 0.);
                    const T next_angle      = channels[i+1].getAzimuth(0., 0., 0.);
                    const T previous_portion= current_angle - previous_angle;
                    const T next_portion    = next_angle - current_angle;
                    m_channels_azimuth_width[channels[i].getIndex()] = (previous_portion + next_portion) * 0.5;
                    m_channels_azimuth_mapped[channels[i].getIndex()]= Math<T>::wrap_twopi((current_angle - previous_portion * 0.5) + m_channels_azimuth_width[channels[i].getIndex()] * 0.5);
                }
                {
                    const ulong index = channels.size() - 1;
                    const T current_angle   = channels[index].getAzimuth(0., 0., 0.);
                    const T previous_angle  = channels[index - 1].getAzimuth(0., 0., 0.);
                    const T next_angle      = channels[0].getAzimuth(0., 0., 0.);
                    const T previous_portion= current_angle - previous_angle;
                    const T next_portion    = (HOA_2PI - current_angle) + next_angle;
                    m_channels_azimuth_width[channels[index].getIndex()] = (previous_portion + next_portion) * 0.5;
                    m_channels_azimuth_mapped[channels[index].getIndex()]= Math<T>::wrap_twopi((current_angle - previous_portion * 0.5) + m_channels_azimuth_width[channels[index].getIndex()] * 0.5);
                }
                channels.clear();
            }
        }

        //! Get the channel mapped azimuth.
        /** Get the channel mapped azimuth.
        @param index    The index of the channel.
        @return The channel mapped azimuth.
         */
        inline T getPlanewaveAzimuthMapped(const ulong index) const noexcept
        {
            return m_channels_azimuth_mapped[index];
        }

        //! Get the channel width.
        /** Get the channel width.
        @param index    The index of the channel.
        @return The channel width.
         */
        inline T getPlanewaveWidth(const ulong index) const noexcept
        {
            return m_channels_azimuth_width[index];
        }

        //! Get the channel energy.
        /** Get the channel energy in dB.
        @param index    The index of the channel.
        @return The channel energy.
         */
        inline T getPlanewaveEnergy(const ulong index) const noexcept
        {
            if(m_channels_peaks[index] > 0.)
            {
                return 20. * log10(m_channels_peaks[index]);
            }
            else
            {
                return -90.;
            }
        }

        //! Get the channel overLed state.
        /** Get the channel overLed state.
        @param index    The index of the channel. Return true if the channel meter have to be drawn.
        @return The channel overLed state.
         */
        inline bool getPlanewaveOverLed(const ulong index) const noexcept
        {
            return m_over_leds[index];
        }

        //! This method update the overLed state of the channels.
        /** This method update the overLed state of the channels.
        @param value    A no-NULL value to activate the overLed state of a channel
         */
        inline void tick(const ulong time) noexcept
        {
            for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                if(m_channels_peaks[i] >= 1.)
                {
                    m_over_leds[i] = time;
                }
                else if(m_over_leds[i])
                {
                    m_over_leds[i]--;
                }
            }
        }

        //! This method update the signal values.
        /** This method update the signal value for every channel to perform the meter calculation.
        @param input  The input samples.
         */
        inline void process(const T* inputs) noexcept
        {
            if(m_ramp++ == m_vector_size)
            {
                m_ramp = 0;
                for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
                {
                    m_channels_peaks[i] = fabs(*inputs++);
                }
            }
            else
            {
                for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
                {
                    const T peak = fabs(*inputs++);
                    if(peak > m_channels_peaks[i])
                    {
                        m_channels_peaks[i] = peak;
                    }
                }
            }
        }

        //! This method update the meter.
        /** This method update the meter.
        @param input  The input samples.
        @param output The output samples.
         */
        void process(const T* input, T* outputs) noexcept override
        {
            process(input);
        }
    };

    template <typename T> class Meter<Hoa3d, T> : public Processor<Hoa3d, T>::Planewaves
    {
    public:
        typedef typename Voronoi<Hoa3d>::Point Point;
        typedef vector<Point> Path;
    private:
        ulong   m_ramp;
        ulong   m_vector_size;
        T*      m_channels_peaks;
        ulong*  m_over_leds;

        Path*   m_top;
        Path*   m_bottom;

    public:
        //! The meter constructor.
        /**	The meter constructor allocates and initialize the base classes.
         @param     numberOfPlanewaves      The number of channels.
         */
        Meter(const ulong numberOfPlanewaves) noexcept : Processor<Hoa3d, T>::Planewaves(numberOfPlanewaves)
        {
            m_ramp                      = 0;
            m_vector_size               = 0;
            m_channels_peaks            = Signal<T>::alloc(Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves());
            m_over_leds                 = new ulong[Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves()];
            m_top                       = new Path[Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves()];
            m_bottom                    = new Path[Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves()];
            for(ulong i = 0; i < Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_channels_peaks[i] = 0;
                m_over_leds[i]      = 0;
            }
        }

        //! The destructor.
        /** The destructor free the memory.
         */
        ~Meter()
        {
            Signal<T>::free(m_channels_peaks);
            delete [] m_over_leds;
            for(ulong i = 0; i < Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_top[i].clear();
                m_bottom[i].clear();
            }
            delete [] m_top;
            delete [] m_bottom;
        }

        //! Set the vector size.
        /** Set the vector size.
        @param vectorSize    The new vector size.
         */
        inline void setVectorSize(const ulong vectorSize) noexcept
        {
            m_vector_size   = vectorSize;
            m_ramp          = 0;
        }

        //! Get the vector size.
        /** Get the vector size.
        @return The vector size.
         */
        inline ulong getVectorSize() const noexcept
        {
            return m_vector_size;
        }

        //! Get the channel energy.
        /** Get the channel energy in dB.
        @param index    The index of the channel
        @return The channel energy.
         */
        inline T getPlanewaveEnergy(const ulong index) const noexcept
        {
            if(m_channels_peaks[index] > 0.)
            {
                return 20. * log10(m_channels_peaks[index]);
            }
            else
            {
                return -90.;
            }
        }

        //! Get the channel overLed state.
        /** Get the channel overLed state.
        @param index    The index of the channel. Return true if the channel meter have to be drawn.
        @return The channel overLed state.
         */
        inline bool getPlanewaveOverLed(const ulong index) const noexcept
        {
            return m_over_leds[index];
        }

        //! This method update the overLed state of the channels.
        /** This method update the overLed state of the channels.
        @param value    A no-NULL value to activate the overLed state of a channel
         */
        inline void tick(const ulong time) noexcept
        {
            for(ulong i = 0; i < Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                if(m_channels_peaks[i] >= 1.)
                {
                    m_over_leds[i] = time;
                }
                else if(m_over_leds[i])
                {
                    m_over_leds[i]--;
                }
            }
        }

        //! This method update the signal values.
        /** This method update the signal value for every channel to perform the meter calculation.
        @param input  The input samples.
         */
        inline void process(const T* inputs) noexcept
        {
            if(m_ramp++ == m_vector_size)
            {
                m_ramp = 0;
                for(ulong i = 0; i < Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(); i++)
                {
                    m_channels_peaks[i] = fabs(*inputs++);
                }
            }
            else
            {
                for(ulong i = 0; i < Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(); i++)
                {
                    const T peak = fabs(*inputs++);
                    if(peak > m_channels_peaks[i])
                    {
                        m_channels_peaks[i] = peak;
                    }
                }
            }
        }

        //! This method update the meter.
        /** This method update the meter.
        @param input  The input samples.
        @param output The output samples.
         */
        void process(const T* input, T* outputs) noexcept override
        {
            process(input);
        }

        //! This method establish a model of a hoa meter.
        /** This method establish a model of a hoa meter.
         */
        void computeRendering()
        {
            Voronoi<Hoa3d> voronoi;

            for(ulong i = 0; i < Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                voronoi.add(Point(Processor<Hoa3d, T>::Planewaves::getPlanewaveAbscissa(i), Processor<Hoa3d, T>::Planewaves::getPlanewaveOrdinate(i), -Processor<Hoa3d, T>::Planewaves::getPlanewaveHeight(i)));
                m_bottom[i].clear();
            }
            voronoi.compute();
            Path const& bottom = voronoi.getPoints();
            for(ulong i = 0; i < Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                for(ulong j = 0; j < bottom[i].bounds.size(); j++)
                {
                    m_bottom[i].push_back(bottom[i].bounds[j]);
                    m_bottom[i][j].z = -m_bottom[i][j].z;
                }
            }

            voronoi.clear();
            for(ulong i = 0; i < Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                voronoi.add(Point(Processor<Hoa3d, T>::Planewaves::getPlanewaveAbscissa(i), Processor<Hoa3d, T>::Planewaves::getPlanewaveOrdinate(i), Processor<Hoa3d, T>::Planewaves::getPlanewaveHeight(i)));
                m_top[i].clear();
            }
            voronoi.compute();
            Path const& top = voronoi.getPoints();
            for(ulong i = 0; i < Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                for(ulong j = 0; j < top[i].bounds.size(); j++)
                {
                    m_top[i].push_back(top[i].bounds[j]);
                }
            }
        }

        //! Get the Voronoi point of a channel.
        /** Get the Voronoi point of a channel.
        @param index      The index of the channel.
        @param top  The meter view wanted (top or bottom).
        @return The path.
         */
        inline Path const& getPlanewavePath(const ulong index, const bool top) const noexcept
        {
            if(top)
            {
                return m_top[index];
            }
            else
            {
                return m_bottom[index];
            }
        }
    };
#endif
}

#endif



