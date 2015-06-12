/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco, Thomas Le Meur & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_PROCESSOR_LIGHT
#define DEF_HOA_PROCESSOR_LIGHT

#include "Harmonics.hpp"
#include "Planewaves.hpp"

namespace hoa
{
    //! The processor.
    /** The processor owns a set of harmonics or planewaves and performs digital signal processing.
     */
    template <Dimension D, typename T> class Processor
    {
    public:

        //! This method performs the processing.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The input array and the outputs array depends of the template and the processing.
         @param     input  The input array.
         @param     outputs The outputs array.
         */
        virtual void process(const T* input, T* outputs) noexcept = 0;

        class Harmonics;

        class Planewaves;
    };

    //! The harmonic processor.
    /** The harmonic processor owns a set of harmonics depending on the order of decomposition.
     */
    template <Dimension D, typename T> class Processor<D, T>::Harmonics : virtual public Processor<D, T>
    {
    private:

        const ulong                 m_order_of_decomposition;
        const ulong                 m_number_of_harmonics;
        vector< Harmonic<D, T> >    m_harmonics;
    public:

        //! The harmonics constructor.
        /** The harmonics constructor allocates and initializes the general member values depending on a order of decomposition \f$N\f$.
         @param order    The order of decomposition \f$N\f$, must be at least 1.
         */
        Harmonics(const ulong order) noexcept :
        m_order_of_decomposition(order),
        m_number_of_harmonics(Harmonic<D, T>::getNumberOfHarmonics(order))
        {
            for(ulong i = 0; i < m_number_of_harmonics; i++)
            {
                m_harmonics.push_back(Harmonic<D, T>(i));
            }
        }

        //! The harmonics destructor.
        /** The harmonics destructor.
         */
        virtual ~Harmonics() noexcept
        {
            m_harmonics.clear();
        }

        //! Retrieve the order of decomposition.
        /** Retrieve the order of decomposition \f$N\f$.
         @return The order.
         */
        inline ulong getDecompositionOrder() const noexcept
        {
            return m_order_of_decomposition;
        }

        //! Retrieve the number of harmonics.
        /** Retrieve the number of harmonics.
         @return The number of harmonics.
         */
        inline ulong getNumberOfHarmonics() const noexcept
        {
            return m_number_of_harmonics;
        }

        //! Retrieve the degree of an harmonic.
        /** The method retrieves the degrees \f$l\f$of the harmonics are in the range \f$0\f$ to \f$N\f$.
         @param     index	The index of an harmonic.
         @return    The method returns the degree \f$l\f$ of the harmonic.
         @see       getHarmonicOrder()
         @see       getHarmonicName()
         */
        inline ulong getHarmonicDegree(const ulong index) const noexcept
        {
            return m_harmonics[index].getDegree();
        }

        //! Retrieve the order of an harmonic.
        /** The method retrieves the orders \f$m\f$ of the harmonic are in the range \f$-l\f$ to \f$l\f$.
         @param     index	The index of an harmonic.
         @return    The method returns the order \f$m\f$ of the harmonic.
         @see       getHarmonicDegree()
         @see       getHarmonicName()
         */
        inline long getHarmonicOrder(const ulong index) const noexcept
        {
            return m_harmonics[index].getOrder();
        }

        //! Retrieve the index of an harmonic.
        /** The method retrieves the index of an harmonic in the range \f$0\f$ to number of harmonics - 1.
         @param     order	The order an harmonic.
         @return    The method returns the index of the harmonic.
         @see       getHarmonicOrder()
         @see       getHarmonicName()
         */
        inline ulong getHarmonicIndex(const ulong degree, long order) const noexcept
        {
            return Harmonic<D, T>::getIndex(degree, order);
        }

        //! Retrieve the name of an harmonic.
        /** Retrieve the name of an harmonic.
         @param     index	The index of an harmonic.
         @return    The method returns the name of the harmonic that contains its degree and its order.
         @see       getHarmonicDegree()
         @see       getHarmonicOrder()
         */
        inline string getHarmonicName(const ulong index) const noexcept
        {
            return m_harmonics[index].getName();
        }
        
        //! Get the normalization of an harmonic.
        /** The method returns the normalization of an harmonics.
         @param     index	The index of an harmonic.
         @return     The normalization of an harmonics.
         */
        inline T getHarmonicNormalization(const ulong index) const noexcept
        {
            return m_harmonics[index].getNormalization();
        }
        
        //! Get the semi-normalization of an harmonic.
        /** The method returns the semi-normalization of an harmonics.
         @param     index	The index of an harmonic.
         @return    The semi-normalization of the harmonics.
         */
        inline T getHarmonicSemiNormalization(const ulong index) const noexcept
        {
            return m_harmonics[index].getSemiNormalization();
        }

        //! This method performs the processing.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The input array and the outputs array depends of the template and the processing.
         @param     input  The input array.
         @param     outputs The outputs array.
         */
        virtual void process(const T* input, T* outputs) noexcept
        {
            ;
        }
    };


    //! The planewave processor.
    /** The planewave processor owns a set of planewaves.
     */
    template <Dimension D, typename T> class Processor<D, T>::Planewaves : virtual public Processor<D, T>
    {
    private:
        const ulong                 m_number_of_planewaves;
        vector<Planewave<D, T> >    m_planewaves;
        T                           m_rotation_z;
        T                           m_rotation_y;
        T                           m_rotation_x;

    public:

        //! The planewaves constructor.
        /** The planewaves constructor allocates and initializes the general member values depending on a number of planewaves. The number of planewaves must be a least 1.
         @param     numberOfPlanewaves	The number of planewaves.
         */
        Planewaves(const ulong numberOfPlanewaves) noexcept :
        m_number_of_planewaves(numberOfPlanewaves),
        m_rotation_z(0.),
        m_rotation_y(0.),
        m_rotation_x(0.)
        {
#ifndef DOXYGEN_SHOULD_SKIP_THIS
            if(D == Hoa2d)
            {
                for(ulong i = 0; i < m_number_of_planewaves; i++)
                {
                    m_planewaves.push_back(Planewave<D, T>(i+1, (T)i / (m_number_of_planewaves) * HOA_2PI, 0.));
                }
            }
            else
            {
                if(m_number_of_planewaves == 4)
                {
                    const double oh = -(sqrt(2. / 3.) / sqrt(3. / 8.) - 1.);
                    const double hc = sqrt(1. - oh * oh);
                    const double el = asin(oh / sqrt(hc*hc + oh*oh));
                    m_planewaves.push_back(Planewave<D, T>(1ul,  0., HOA_PI2));
                    m_planewaves.push_back(Planewave<D, T>(2ul, 0., el));
                    m_planewaves.push_back(Planewave<D, T>(3ul, HOA_2PI / 3., el));
                    m_planewaves.push_back(Planewave<D, T>(4ul, 2. * HOA_2PI / 3., el));
                }
                else if(m_number_of_planewaves == 6)
                {
                    m_planewaves.push_back(Planewave<D, T>(1ul, 0., HOA_PI2));
                    m_planewaves.push_back(Planewave<D, T>(2ul, 0., 0.));
                    m_planewaves.push_back(Planewave<D, T>(3ul, HOA_PI2, 0.));
                    m_planewaves.push_back(Planewave<D, T>(4ul, 2. * HOA_PI2, 0.));
                    m_planewaves.push_back(Planewave<D, T>(5ul, 3. * HOA_PI2, 0.));
                    m_planewaves.push_back(Planewave<D, T>(6ul, 0., -HOA_PI2));
                }
                else if(m_number_of_planewaves == 8)
                {
                    m_planewaves.push_back(Planewave<D, T>(1ul, 1., 1., 1.));
                    m_planewaves.push_back(Planewave<D, T>(2ul, -1., 1., 1.));
                    m_planewaves.push_back(Planewave<D, T>(3ul, -1., -1., 1.));
                    m_planewaves.push_back(Planewave<D, T>(4ul, 1., -1., 1.));

                    m_planewaves.push_back(Planewave<D, T>(5ul, 1., 1., -1.));
                    m_planewaves.push_back(Planewave<D, T>(6ul, -1., 1., -1.));
                    m_planewaves.push_back(Planewave<D, T>(7ul, -1., -1., -1.));
                    m_planewaves.push_back(Planewave<D, T>(8ul, 1., -1., -1.));
                }
                else if(m_number_of_planewaves == 12)
                {
                    m_planewaves.push_back(Planewave<D, T>(1, 0., HOA_PI2));
                    for(ulong i = 1; i < 6; i++)
                    {
                        m_planewaves.push_back(Planewave<D, T>(i*2, T(i - 1.) / 5. * HOA_2PI, atan(0.5)));
                        m_planewaves.push_back(Planewave<D, T>(i*2+1, T(i - 1.) / 5. * HOA_2PI - HOA_PI / 5., -atan(0.5)));
                    }
                    m_planewaves.push_back(Planewave<D, T>(12, 0., -HOA_PI2));
                }
                else if(m_number_of_planewaves == 20)
                {
                    const T phi = (sqrt(5.) - 1.) / 2.;
                    const T R = 1. / sqrt(3.);
                    const T a = R;
                    const T b = R / phi;
                    const T c = R * phi;
                    ulong index = 1;
                    for(long i = -1; i < 2; i += 2)
                    {
                        for(long j = -1; j < 2; j += 2)
                        {
                            m_planewaves.push_back(Planewave<D, T>(index++, 0., i * c * R, -j * b * R));
                            m_planewaves.push_back(Planewave<D, T>(index++, i * c * R, j * b * R, 0.));
                            m_planewaves.push_back(Planewave<D, T>(index++, i * b * R, 0., -j * c * R));
                            for(long k = -1; k < 2; k += 2)
                            {
                                m_planewaves.push_back(Planewave<D, T>(index++, i * a * R, j * a * R, -k * a * R));
                            }
                        }
                    }
                }
                else
                {
                    if(m_number_of_planewaves % 2)
                    {
                        m_planewaves.push_back(Planewave<D, T>(1, 0., HOA_PI2));
                    }
                    const T     phi     = (sqrt(5.) - 1.) / 4.;
                    const ulong limit   = (m_number_of_planewaves - (m_number_of_planewaves % 2)) / 2;
                    const T     offset  = 1. / T(limit) * HOA_PI;
                    for(ulong i = 0; i < m_number_of_planewaves - (m_number_of_planewaves % 2); i++)
                    {
                        if(i < limit)
                        {
                            m_planewaves.push_back(Planewave<D, T>(i+1, T(i) / T(limit) * HOA_2PI, phi));
                        }
                        else
                        {
                            m_planewaves.push_back(Planewave<D, T>(i+1, T(i - (limit)) / T(limit) * HOA_2PI + offset, phi * 2));
                        }
                    }
                }
            }
#endif
        }

        //! The planewaves destructor.
        /** The planewaves destructor.
         */
        virtual ~Planewaves() noexcept
        {
            m_planewaves.clear();
        }

        //! Retrieve the order of decomposition.
        /** Retrieve the order of decomposition.
         @return The order.
         */
        inline ulong getNumberOfPlanewaves() const noexcept
        {
            return m_number_of_planewaves;
        }

        //! Set the offset of the planewaves.
        /**	Set the azimuth offset of the planewaves in radian.
         @param     x_axe		An azimuth value.
         @param     y_axe		An azimuth value.
         @param     z_axe		An azimuth value.
         */
        inline void setPlanewavesRotation(const T x_axe, const T y_axe, const T z_axe) noexcept
        {
            m_rotation_x = Math<T>::wrap_twopi(x_axe);
            m_rotation_y = Math<T>::wrap_twopi(y_axe);
            m_rotation_z = Math<T>::wrap_twopi(z_axe);
        }

        //! Get the offset of the planewaves.
        /**	Retrieve the azimuth offset of the planewaves in radian.
         @return    The offset of the planewaves.
         */
        inline T getPlanewavesRotationX() const noexcept
        {
            return m_rotation_x;
        }

        //! Get the offset of the planewaves.
        /**	Retrieve the azimuth offset of the planewaves in radian.
         @return    The offset of the planewaves.
         */
        inline T getPlanewavesRotationY() const noexcept
        {
            return m_rotation_y;
        }

        //! Get the offset of the planewaves.
        /**	Retrieve the azimuth offset of the planewaves in radian.
         @return    The offset of the planewaves.
         */
        inline T getPlanewavesRotationZ() const noexcept
        {
            return m_rotation_z;
        }

        //! Get the index of a planewave.
        /** Get the index of a planewave.
         @param      index   The index of the planewave.
         @return     The index of the planewave.
         */
        inline ulong getPlanewaveIndex(const ulong index) noexcept
        {
            return m_planewaves[index].getIndex();
        }

        //! Set the azimuth of a planewave.
        /** Set the azimuth of a planewave. The azimuth of the planewaves is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of planewaves - 1.

         @param     index   The index of the planewave.
         @param     azimuth The new azimuth of the planewave.
         */
        inline void setPlanewaveAzimuth(const ulong index, const T azimuth) noexcept
        {
            m_planewaves[index].setAzimuth(Math<T>::wrap_twopi(azimuth));
        }

        //! Get the azimuth of a planewave.
        /** Get the azimuth of a planewave. The azimuth of the planewaves is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of planewaves - 1.

         @param      index   The index of the planewave.
         @param      rotation   False if you don't want to consider the rotation, otherwise true (default).
         @return     The azimuth of the planewave.
         */
        inline T getPlanewaveAzimuth(const ulong index, const bool rotation = true) const noexcept
        {
            return m_planewaves[index].getAzimuth(rotation ? m_rotation_x : T(0.), rotation ? m_rotation_y : T(0.), rotation ? m_rotation_z : T(0.));
        }

        //! Set the elevation of a planewave.
        /** Set the elevation of a planewave. The elevation of the planewaves is in radian. The maximum index must be the number of planewaves - 1.

         @param      index    The index of the planewave.
         @param      azimuth  The azimuth of the planewave.
         */
        inline void setPlanewaveElevation(const ulong index, const T azimuth) noexcept
        {
            m_planewaves[index].setElevation(Math<T>::wrap_pi(azimuth));
        }

        //! Get the elevation of a planewave.
        /** Get the elevation of a planewave. The elevation of the planewaves is in radian. The maximum index must be the number of planewaves - 1.

         @param      index      The index of the planewave.
         @param      rotation   False if you don't want to consider the rotation, otherwise true (default).
         @return     The elevation of the planewave.
         */
        inline T getPlanewaveElevation(const ulong index, const bool rotation = true) const noexcept
        {
            return m_planewaves[index].getElevation(rotation ? m_rotation_x : T(0.), rotation ? m_rotation_y : T(0.), rotation ? m_rotation_z : T(0.));
        }

        //! Get the abscissa of a planewave.
        /** Get the abscissa of a planewave. The abscissa is between -1 and 1, -1 is the left of the soundfield, 0 is the center of the soundfield and 1 is the right of the soundfield. The maximum index must be the number of planewaves - 1.
         @param     index    The index of the planewave.
         @param      rotation   False if you don't want to consider the rotation, otherwise true (default).
         @return    The abscissa of the planewave.
         */
        inline T getPlanewaveAbscissa(const ulong index, const bool rotation = true) const noexcept
        {
            return m_planewaves[index].getAbscissa(rotation ? m_rotation_x : T(0.), rotation ? m_rotation_y : T(0.), rotation ? m_rotation_z : T(0.));
        }

        //! Get the ordinate of a planewave.
        /** Get the ordinate of a planewave. The ordinate is between -1 and 1, -1 is the back of the soundfield, 0 is the center of the soundfield and 1 is the front of the soundfield. The maximum index must be the number of planewaves - 1.
         @param     index	The index of the planewave.
         @param      rotation   False if you don't want to consider the rotation, otherwise true (default).
         @return    The ordinate of the planewave.
         */
        inline T getPlanewaveOrdinate(const ulong index, const bool rotation = true) const noexcept
        {
            return m_planewaves[index].getOrdinate(rotation ? m_rotation_x : T(0.), rotation ? m_rotation_y : T(0.), rotation ? m_rotation_z : T(0.));
        }

        //! Get the height of a planewave.
        /** Get the ordinate of a planewave. The height is between -1 and 1, -1 is the back of the soundfield, 0 is the center of the soundfield and 1 is the front of the soundfield. The maximum index must be the number of planewaves - 1.
         @param     index	The index of the planewave.
         @param      rotation   False if you don't want to consider the rotation, otherwise true (default).
         @return    The height of the planewave.
         */
        inline T getPlanewaveHeight(const ulong index, const bool rotation = true) const noexcept
        {
            return m_planewaves[index].getHeight(rotation ? m_rotation_x : T(0.), rotation ? m_rotation_y : T(0.), rotation ? m_rotation_z : T(0.));
        }

        //! Get a name for a planewave.
        /** Get a name for a planewave in a string format that will be "Planewave index azimuth (in degrees)".
         @param     index	The index of a planewave.
         @return    The method returns a name for the planewave.
         */
        inline string getPlanewaveName(const ulong index) const noexcept
        {
            return m_planewaves[index].getName();
        }

        //! This method performs the processing.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The input array and the outputs array depends of the template and the processing.
         @param     input  The input array.
         @param     outputs The outputs array.
         */
        virtual void process(const T* input, T* outputs) noexcept
        {
            ;
        }
    };
}

#endif


