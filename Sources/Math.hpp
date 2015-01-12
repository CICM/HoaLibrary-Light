/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_MATH_LIGHT
#define DEF_HOA_MATH_LIGHT

#include "Defs.hpp"

namespace hoa
{
    template <typename T> class Math
    {
    public:
        
        //! The clipping function
        /** The function clips a number between boundaries. \n
         If \f$x < min\f$, \f$f(x) = min\f$ else if \f$x > max\f$, \f$f(x) = max\f$ else \f$f(x) = x\f$.
         
         @param     value   The value to clip.
         @param     min     The low boundary.
         @param     max     The high boundary.
         @return    The function return the clipped value.
         */
        static inline T clip(const T& n, const T& lower, const T& upper)
        {
            return std::max(lower, std::min(n, upper));
        }
        
        //! The wrapping function
        /** The function wraps a number between boundarys.
         
         @param     value   The value to wrap.
         @param     low     The low boundary.
         @param     high    The high boundary.
         @return    The function return the wrapped value.
         
         @see    wrap_twopi
         */
        static inline T wrap(const T value, const T low, const T high)
        {
            const T increment = high - low;
            T new_value = value;
            while(new_value < low)
            {
                new_value += increment;
            }
            while(new_value > high)
            {
                new_value -= increment;
            }
            return new_value;
        }
        
        //! The wrapping function over \f$2\pi\f$
        /** The function wraps a number between \f$0\f$ and \f$2\pi\f$.
         
         @param     value   The value to wrap.
         @return    The function return the wrapped value.
         */
        static inline T wrap_twopi(const T& value)
        {
            T new_value = value;
            while(new_value < 0.)
            {
                new_value += HOA_2PI;
            }
            while(new_value >= HOA_2PI)
            {
                new_value -= HOA_2PI;
            }
            return new_value;
        }
        
        //! The wrapping function over \f$2\pi\f$
        /** The function wraps a number between \f$0\f$ and \f$2\pi\f$.
         @param     value   The value to wrap.
         @return    The function return the wrapped value.
         */
        static inline T wrap_pi(const T& value)
        {
            T new_value = value;
            while(new_value < -HOA_PI)
            {
                new_value += HOA_2PI;
            }
            while(new_value >= HOA_PI)
            {
                new_value -= HOA_2PI;
            }
            return new_value;
        }
        
        static inline ulong wrap_ptr(long index, const ulong size)
        {
            while(index < 0)
            {
                index += size;
            }
            while(index >= size)
            {
                index -= size;
            }
            return index;
        }
        
        //! The abscissa converter function.
        /** This function takes a radius and an azimuth and convert them to an abscissa.
         @param     radius		The radius (greather than 0).
         @param     azimuth		The azimuth (between \f$0\f$ and \f$2\pi\f$).
         @param     elevation   The elevation (between \f$-\pi\f$ and \f$\pi\f$).
         @return    The abscissa.
         */
        static inline T abscissa(const T radius, const T azimuth, const T elevation = 0.)
        {
            return radius * cos(azimuth + HOA_PI2) * cos(elevation);
        }
        
        //! The ordinate converter function.
        /** This function takes a radius and an azimuth and convert them to an ordinate.
         @param     radius		The radius (greather than 0).
         @param     azimuth		The azimuth (between \f$0\f$ and \f$2\pi\f$).
         @param     elevation   The elevation (between \f$-\pi\f$ and \f$\pi\f$).
         @return    The ordinate.
         */
        static inline T ordinate(const T radius, const T azimuth, const T elevation = 0.)
        {
            return radius * sin(azimuth + HOA_PI2) * cos(elevation);
        }
        
        //! The height converter function.
        /** This function takes a radius and an azimuth and convert them to an ordinate.
         @param     radius		The radius (greather than 0).
         @param     azimuth		The azimuth (between \f$0\f$ and \f$2\pi\f$).
         @param     elevation   The elevation (between \f$-\pi\f$ and \f$\pi\f$).
         @return    The height.
         */
        static inline T height(const T radius, const T azimuth, const T elevation = 0.)
        {
            return radius * sin(elevation);
        }
        
        //! The azimuth converter function.
        /** This function takes an abscissa and an ordinate and convert them to an azimuth.
         @param     x		The abscissa.
         @param     y		The ordinate.
         @param     z		The height.
         @return    The azimuth.
         */
        static inline T azimuth(const T x, const T y, const T z = 0.)
        {
            if (x == 0 && y == 0)
                return 0;
            return atan2(y, x) - HOA_PI2;
        }
        
        //! The radius converter function.
        /** This function takes an abscissa and an ordinate and convert them to a radius.
         @param     x		The abscissa.
         @param     y		The ordinate.
         @param     z		The height.
         @return    The radius.
         */
        static inline T radius(const T x, const T y, const T z = 0.)
        {
            return sqrt(x*x + y*y + z*z);
        }
        
        //! The elevation converter function.
        /** This function takes an abscissa and an ordinate and convert them to a elevation.
         @param     x		The abscissa.
         @param     y		The ordinate.
         @param     z		The height.
         @return    The elevation.
         */
        static inline double elevation(const T x, const T y, const T z = 0.)
        {
            if(z == 0)
                return 0;
            return asin(z / sqrt(x*x + y*y + z*z));
        }
        
        //! The factorial
        /** The function computes the factorial, the product of all positive integers less than or equal to an integer.
         \f[n! = n \times (n - 1) \times (n - 2) \times {...} \f]
         @param     n     The interger.
         @return    The function return the factorial of n.
         */
        static inline long double factorial(long n)
        {
            long double result = n;
            if(n == 0)
                return 1;
            
            while(--n > 0)
                result *= n;
            
            return result;
        }
        
        //! The great-circle distance
        /** This function compute the great circle distance of two points in radians on a unit sphere.
         @param     azimuth1		The azimuth of the first point in radian.
         @param     elevation1   The elevation of the first point in radian.
         @param     azimuth2		The azimuth of the second point in radian.
         @param     elevation2   The elevation of the second point in radian.
         @return    The great-circle distance.
         */
        static inline T distance_spherical(const T azimuth1, const T elevation1, const T azimuth2, const T elevation2)
        {
            return acos(sin(elevation1) * sin(elevation2) + cos(elevation1) * cos(elevation2) * cos(azimuth1 - azimuth2));
        }
    };
}

#endif
