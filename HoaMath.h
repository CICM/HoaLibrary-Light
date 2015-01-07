/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_MATHS__
#define __DEF_HOA_MATHS__

#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string>
#include <assert.h>
#include <string.h>

#ifdef _WINDOWS
static inline double round(double val)
{    
    return floor(val + 0.5);
}

#undef min
#undef max

#endif

//! The high order ambisonic namespace.
/**
 This namespace have all the standard methods and functions necessary for ambisonic processing.
*/
namespace Hoa
{
	//! The minimum function
    /** This function returns the minimum value between two values

	 @param     v1   The first value.
	 @param     v2   The low boundary.
	 @return    The minimum value.

	 @see    max
     */
	inline double min(const double v1, const double v2)
	{
		if (v1 <= v2)
			return v1; 
		else
			return v2;
	}

	//! The maximum function
    /** This function returns the maximum value between two values

	 @param     v1   The first value.
	 @param     v2   The low boundary.
	 @return    The maximum value.

	 @see    min
     */
	inline double max(const double v1, const double v2)
	{
		if (v1 >= v2)
			return v1;
		else
			return v2;
	}

	//! The clipping function
    /** The function clips a number between boundaries. \n
	 If \f$x < min\f$, \f$f(x) = min\f$ else if \f$x > max\f$, \f$f(x) = max\f$ else \f$f(x) = x\f$.

	 @param     value   The value to clip.
	 @param     min     The low boundary.
	 @param     max     The high boundary.
	 @return    The function return the clipped value.

	 @see    clip_min
	 @see    clip_max
     */
    inline double clip_minmax(const double value, const double min, const double max)
    {
        if(value < min)
            return min;
        else if(value > max)
            return max;
        else
            return value;
    }

    //! The minimum clipping function
    /** The function clips a number at a minimum value. \n
	 If \f$x < min\f$, \f$f(x) = min\f$ else \f$f(x) = x\f$.

	 @param     value   The value to clip.
	 @param     low     The low boundary.
	 @return    The function return the clipped value.

	 @see    clip_minmax
	 @see    clip_max
     */
    inline double clip_min(const double value, const double low)
    {
        if(value < low)
            return low;
        else
            return value;
    }

    //! The maximum clipping function
    /** The function clips a number at a maximum value. \n
	 If \f$x > max\f$, \f$f(x) = max\f$ else \f$f(x) = x\f$.

	 @param     value   The value to clip.
	 @param     high    The high boundary.
	 @return    The function return the clipped value.

        @see    clip_minmax
        @see    clip_min
     */
    inline double clip_max(const double value, const double high)
    {
        if(value > high)
            return high;
        else
            return value;
    }

	//! The factorial
    /** The function computes the factorial, the product of all positive integers less than or equal to an integer.
        \f[n! = n \times (n - 1) \times (n - 2) \times {...} \f]

        @param     n     The interger.
        @return    The function return the factorial of n.

        @see    double_factorial
     */
    inline long double factorial(long n)
	{
        long double result = n;
		if(n == 0)
			return 1;

		while(--n > 0)
			result *= n;

		return result;
	}

    //! The double factorial
    /**	The function computes the double factorial, the product of all the odd integers up to some odd positive integer :\n
	 \f[n!! = n \times (n - 2) \times (n - 4) \times {...} \f]

	 @param     n     The interger.
	 @return    The function return the double factorial of n.

	 @see    factorial
     */
    inline long double double_factorial(long n)
	{
		if (n == 0 || n == -1) {
			return 1;
		}

		long result = n;
		while ((n -= 2) > 0) {
			result *= n;
		}

		return result;
	}

    //! The associated Legendre polynomials
    /**	The function computes the associated Legendre polynomial \f$P(l, m)\f$ that is a part of the formula that compute the spherical harmonic coefficient where l is the band and the m is the argument of a spherical harmonic and x is the cosinus of the elevation. It uses three recurrence formulas :
        \f[P(l, l)(x) = (-1)^l \times (2l - 1)!! \times (1 - x^2)^{0.5l}\f]
        \f[P(l + 1, l)(x) = x \times (2l + 1) \times P(l, l)\f]
        \f[P(l + 1, m)(x) = \frac{(2l + 1) \times x \times P(m, l) - (l + m) \times P(m, l - 1)}{(l - m + 1)}\f]
        with \f$0 \leq l\f$ and \f$-l \leq m \leq +l\f$

	 @param     l   The degree of the spherical harmonic.
	 @param     m	The order of the spherical harmonic.
	 @param     x    The cosinus of the elevation.
	 @return    The function return the associated Legendre polynomial of x for l and m.

	 @see    legendre_normalization
	 @see    spherical_harmonics
     */
	inline double associated_legendre(int l, int m, const double x)
	{
        l = abs(l);
        m = abs(m);

		if(l == m)
        {
			return pow((double)-1.0, (double)m) * pow((double)(1. - x * x), (double)(0.5 * m)) * double_factorial(2. * m - 1);
		}
		else if(l == m + 1)
        {
			return x * associated_legendre(m, m, x) * (2 * m + 1);
		}
        else
        {
            return ((double)(2 * l - 1) * x *  associated_legendre(l - 1, m, x) - (double)(l + m - 1.) * associated_legendre(l - 2, m, x)) / (double)(l - m);
        }
	}

    //! The legendre normalization
    /**	The function normalizes the associated Legendre polynomial over \f$2\pi\f$ : \n
	 if \f$ m = 0\f$
	 \f[K(l) = \sqrt{\frac{2 \times l + 1}{4\pi}}\f]
	 else
	 \f[K(l, m) = \sqrt{2 \times \frac{2 \times l + 1}{4\pi} \times \frac{(l - |m|)!}{(l + |m|)!}}\f]
	 with \f$0 \leq l\f$ and \f$-l \leq m \leq +l\f$

	 @param     l    The band of the spherical harmonic.
	 @param     m    The argument of the spherical harmonic.
	 @return    The function return the normalization of the associated Legendre polynomial for l and m.

	 @see    associated_legendre
	 @see    spherical_harmonics
     */
    inline double legendre_normalization(const int l, const int m)
	{
        if(m == 0)
            return sqrt((2. * l + 1.) / (4. * HOA_PI));
        else
            return sqrt((2. * l + 1.) / (4. * HOA_PI) * (long double)factorial(l - abs(m)) / (long double)factorial(l + abs(m))) * sqrt(2.);
	}

    /*
	 inline double hoa_normalization(int l, int m)
	 {
	 l = abs(l);
	 m = abs(m);
	 if(m == l || m == 0)
	 return 1. / double_factorial(2. * m - 1);
	 else if(l == m + 1)
	 return 1. / associated_legendre(l, m , cos(HOA_PI / 4.));
	 else if(l == m + 2)
	 return 1. / associated_legendre(l, m , cos(HOA_PI / 2.));
	 else
	 return 1.;
	 }*/

    //! The azimuth part of the spherical harmonics function
    /**	The function computes the azimuth coefficient of the spherical harmonic \f$[l, m]\f$ for an angle \f$\phi\f$ in radian :\n
	 if \f$ m \geq 0\f$
	 \f[Y_{azimuth}(l, m, \phi) = cos(m \times \phi)\f]
	 else
	 \f[Y_{azimuth}(l, m, \phi) = sin(-m \times \phi)\f]
	 with \f$0 \leq l\f$ and \f$-l \leq m \leq +l\f$

	 @param     l    The band of the spherical harmonic.
	 @param     m    The argument of the spherical harmonic.
	 @param     phi  The azimuth.
	 @return    The function return the azimuth coefficient for phi of the spherical harmonic of band l and argument m.

	 @see    spherical_harmonics_elevation
	 @see    spherical_harmonics
     */
    inline double spherical_harmonics_azimuth(const int l, const int m, const double phi)
	{
        if(m >= 0)
            return cos((double)m * phi);
        else
            return sin((double)-m * phi);
    }

    //! The elevation part of the spherical harmonics function
    /**	The function computes the elevation coefficient of the the spherical harmonic \f$[l, m]\f$  for an angle \f$\theta\f$ in radian. It uses the associated Legendre polynomial and applies the Legendre normalization :
	 \f[Y_{elevation}(l, m, \theta) = N(l, m) \times P(l, m, cos(/theta)\f]
	 with \f$0 \leq l\f$, \f$-l \leq m \leq +l\f$ and \f$N\f$ the normalization and \f$P\f$ the associated Legendre polynomial.

	 @param     l        The band of the spherical harmonic.
	 @param     m        The argument of the spherical harmonic.
	 @param     theta    The elevation.
	 @return    The function return the elevation coefficient for theta of the spherical harmonic of band l and argument m.

	 @see    legendre_normalization
	 @see    associated_legendre
	 @see    spherical_harmonics_azimuth
	 @see    spherical_harmonics
     */
    inline double spherical_harmonics_elevation(const int l, const int m, const double theta)
	{
        return associated_legendre(l, m, cos(theta)) * legendre_normalization(l, m);
    }

    //! The spherical harmonics function
    /** The function computes the spherical harmonics coefficient for the angles \f$\phi\f$ and \f$\theta\f$ in radian.
	 \f[Y(l, m, \phi, \theta) = Y_{azimuth}(l, m, \phi) \times Y_{elevation}(l, m, \theta)\f]

	 @param     l        The band of the spherical harmonic.
	 @param     m        The argument of the spherical harmonic.
	 @param     phi      The azimuth.
	 @param     theta    The elevation.
	 @return    The function return the coefficient for phi and theta for the spherical harmonics of band l and argument m.

	 @see    spherical_harmonics_azimuth
	 @see    spherical_harmonics_elevation
     */
    inline double spherical_harmonics(const int l, const int m, const double phi, const double theta)
	{
        return spherical_harmonics_elevation(l, m, theta) * spherical_harmonics_azimuth(l, m, phi);
    }

    //! The wrapping function
    /** The function wraps a number between boundarys.

	 @param     value   The value to wrap.
	 @param     low     The low boundary.
	 @param     high    The high boundary.
	 @return    The function return the wrapped value.

	 @see    wrap_twopi
     */
    inline double wrap(double value, const double low, const double high)
    {
        double increment = high - low;
        while(value < low)
        {
            value += increment;
        }

        while(value > high)
        {
            value -= increment;
        }

        return value;
    }

    //! The wrapping function over \f$2\pi\f$
    /** The function wraps a number between \f$0\f$ and \f$2\pi\f$.

	 @param     value   The value to wrap.
	 @return    The function return the wrapped value.

	 @see    wrap, wrap_360
     */
    inline double wrap_twopi(const double value)
    {
        double new_value = value;
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

	//! The wrapping function in degrees.
    /** The function wraps a number between \f$0\f$ and \f$360°\f$.

	 @param     value   The value to wrap.
	 @return    The wrapped value.

	 @see    wrap, wrap_twopi
     */
	inline double wrap_360(const double value)
    {
		double new_value = value;
        while(new_value < 0.)
            new_value += 360.;
        while(new_value >= 360.)
            new_value -= 360.;

        return new_value;
    }

	//! The abscissa converter function.
    /** This function takes a radius and an azimuth value and convert them to an abscissa value.

	 @param     radius		The radius value (greather than 0).
	 @param     azimuth		The azimuth value (between \f$0\f$ and \f$2\pi\f$).
	 @return    The abscissa value.

	 @see    ordinate
     */
    inline double abscissa(const double radius, const double azimuth, const double elevation = 0.)
	{
		return radius * cos(azimuth + HOA_PI2) * cos(elevation);
	}
	
	//! The ordinate converter function.
    /** This function takes a radius and an azimuth value and convert them to an ordinate value.
	 
	 @param     radius		The radius value (greather than 0).
	 @param     azimuth		The azimuth value (between \f$0\f$ and \f$2\pi\f$).
	 @return    The ordinate value.
	 
	 @see    abscissa
     */
    inline double ordinate(const double radius, const double azimuth, const double elevation = 0.)
	{
		return radius * sin(azimuth + HOA_PI2) * cos(elevation);
	}
    
    inline double height(const double radius, const double azimuth, const double elevation)
	{
		return radius * sin(elevation);
	}
	
	//! The azimuth converter function.
    /** This function takes an abscissa and an ordinate value and convert them to an azimuth value.
     
	 @param     x		The abscissa value.
	 @param     y		The ordinate value.
	 @return    The azimuth value.
     
	 @see    radius
     */
	inline double azimuth(const double x, const double y, const double z = 0.)
	{
		if (x == 0 && y == 0)
			return 0;
		return atan2(y, x) - HOA_PI2;
	}
	
	//! The radius converter function.
    /** This function takes an abscissa and an ordinate value and convert them to a radius value.
	 
	 @param     x		The abscissa value.
	 @param     y		The ordinate value.
	 @return    The radius value.
	 
	 @see    azimuth
     */
    inline double radius(const double x, const double y, const double z = 0.)
	{
		return sqrt(x*x + y*y + z*z);
	}
    
    inline double elevation(const double x, const double y, const double z)
	{
		if(z == 0)
			return 0;
		return asin(z / sqrt(x*x + y*y + z*z));
	}
    
    inline double distance(double x1, double y1, double x2, double y2)
	{
		return sqrt((x1-x2) * (x1-x2) + (y1-y2) * (y1-y2));
	}
    
    inline double distance(double x1, double y1, double z1, double x2, double y2, double z2)
	{
		return sqrt((x1-x2) * (x1-x2) + (y1-y2) * (y1-y2) + (z1-z2) * (z1-z2));
	}
    
	inline double distance_radian(double angle1, double angle2)
    {
        angle1 = wrap_twopi(angle1);
        angle2 = wrap_twopi(angle2);
        if(angle1 > angle2)
            return (angle1 - angle2);
        else
            return (angle2 - angle1);
    }
    
    //! The great-circle distance
    /** This function compute the great circle distance of two points in radians on a unit sphere.
     
        @param     azimuth1		The azimuth of the first point in radian.
        @param     elevation1   The elevation of the first point in radian.
        @param     azimuth2		The azimuth of the second point in radian.
        @param     elevation2   The elevation of the second point in radian.
        @return    The great-circle distance.
     
        @see    radius
     */
    inline double distance_spherical(const double azimuth1, const double elevation1, const double azimuth2, const double elevation2)
    {
        return acos(sin(elevation1) * sin(elevation2) + cos(elevation1) * cos(elevation2) * cos(azimuth1 - azimuth2));
    }

	inline double spherical_azimuth_interpolation(const double azimuth1, const double elevation1, const double azimuth2, const double elevation2, double mu)
	{
		double distance;
		double angle1 = wrap_twopi(azimuth1);
        double angle2 = wrap_twopi(azimuth2);
        if(angle1 > angle2)
            distance = (angle1 - angle2);
        else
            distance = (angle2 - angle1);
		if(azimuth1 == azimuth2)
			return azimuth1;
		if(azimuth1 < azimuth2)
		{
			if(distance > HOA_PI)
				return wrap_twopi(azimuth1 - (HOA_2PI - distance) * mu);
			else
				return azimuth1 + distance * mu;
		}
		else
		{
			if(distance > HOA_PI)
				return wrap_twopi(azimuth1 + (HOA_2PI - distance) * mu);
			else
				return azimuth1 - distance * mu;
		}
	}

	inline double spherical_elevation_interpolation(const double azimuth1, const double elevation1, const double azimuth2, const double elevation2, double mu)
	{
		double distance = distance_radian(elevation1, elevation2);
		if(round(azimuth1 / HOA_2PI *360.) == round(wrap_twopi(azimuth2 + HOA_PI) / HOA_2PI * 360.))
			distance = elevation1 + elevation2;

		if(distance > HOA_PI)
			distance  = HOA_2PI - distance;
		if(elevation1 < elevation2)
			return elevation1 + distance * mu;
		else
			return elevation2 + distance * (1. - mu);
	}

	inline double radianClosestDistance(double angle1, double angle2)
    {
        double minRad, maxRad;
        angle1 = wrap_twopi(angle1);
        angle2 = wrap_twopi(angle2);
        minRad = min(angle1, angle2);
        maxRad = max(angle1, angle2);

        if (maxRad - minRad <= HOA_PI)
            return maxRad - minRad;
        else
            return HOA_2PI - (maxRad - minRad);
    }

	inline double radToDeg(const double radian)
    {
        return radian * 360. / HOA_2PI;
    }

    inline double degToRad(const double degree)
    {
        return degree / 360. * HOA_2PI;
    }

	inline long mstosamps(const double ms, const double samplerate = 44100.)
    {
        return (long)(samplerate * ms * 0.001);
    }

    inline double sampstoms(const double s, const double samplerate=44100.)
    {
        return 1000. * s / samplerate;
    }

    inline double atodb(const double amp)
    {
        return (amp <= 0.) ? -999.f : (20. * log10(amp));
    }

    inline double dbtoa(const double dB)
    {
        return pow(10., dB * 0.05);
    }

	inline double safediv(const double num, const double denom)
    {
        return denom == 0. ? 0. : num/denom;
    }

	inline double scale(const double in, const double inlow, const double inhigh, const double outlow, const double outhigh, const double power)
    {
        double value;
        double inscale = safediv(1., inhigh - inlow);
        double outdiff = outhigh - outlow;

        value = (in - inlow) * inscale;
        if (value > 0.0)
            value = pow(value, power);
        else if (value < 0.0)
            value = -pow(-value, power);
        value = (value * outdiff) + outlow;

        return value;
    }

    inline double scale(const double in, const double inlow, const double inhigh, const double outlow, const double outhigh)
    {
        return ( (in - inlow) * safediv(1., inhigh - inlow) * (outhigh - outlow) ) + outlow;
    }

	inline bool isInside(const double val, const double v1, const double v2)
	{
        return (v1 <= v2) ? (val >= v1 && val <= v2) : (val >= v2 && val <= v1);
	}

	inline bool isInsideRad(const double radian, const double loRad, const double hiRad)
	{
        return isInside(wrap_twopi(radian-loRad), 0, wrap_twopi(hiRad-loRad));
	}

    inline bool isInsideDeg(const double degree, const double loDeg, const double hiDeg)
	{
        return isInside(wrap_360(degree-loDeg), 0, wrap_360(hiDeg-loDeg));
	}

	inline void vector_add(unsigned int size, double *vec, double inc)
	{
		for (unsigned int i=0; i < size; i++)
			vec[i] += inc;
	}

	inline void vector_clip_minmax(unsigned int size, double *vec, double min, double max)
	{
		for (unsigned i=0; i < size; i++)
			vec[i] = clip_minmax(vec[i], min, max);
	}

    inline void vector_sort(unsigned int size, double* vector)
	{
        int index;
        double* temp = new double[size];
        if(temp && size)
        {
            index  = 0;
            temp[0] = vector[0];
            temp[size - 1] = vector[0];

            for(unsigned int i = 1; i < size; i++)
            {
                if(vector[i] < temp[0]) // Get the minimum
                {
                    temp[0] = vector[i];
                    index = i;
                }
                if(vector[i] > temp[size - 1]) // Get the maximum
                    temp[size - 1] = vector[i];
            }
            vector[index] -= 1;

            for(unsigned int i = 1; i < size - 1; i++)
            {
                temp[i] = temp[size - 1];
                index   = -1;
                for(unsigned int j = 0; j < size; j++)
                {
                    if(vector[j] >= temp[i-1] && vector[j] <= temp[i])
                    {
                        temp[i] = vector[j];
                        index = j;
                    }
                }
                if(index > -1)
                {
                    vector[index] -= 1;
                }
            }
            memcpy(vector, temp, size * sizeof(double));
            delete [] temp;
        }
	}

	inline void vector_sort_coordinates(unsigned int size, double* azimuths, double* elevations, double azymuth, double elevation)
	{
        double* abs	= new double[size];
		double* ord	= new double[size];
		double* azi	= new double[size];
		double* cpa	= new double[size];
		double* cpe	= new double[size];
		int* index	= new int[size];
		double g_x = abscissa(1, azymuth, elevation), g_y = ordinate(1, azymuth, elevation);
		memcpy(cpa, azimuths, size * sizeof(double));
		memcpy(cpe, elevations, size * sizeof(double));
	
		for(unsigned int i = 0; i < size; i++)
		{
			abs[i] = abscissa(1., azimuths[i], elevations[i]);
			ord[i] = ordinate(1., azimuths[i], elevations[i]) - g_y;
			if(elevation >= 0 && elevations[i] < 0.)
			{
				if(abs[i] < 0)
					abs[i] *= -1;
				abs[i] += 2.;
			}
			else if(elevation < 0 && elevations[i] > 0.)
			{
				if(abs[i] > 0)
					abs[i] *= -1;
				abs[i] += 2.;
			}
			abs[i] -= g_x;
			
		}
		double max = -100;
		int inc = size - 1;
		for(unsigned int i = 0; i < size; i++)
		{
			azi[i] = wrap_twopi(azimuth(abs[i], ord[i]));
			if(azi[i] > max)
			{
				max = azi[i];
				index[inc] = i;
			}
			
		}
		
		azi[index[inc]] = -100;
		inc--;
		while(inc > -1)
		{
			max = -1;
			for(unsigned int j = 0; j < size; j++)
			{
				if(azi[j] > max)
				{
					max = azi[j];
					index[inc] = j;
				}
			}
			azi[index[inc]] = -100;
			inc--;
		}

		for(unsigned int i = 0; i < size; i++)
		{
			azimuths[i] = cpa[index[i]];
			elevations[i] = cpe[index[i]];
		}
		
		delete [] abs;
		delete [] ord;
		delete [] azi;
		delete [] cpa;
		delete [] cpe;
		delete [] index;

	}

	inline void vector_sort_byone(unsigned int size, double* vector, double* vector2)
	{
        int index;
        double* temp	= new double[size];
		double* temp2	= new double[size];
        if(temp && size)
        {
            index  = 0;
            temp[0] = vector[0];
            temp[size - 1] = vector[0];
			temp2[0] = vector2[0];
			temp2[size - 1] = vector2[0];
            for(unsigned int i = 1; i < size; i++)
            {
                if(vector[i] < temp[0]) // Get the minimum
                {
                    temp[0] = vector[i];
					temp2[0] = vector2[i];
                    index = i;
                }
                if(vector[i] > temp[size - 1]) // Get the maximum
				{
                    temp[size - 1] = vector[i];
					temp2[size - 1] = vector2[i];
				}
            }
            vector[index] -= 1;

            for(unsigned int i = 1; i < size - 1; i++)
            {
                temp[i] = temp[size - 1];
				temp2[i] = temp2[size - 1];
                index   = -1;
                for(unsigned int j = 0; j < size; j++)
                {
                    if(vector[j] >= temp[i-1] && vector[j] <= temp[i])
                    {
                        temp[i] = vector[j];
						temp2[i] = vector2[j];
                        index = j;
                    }
                }
                if(index > -1)
                {
                    vector[index] -= 1;
                }
            }
            memcpy(vector, temp, size * sizeof(double));
			memcpy(vector2, temp2, size * sizeof(double));
            delete [] temp;
			delete [] temp2;
        }
	}
}

#endif
