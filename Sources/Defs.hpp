/*
// Copyright (c) 2012-2015 Eliott Paris & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_DEFS_LIGHT
#define DEF_HOA_DEFS_LIGHT

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <malloc.h>
#endif

#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

// avoid min and max to be defined when compiling on VS
#ifdef _MSC_VER
#define NOMINMAX
#endif

#ifdef PD_DEBUG
#include "../../CicmWrapper/Sources/cicm_wrapper.h"
#endif
#ifdef MAX_DEBUG
#include <ext.h>
#include <ext_obex.h>
#include <z_dsp.h>
#include <ext_common.h>
#endif

#if (__cplusplus <= 199711L)
#define noexcept
#define nullptr NULL
#define override
#endif

#define HOA_PI  3.14159265358979323846264338327950288
#define HOA_2PI 6.283185307179586476925286766559005
#define HOA_PI2 1.57079632679489661923132169163975144
#define HOA_PI4 0.785398163397448309615660845819875721
#define HOA_EPSILON 1e-6

//! The namespace of the hoa library.
/** All the classes of the hoa library are inside this namespace.
 */
namespace hoa
{
#if (__cplusplus <= 199711L)
    inline std::string to_string(const size_t& val)
    {
        char buffer[1024];
		sprintf(buffer, "%lu", val);
        return buffer;
    }

    inline std::string to_string(const long& val)
    {
        char buffer[1024];
        sprintf(buffer, "%ld", val);
        return buffer;
    }

    inline std::string to_string(const float& val)
    {
        char buffer[1024];
        sprintf(buffer, "%f", val);
        return buffer;
    }

    inline std::string to_string(const double& val)
    {
        char buffer[1024];
        sprintf(buffer, "%lf", val);
        return buffer;
    }

#endif

    //! The dimension of class.
    /** Most of the classes are specialized for 2d or 3d.
     */
    enum Dimension
    {
        Hoa2d = 0, /*!<  The 2d dimension. */
        Hoa3d = 1  /*!<  The 3d dimension. */
    };
}

#endif


