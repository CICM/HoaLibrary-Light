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

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <malloc.h>
#endif

#include <string>
#include <cmath>
#include <cstring>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>

// avoid min and max to be defined when compiling on VS
#ifdef _MSC_VER
#define NOMINMAX
#endif

#include <type_traits>
#define hoa_unused(expr) (void)(expr)

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
    //! The dimension of class.
    //! Most of the classes are specialized for 2d or 3d.
    enum Dimension
    {
        Hoa2d = 0, /*!<  The 2d dimension. */
        Hoa3d = 1  /*!<  The 3d dimension. */
    };
}
