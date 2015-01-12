/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_DEFS_LIGHT
#define DEF_HOA_DEFS_LIGHT

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#elif _WINDOWS
#include <gsl_cblas.h>
#else
#include <cblas.h>
#endif

#include <stdio.h>
#include <stdarg.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#ifdef PD_DEBUG
#include "../../CicmWrapper/Sources/cicm_wrapper.h"
#endif
#ifdef MAX_DEBUG
#include <ext.h>
#include <ext_obex.h>
#include <z_dsp.h>
#include <ext_common.h>
#endif


#define HOA_PI  3.14159265358979323846264338327950288
#define HOA_2PI 6.283185307179586476925286766559005
#define HOA_PI2 1.57079632679489661923132169163975144
#define HOA_PI4 0.785398163397448309615660845819875721
#define HOA_1DEG (HOA_PI / 180.)

using namespace std;

namespace hoa
{
    typedef unsigned long ulong;
#ifdef _WINDOWS
    static inline double round(double val)
    {
        return floor(val + 0.5);
    }
#endif
    
    enum Dimension
    {
        Hoa2d = 0,
        Hoa3d = 1
    };
}

#endif


