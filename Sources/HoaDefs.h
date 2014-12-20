/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_DEFS__
#define __DEF_HOA_DEFS__

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
#include <algorithm>

#ifdef __APPLE__

#include <Accelerate/Accelerate.h>

#endif

#ifdef _WINDOWS

#include <gsl_cblas.h>

#endif

#ifdef _LINUX
#include <cblas.h>
#endif

#ifdef PD_DEBUG
#include "../ThirdParty/PureData/Sources/cicm_wrapper.h"
#endif
#ifdef MAX_DEBUG
#include <ext.h>
#include <ext_obex.h>
#include <z_dsp.h>
#include <ext_common.h>
#endif

#define HOA_PI  (3.141592653589793238462643383279502884)
#define HOA_2PI (6.283185307179586476925286766559005)
#define HOA_PI2 (1.57079632679489661923132169163975144)
#define HOA_PI4 (0.785398163397448309615660845819875721)

#define NUMBEROFLINEARPOINTS 10000
#define NUMBEROFCIRCLEPOINTS 36000
#define NUMBEROFCIRCLEPOINTS_UI 360
#define NUMBEROFCIRCLEPOINTS_UI2 180

namespace Hoa
{
};

#endif


