#ifndef INCLUDEFILE_H
#define INCLUDEFILE_H


#include <malloc.h>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

#include "nmsimplex.h"


double R_G = 8.31446261815324;
double T_0 = 273;

#include "TEMP_DEPENDENCE.h"
#include "JENSEN_SEATON.h"
#include "REDLICH_PETERSON.h"
#include "LANGMUIR.h"
#include "FREUNDLICH.h"
#include "SIPS.h"
#include "TOTH.h"
#include "UNILAN.h"
#include "KELLER_STAUDT_TOTH.h"
#include "ERRORS.h"


#endif