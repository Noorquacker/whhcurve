#define DATAFITTING_INCLUDE_ALL

#ifdef DATAFITTING_INCLUDE_ALL
#include "tntnr.h"
#include "svdfit.h"
#include "gaussj.h"
#include "mrqmin.h"
#endif

const char* SVDFIT_INDEX_OUT_OF_BOUNDS = "Fit index out of bounds.";
const char* SVDFIT_TOO_MANY_FUNCTIONS = "Too many functions for declared SVDFit.";
const char* SVDFIT_TOO_MANY_POINTS = "Too many points for declared SVDFit.";
const char* SVDFIT_DIM_MISMATCH = "Dimension mismatch in SVDFit.";

const char* GAUSSJ_DIM_MISMATCH = "Dimension mismatch in gaussj.";
const char* GAUSSJ_SINGULAR_1 = "Singular matrix in gaussj (#1).";
const char* GAUSSJ_SINGULAR_2 = "Singular matrix in gaussj (#2).";

const char* COVSRT_DIM_MISMATCH = "Dimension mismatch in covsrt.";
const char* MRQCOF_DIM_MISMATCH = "Dimension mismatch in mrqcof.";
const char* FGAUSS_DIM_MISMATCH = "Dimension mismatch in fgauss.";

const char* MRQMIN_DIM_MISMATCH1 = "Dimension mismatch in mrqmin (x/y/sig).";
const char* MRQMIN_DIM_MISMATCH2 = "Dimension mismatch in mrqmin (a/covar/alpha).";
const char* MRQMIN_BAD_STATE_START = "LevenbergMarquardtFit is in the wrong state to start.";
const char* MRQMIN_BAD_STATE_STEP = "LevenbergMarquardtFit is in the wrong state to step.";
const char* MRQMIN_BAD_STATE_STOP = "LevenbergMarquardtFit is in the wrong state to stop.";
const char* MRQMIN_BAD_STATE_READ = "LevenbergMarquardtFit is in the wrong state to read parameters.";
