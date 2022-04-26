#include "packed.h"
#include "misc.h"
#include "v1.h"
#include "v1_tablegen.cpp"
#include "v1_basecases.cpp"
#include "v1_fft.cpp"
#include "v1_mul.cpp"
template class pd_fft_ctx<4>;   // in case we forgot to emit the code
