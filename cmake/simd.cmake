include(CheckIncludeFile)
include(CheckCSourceRuns)

include("${CMAKE_CURRENT_LIST_DIR}/CompareVersionStrings.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

set(HAVE_CHECK_COMPILER_FLAGS 0)
COMPARE_VERSION_STRINGS("${CMAKE_VERSION}" "3.0.2" HAVE_CHECK_COMPILER_FLAGS)

# KitWare suport test for presence of compiler flags (Really slow)
if(HAVE_CHECK_COMPILER_FLAGS LESS 0)
  # TODO: Establish flags needed for testing availability of instructions
  #       based on platform
else()
  if (0)
    check_c_compiler_flag(-msse3  HAS_MSSE3)
    check_c_compiler_flag(-msse4  HAS_MSSE4)
    check_c_compiler_flag(-msse4a HAS_MSSE4a)
    check_c_compiler_flag(-maes   HAS_MAES)
    check_c_compiler_flag(-mavx   HAS_MAVX)
    check_c_compiler_flag(-mavx2  HAS_MAVX2)
    check_c_compiler_flag(-mfma   HAS_MFMA)
    check_c_compiler_flag(/arch:SSE2 HAS_ARCH_SSE2)
    check_c_compiler_flag(/arch:AVX  HAS_ARCH_AVX)
    check_c_compiler_flag(/arch:AVX2  HAS_ARCH_AVX2)
    # TODO: Establish flags needed for testing availability of instructions
  else()
    # TODO: Establish flags needed for testing availability of instructions
  endif()
endif()

if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|^i[3,6,9]86$")
  message("System processor is ${CMAKE_SYSTEM_PROCESSOR}")

  # SSE2
  check_include_file(emmintrin.h HAVE_EMMINTRIN_H)

  # SSE3
  check_include_file(pmmintrin.h HAVE_PMMINTRIN_H "-msse3")
  
  # SSE4.1
  check_include_file(smmintrin.h HAVE_SMMINTRIN_H "-msse4")
  
  # SSE4.2
  check_include_file(nmmintrin.h HAVE_NMMINTRIN_H "-msse4")
  
  # SSE4A - AMD only XTRQ/INSERTQ (Combined mask-shift), MOVNTSD/MOVNTSS (Scalar stream/store)
  check_include_file(ammintrin.h HAVE_AMMINTRIN_H "-msse4a")
  
  # AES
  check_include_file(wmmintrin.h HAVE_WMMINTRIN_H "-maes")
  
  # AVX
  check_include_file(immintrin.h HAVE_IMMINTRIN_H "-mavx")
  
  # AVX2
  if (HAVE_IMMINTRIN_H)
    set(CMAKE_REQUIRED_FLAGS "-mavx2")
    check_c_source_runs("
    #include <immintrin.h>
    int main()
    {
      __m128 a = _mm_setzero_ps();
      __m128 b = _mm_broadcastss_ps(a);
      return 0;
    }"
    HAVE_ZMMINTRIN_H)
  endif()

  # FMA (cannot be included directly)
  set(CMAKE_REQUIRED_FLAGS "-mfma")
  check_c_source_runs("
      #include <immintrin.h>
      int main()
      {
        
        __m128 a = _mm_setzero_ps();
        __m128 b = _mm_setzero_ps();
        __m128 c = _mm_setzero_ps();
        __m128 d = _mm_fmadd_ps(a,b,c);
        
        return 0;
      }"
      HAVE_FMAINTRIN_H)

  # Hack for CYGWIN at work
  if (CYGWIN)
    set(HAVE_ZMMINTRIN_H 1)
    set(HAVE_FMAINTRIN_H 1)
  endif()
    
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64")
  # SSE2
  check_include_file(emmintrin.h HAVE_EMMINTRIN_H)
  
  # SSE3
  check_include_file(pmmintrin.h HAVE_PMMINTRIN_H "/arch:SSE2")
  
  # SSE4.1
  check_include_file(smmintrin.h HAVE_SMMINTRIN_H "/arch:SSE2")
  
  # SSE4.2
  check_include_file(nmmintrin.h HAVE_NMMINTRIN_H "/arch:SSE2")
  
  # SSE4A - AMD only XTRQ/INSERTQ (Combined mask-shift), MOVNTSD/MOVNTSS (Scalar stream/store)
  check_include_file(ammintrin.h HAVE_AMMINTRIN_H "/arch:SSE2")
  
  # AES
  check_include_file(wmmintrin.h HAVE_WMMINTRIN_H "/arch:SSE2")
  
  # AVX
  check_include_file(immintrin.h HAVE_IMMINTRIN_H "/arch:AVX")
  
  # Intel doesn't allow inclusion of zmmintrin.h
  if (HAVE_IMMINTRIN_H)
    # Make example using _mm_broadcastss_ps
    set(CMAKE_REQUIRED_FLAGS "/arch:AVX2")
    check_c_source_runs("
    #include <immintrin.h>
    int main()
    {
      __m128 a = _mm_setzero_ps();
      __m128 b = _mm_broadcastss_ps(a);
      return 0;
    }"
    HAVE_ZMMINTRIN_H)
  endif()

  # Microsoft screws up FMA instructions
  if (1)
    # Intel doesn't allow inclusion of fmaintrin.h
    if(HAVE_IMMINTRIN_H)
      # Make example using _mm_fmadd_ps
      set(CMAKE_REQUIRED_FLAGS "/arch:AVX2")
      check_c_source_runs("
      #include <immintrin.h>
      int main()
      {
        __m128 a = _mm_setzero_ps();
        __m128 b = _mm_setzero_ps();
        __m128 c = _mm_setzero_ps();
        __m128 d = _mm_fmadd_ps(a,b,c);
        return 0;
      }"
      HAVE_FMAINTRIN_H)
    endif()
  endif()
else()
  message(FATAL_ERROR "Unsupported CMake system processor: ${CMAKE_SYSTEM_PROCESSOR}")
endif()

if ((CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64") AND HAVE_EMMINTRIN_H)
  set(ENABLE_SSE2_DEFAULT ON)
else()
  set(ENABLE_SSE2_DEFAULT OFF)
endif()

if ((CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64") AND HAVE_PMMINTRIN_H)
  set(ENABLE_SSE3_DEFAULT ON)
else()
  set(ENABLE_SSE3_DEFAULT OFF)
endif()

if ((CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64") AND HAVE_SMMINTRIN_H AND HAVE_NMMINTRIN_H)
  set(ENABLE_SSE4_DEFAULT ON)
else()
  set(ENABLE_SSE4_DEFAULT OFF)
endif()

if ((CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64") AND HAVE_IMMINTRIN_H)
  set(ENABLE_AVX_DEFAULT ON)
else()
  set(ENABLE_AVX_DEFAULT OFF)
endif()

if ((CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64") AND HAVE_ZMMINTRIN_H)
  set(ENABLE_AVX2_DEFAULT ON)
else()
  set(ENABLE_AVX2_DEFAULT OFF)
endif()

if ((CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64") AND HAVE_FMAINTRIN_H)
  set(ENABLE_FMA_DEFAULT ON)
else()
  set(ENABLE_FMA_DEFAULT OFF)
endif()

option(ENABLE_SSE2 "Enable compile-time SSE2 support." ${ENABLE_SSE2_DEFAULT})
option(ENABLE_SSE3 "Enable compile-time SSE3 support." ${ENABLE_SSE3_DEFAULT})
option(ENABLE_SSE4 "Enable compile-time SSE4 support." ${ENABLE_SSE4_DEFAULT})
option(ENABLE_FMA  "Enable compile-time FMA support."  ${ENABLE_FMA_DEFAULT})
option(ENABLE_AVX  "Enable compile-time AVX support."  ${ENABLE_AVX_DEFAULT})
option(ENABLE_AVX2 "Enable compile-time AVX2 support." ${ENABLE_AVX2_DEFAULT})

if(CMAKE_COMPILER_IS_GNUCXX)
  if (ENABLE_SSE2)
    if (!HAVE_EMMINTRIN_H)
      message( SEND_ERROR "SSE2 is enabled, but is not supported by compiler.")
    else()
      add_compile_flags("C;CXX" "-msse2")
      message(STATUS "SSE2 is enabled - target CPU must support it")
    endif()
  endif()
  
  if (ENABLE_SSE3)
    if (!HAVE_PMMINTRIN_H)
      message( SEND_ERROR "SSE3 is enabled, but is not supported by compiler.")
    else()
      add_compile_flags("C;CXX" "-msse3")
      message(STATUS "SSE3 is enabled - target CPU must support it")
    endif()
  endif()
  
  if (ENABLE_SSE4)
    if (!HAVE_SMMINTRIN_H OR !HAVE_NMMINTRIN_H)
      message( SEND_ERROR "SSE4 is enabled, but is not supported by compiler.")
    else()
      add_compile_flags("C;CXX" "-msse4")
      message(STATUS "SSE4 is enabled - target CPU must support it")
    endif()
  endif()
  
  if (ENABLE_AVX)
    if (!HAVE_IMMINTRIN_H)
      message(SEND_ERROR "AVX is enabled")
    else()
      add_compile_flags("C;CXX" "-mavx")
      message(STATUS "AVX is enabled - target CPU must support it")
    endif()
  endif()
  
  if (ENABLE_AVX2)
    if (!HAVE_ZMMINTRIN_H)
      message(SEND_ERROR "AVX2 is enabled")
    else()
      add_compile_flags("C;CXX" "-mavx2")
      message(STATUS "AVX2 is enabled - target CPU must support it")
    endif()
  endif()
  
  if (ENABLE_FMA)
    if (!HAVE_FMAINTRIN_H)
      message(SEND_ERROR "FMA is enabled")
    else()
      add_compile_flags("C;CXX" "-mfma")
      message(STATUS "FMA is enabled - target CPU must support it")
    endif()
  endif()
endif()

if (MSVC)
  if (ENABLE_AVX2)
    if (!HAVE_ZMMINTRIN_H)
      message(SEND_ERROR "AVX2 is enabled")
    else()
      add_compile_flags("C;CXX" "/arch:AVX2")
      message(STATUS "AVX2 is enabled - target CPU must support it")
    endif()
  endif()
endif()

# TODO:
#   #if defined(_MSC_VER)
#        /* Microsoft C/C++-compatible compiler */
#        #include <intrin.h>
#   #elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#        /* GCC-compatible compiler, targeting x86/x86-64 */
#        #include <x86intrin.h>
#   #elif defined(__GNUC__) && defined(__ARM_NEON__)
#        /* GCC-compatible compiler, targeting ARM with NEON */
#        #include <arm_neon.h>
#   #elif defined(__GNUC__) && defined(__IWMMXT__)
#        /* GCC-compatible compiler, targeting ARM with WMMX */
#        #include <mmintrin.h>
#   #elif (defined(__GNUC__) || defined(__xlC__)) && (defined(__VEC__) || defined(__ALTIVEC__))
#        /* XLC or GCC-compatible compiler, targeting PowerPC with VMX/VSX */
#        #include <altivec.h>
#   #elif defined(__GNUC__) && defined(__SPE__)
#        /* GCC-compatible compiler, targeting PowerPC with SPE */
#        #include <spe.h>
#   #endif
# UltraSPARC+VIS with gcc, use visintrin.h or if you have VSDK, use vis.h

