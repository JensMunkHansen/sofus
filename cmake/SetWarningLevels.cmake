# Threat warnings as errors 
option(WARNING_TREAT_AS_ERRORS "Warnings as errors" OFF)

if (MSVC)
  option(VERIFICATION_CODE_ANALYSIS "Code analysis" OFF)
  if (VERIFICATION_CODE_ANALYSIS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -analyze")
  endif()

  # Set warning level to 4
  set(WARNING_LEVEL "-W4" CACHE STRING "Warnings level")
  set_property(CACHE WARNING_LEVEL PROPERTY STRINGS "-W4" "-W3" "-W2" "-W1" "-W0")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_LEVEL}")
  if (WARNING_TREAT_AS_ERRORS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX")
  endif()
else()
  set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -Wno-overflow")
  if (WARNING_TREAT_AS_ERRORS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  endif()
endif(MSVC)

if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4661")
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -Wno-overflow")
endif()

