
# Warning flags
# -------------
# Warning flags common to C/C++
set(WFLAGSmini
    "-Wpedantic -Wall -Wextra"
)

set(WFLAGS
    "${WFLAGSmini} -Wformat=2 -Wconversion -Wsign-conversion -Wtrampolines -Wimplicit-fallthrough -Wcast-align -Wcast-qual -Wpointer-arith -Wshadow"
)

set(WFLAGS1
    "${WFLAGS} -Wctor-dtor-privacy -Wdisabled-optimization -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Woverloaded-virtual -Wredundant-decls -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef"
)

# Warning flags for C
set(WCFLAGS
    "${WFLAGS} -Wnested-externs -Wmissing-prototypes -Wstrict-prototypes"
)

# Warning flags for C++
set(WCXXFLAGS
    "${WFLAGS} -Wnon-virtual-dtor -Wzero-as-null-pointer-constant -Wunused -Woverloaded-virtual -Wfloat-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference"
)
    

# Optimization flags
set(OFLAGS
    "-gdwarf-2 -g3"#-O2
)

# C compiler flags
set(CFLAGS
    "${WCFLAGS}"
)

# C++ compiler flags
# fpermissive is not satisfactory but at the present time it is needed
# because most of the implementations adopt a C-style which are not 
# allowed in C++.
set(CXXFLAGS
    "${WCXXFLAGS} -fpermissive"
)

set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${CFLAGS}")
set(CMAKE_C_FLAGS_DEBUG          "-g3 -DDEBUG")
set(CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -g3 -gdwarf-3")


set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${CXXFLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "-g3 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g3 -gdwarf-3")


set(CMAKE_Fortran_FLAGS "-cpp")


#[[
if( NOT C_FLAGS_INITIALIZED )
  # only do this on the first pass through to avoid overwriting user added options.
  set( C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Are compiler flags already set?" )

  # Overwrite CMake's defaults...
  if( "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set( CMAKE_C_FLAGS                "-pedantic -fno-common ${WCFLAGS}" )
    set( CMAKE_C_FLAGS_DEBUG          "-g -DDEBUG")
    set( CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG" )
    set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
    set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -g -gdwarf-3" )
    # in a similar fashion, provide CXX_FLAGS and Fortran_FLAGS

  elseif( "${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC" )
    set( CMAKE_C_FLAGS "/fp:precise /DWIN32 /D_WINDOWS /MP" )
    ...
  endif()
endif()

# Save the current compiler flags to the cache every time cmake configures the project.
set(CMAKE_C_FLAGS                "${CMAKE_C_FLAGS}"                CACHE
     STRING "compiler flags" FORCE)
set(CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}"          CACHE
     STRING "compiler flags" FORCE)
set(CMAKE_C_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}"        CACHE
     STRING "compiler flags" FORCE)
set(CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_MINSIZEREL}"     CACHE
     STRING "compiler flags" FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" CACHE
     STRING "compiler flags" FORCE)
# and similar for CXX_FLAGS and Fortran_FLAGS...
]]
