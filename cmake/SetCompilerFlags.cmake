
# Flags
# -----
# Warning flags
set(WFLAGS
    "-Wall -Wcast-align -Wcast-qual -Wnested-externs -Wpointer-arith -Wmissing-prototypes -Wstrict-prototypes -Wshadow -fno-common -Dinline= -Wvariadic-macros" 
    #-W -Wtraditional -Wconversion -Wwrite-strings -Werror -fshort-enums -Wunused-macros
)
    

# Optimization flags
set(OFLAGS
    "-gdwarf-2 -g3"#-O2
)

set(CMAKE_C_FLAGS  "-pedantic -fno-common ${WFLAGS}")
set(CMAKE_C_FLAGS_DEBUG          "-g3 -DDEBUG")
set(CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -g3 -gdwarf-3")
set(CMAKE_Fortran_FLAGS "-cpp")


#[[
if( NOT C_FLAGS_INITIALIZED )
  # only do this on the first pass through to avoid overwriting user added options.
  set( C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Are compiler flags already set?" )

  # Overwrite CMake's defaults...
  if( "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set( CMAKE_C_FLAGS                "-pedantic -fno-common ${WFLAGS}" )
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
