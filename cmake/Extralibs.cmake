
# Re-build the library file
message("Updating the library file: ${BILEXTRALIBS_H}")
file(REMOVE ${BILEXTRALIBS_H})
file(TOUCH  ${BILEXTRALIBS_H})

file(APPEND ${BILEXTRALIBS_H} "#ifndef BILEXTRALIBS_H\n")
file(APPEND ${BILEXTRALIBS_H} "#define BILEXTRALIBS_H\n")

file(STRINGS EXTRALIBS libs REGEX "^[^# ]")

foreach(lib ${libs})
  #message("lib = ${lib}")
  string(REGEX MATCH "^[A-Z]+" lib_name "${lib}")
  string(REGEX REPLACE "^[ ]*[A-Z]+[ ]+=[ ]+|[ ]+$" "" lib_path "${lib}")

  #message("lib_name = ${lib_name}")
  #message("lib_path = ${lib_path}")
  
  if(EXISTS "${lib_path}")
  
    get_filename_component(lib_full_path ${lib_path} ABSOLUTE)
    
    if(NOT DEFINED ${lib_name}_DIR)
      get_filename_component(${lib_name}_DIR ${lib_full_path} DIRECTORY CACHE)
      message("The ${lib_name} library found at " ${lib_full_path})
      message("is added to the cache")
    endif()
    
    file(APPEND ${BILEXTRALIBS_H} "#define ${lib_name}LIB  ${${lib_name}_DIR}\n")
    set(BIL_EXTRALIBS ${BIL_EXTRALIBS} ${lib_full_path})
    
  endif()
endforeach()


file(APPEND ${BILEXTRALIBS_H} "#endif")


if(NOT "${BIL_EXTRALIBS}" STREQUAL "")
  message("Some extra-libraries are used.")
  message("Full paths to these extra-libraries:")
  message("${BIL_EXTRALIBS}")
else()
  message("No extra-libraries are used")
endif()


#welcome:
#       read -p "Do you use BLAS library (Y/N)?:" BLAS_USE
