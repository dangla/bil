macro(dirs_and_sources_of_bil dirsall srcall)
#[[
  The list of directories to consider in the executable:
  i.e. the children and grandchildren of src
]]

set(BIL_DIRMODELS  ${CMAKE_CURRENT_SOURCE_DIR}/Models/ModelFiles)
set(BIL_DIRMODULES ${CMAKE_CURRENT_SOURCE_DIR}/Modules/ModuleFiles)

include(AppendBilDirectories)
append_bil_directories("${CMAKE_CURRENT_SOURCE_DIR}" ${dirsall})
append_bil_directories("${${dirsall}}" ${dirsall})

set(BIL_DIRSOTHER "${${dirsall}}")
list(REMOVE_ITEM BIL_DIRSOTHER ${BIL_DIRMODELS} ${BIL_DIRMODULES})

#message(STATUS "directories under src within two level depth:")
#message("directory of the model files:")
#message("${BIL_DIRMODELS}")
#message("directory of the module files:")
#message("${BIL_DIRMODULES}")
#message("other directories:")
#message("${BIL_DIRSOTHER}")


#[[
  The list of source files to be considered in the library and the executable
 ]]
 
# The main entry
set(BIL_SRCMAIN ${CMAKE_CURRENT_SOURCE_DIR}/Main/Main.c)

# The model files
include(AppendSelectedModels)
append_selected_models(${CMAKE_CURRENT_SOURCE_DIR}/Models/ListOfModels.inc BIL_SRCMODELS)
list(TRANSFORM BIL_SRCMODELS PREPEND "${BIL_DIRMODELS}/")

# The module files
append_selected_models(${CMAKE_CURRENT_SOURCE_DIR}/Modules/ListOfModules.inc BIL_SRCMODULES)
list(TRANSFORM BIL_SRCMODULES PREPEND "${BIL_DIRMODULES}/")

# The list of source files (f) in the directories built above.
include(AppendBilFiles)
append_bil_files("${BIL_DIRSOTHER}" "f" BIL_SRCFORTRAN)
# The list of source files (c) in the directories built above.
append_bil_files("${BIL_DIRSOTHER}" "c" BIL_SRCOTHER)
# Remove the main entry
list(REMOVE_ITEM BIL_SRCOTHER ${BIL_SRCMAIN})

# All source files
set(${srcall} ${BIL_SRCOTHER} ${BIL_SRCMODELS} ${BIL_SRCMODULES} ${BIL_SRCFORTRAN})


#message("BIL_SRCMODELS = ${BIL_SRCMODELS}")
#message("BIL_SRCMODULES = ${BIL_SRCMODULES}")
#message(STATUS "source c and f files under src within two level depth:")
#message("${${srcall}}")
#message(STATUS "header files under src within two level depth:")
#message("${BIL_HDRALL}")


# Compile definitions
foreach(i IN ITEMS ${BIL_SRCMODELS} ${BIL_SRCMODULES})
  get_filename_component(b ${i} NAME_WE)
  set_source_files_properties(${i} PROPERTIES COMPILE_DEFINITIONS "BASENAME=${b}")
endforeach()
endmacro()
