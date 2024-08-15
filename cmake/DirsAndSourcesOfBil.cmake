macro(dirs_and_sources_of_bil dirsall srcall srcmain)
#[[
  The list of directories to consider in the executable:
  i.e. the children, grandchildren, grandgrandchildren of
  ${CMAKE_CURRENT_SOURCE_DIR} = src
#]]

set(BIL_DIRMODELS  ${CMAKE_CURRENT_SOURCE_DIR}/Models/ModelFiles)
set(BIL_DIRMODULES ${CMAKE_CURRENT_SOURCE_DIR}/Modules/ModuleFiles)
set(BIL_DIRPLASTICITYMODELS ${CMAKE_CURRENT_SOURCE_DIR}/Models/ConstitutiveLaws/PlasticityModels)
set(BIL_DIRDAMAGEMODELS ${CMAKE_CURRENT_SOURCE_DIR}/Models/ConstitutiveLaws/DamageModels)

include(AppendBilDirectories)

# Set the list of directories within src (one level depth)
append_bil_directories("${CMAKE_CURRENT_SOURCE_DIR}" ${dirsall})
# Append with the list of directories within the previous ones (two level depth)
append_bil_directories("${${dirsall}}" ${dirsall})
# Append with the list of directories within the previous ones (three level depth)
append_bil_directories("${${dirsall}}" ${dirsall})

set(BIL_DIRSOTHER "${${dirsall}}")
list(REMOVE_ITEM BIL_DIRSOTHER ${BIL_DIRMODELS} ${BIL_DIRMODULES})

#[[
message(STATUS "directories under ${CMAKE_CURRENT_SOURCE_DIR} within three level depth:")
message("directory of the model files:")
message("${BIL_DIRMODELS}")
message("directory of the module files:")
message("${BIL_DIRMODULES}")
message("other directories:")
message("${BIL_DIRSOTHER}")
message("all directories:")
message("${${dirsall}}")
#]]


#[[
  The list of source files to be considered in the library and the executable
#]]
 
# The main entry
#set(BIL_SRCMAIN ${CMAKE_CURRENT_SOURCE_DIR}/Main/Main.[c,cpp])
file(GLOB BIL_SRCMAIN LIST_DIRECTORIES false ${CMAKE_CURRENT_SOURCE_DIR}/Main/Main.c*)
set(${srcmain} ${BIL_SRCMAIN})
#[[
message("Main sources under ${CMAKE_CURRENT_SOURCE_DIR}/Main:")
message("BIL_SRCMAIN = ${BIL_SRCMAIN}")
#]]


include(AppendSelectedModels)
include(AppendBilFiles)

# The model files
append_bil_files("${BIL_DIRMODELS}" "c" BIL_SRCMODELS)
append_bil_files("${BIL_DIRMODELS}" "cpp" BIL_SRCMODELS)
append_selected_models(${BIL_DIRMODELS}/../ListOfModels.inc
                       BIL_SRCMODELS)
#list(TRANSFORM BIL_SRCMODELS PREPEND "${BIL_DIRMODELS}/")

# The module files
append_bil_files("${BIL_DIRMODULES}" "c" BIL_SRCMODULES)
append_bil_files("${BIL_DIRMODULES}" "cpp" BIL_SRCMODULES)
append_selected_models(${BIL_DIRMODULES}/../ListOfModules.inc
                       BIL_SRCMODULES)
#list(TRANSFORM BIL_SRCMODULES PREPEND "${BIL_DIRMODULES}/")

# The plasticity model files
append_bil_files("${BIL_DIRPLASTICITYMODELS}" "c" BIL_SRCPLASTICITYMODELS)
append_bil_files("${BIL_DIRPLASTICITYMODELS}" "cpp" BIL_SRCPLASTICITYMODELS)
# We remove the file so that all the plasticity models are included
file(REMOVE ${BIL_DIRPLASTICITYMODELS}/ListOfPlasticityModels.inc)
append_selected_models(${BIL_DIRPLASTICITYMODELS}/ListOfPlasticityModels.inc
                       BIL_SRCPLASTICITYMODELS)

# The damage model files
append_bil_files("${BIL_DIRDAMAGEMODELS}" "c" BIL_SRCDAMAGEMODELS)
append_bil_files("${BIL_DIRDAMAGEMODELS}" "cpp" BIL_SRCDAMAGEMODELS)
# We remove the file so that all the damage models are included
file(REMOVE ${BIL_DIRDAMAGEMODELS}/ListOfDamageModels.inc)
append_selected_models(${BIL_DIRDAMAGEMODELS}/ListOfDamageModels.inc
                       BIL_SRCDAMAGEMODELS)


# The list of source files (f) in the directories built above.
append_bil_files("${BIL_DIRSOTHER}" "f" BIL_SRCFORTRAN)
# The list of source files (c,cpp) in the directories built above.
append_bil_files("${BIL_DIRSOTHER}" "c" BIL_SRCOTHER)
append_bil_files("${BIL_DIRSOTHER}" "cpp" BIL_SRCOTHER)
# Remove the main entry
list(REMOVE_ITEM BIL_SRCOTHER ${BIL_SRCMAIN})

# All source files
set(${srcall} ${BIL_SRCOTHER} ${BIL_SRCMODELS} ${BIL_SRCMODULES} ${BIL_SRCFORTRAN})

#[[
message(STATUS "source files under ${CMAKE_CURRENT_SOURCE_DIR} within three level depth:")
message("sources of models under ${BIL_DIRMODELS}:")
message("BIL_SRCMODELS = ${BIL_SRCMODELS}")
message("sources of modules under ${BIL_DIRMODULES}:")
message("BIL_SRCMODULES = ${BIL_SRCMODULES}")
message("sources of plasticity models:")
message("BIL_SRCMODULES = ${BIL_SRCPLASTICITYMODELS}")
message("all source files:")
message("${srcall} = ${${srcall}}")
#]]

unset(BIL_DIRMODELS)
unset(BIL_DIRMODULES)
unset(BIL_DIRPLASTICITYMODELS)
unset(BIL_DIRDAMAGEMODELS)


# Compile definition: BASENAME, __CPLUSPLUS
foreach(i IN ITEMS ${BIL_SRCOTHER} ${BIL_SRCMODELS} ${BIL_SRCMODULES})
  get_filename_component(b ${i} NAME_WE)
  set_source_files_properties(${i} PROPERTIES COMPILE_DEFINITIONS "BASENAME=${b}")
endforeach()


# Force the use of c++ or not
# [[
foreach(i IN ITEMS ${BIL_SRCMAIN} ${BIL_SRCOTHER} ${BIL_SRCMODELS} ${BIL_SRCMODULES})
  # We don't
  #set_properties(SOURCE ${i} APPEND PROPERTIES COMPILE_DEFINITIONS "__CPLUSPLUS=__cplusplus")
  # We do
  set_source_files_properties(${i} PROPERTIES LANGUAGE CXX)
endforeach()
# ]]

endmacro()
