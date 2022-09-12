# Append selected models found in a file
macro(append_selected_models file selectedmodelfiles)
  #message("Entry in append_selected_models")
  #[[
  On input ${selectedmodelfiles} contains a list of model files from which
  we should select some of them. Two cases occur:
    1. if ${file} exists, we read the selected models from it. 
    2. if ${file} doesn't exist, we consider all the models found in
       ${selectedmodelfiles}. At the same time we create the file ${file}
       with all the models in it.
  On output ${selectedmodelfiles} contains the selected models.
  #]]
  if(NOT EXISTS ${file})
    set(selectednames)
    foreach(i IN ITEMS ${${selectedmodelfiles}})
      get_filename_component(b ${i} NAME_WE)
      list(APPEND selectednames "${b}")
    endforeach()
    string(REPLACE ";" " " selectednames "${selectednames}")
    file(WRITE ${file} "SELECTEDMODELS = ${selectednames}")
    unset(selectednames)
  else()
    set(selectednames)
    file(STRINGS ${file} selectednames)
    string(REGEX REPLACE "^[ ]*[A-Z]+[ ]+=[ ]+|[ ]+$" "" selectednames "${selectednames}")
    string(REGEX REPLACE "[ ]+" ";" selectednames "${selectednames}")
    #message("SELECTEDNAMES = ${selectednames}")
    set(selection)
    string(REGEX REPLACE ";" "|" selection "${selectednames}")
    #message("selection = ${selection}")
    #message("Before FILTER: ${${selectedmodelfiles}}")
    list(FILTER ${selectedmodelfiles} INCLUDE REGEX ".*(${selection})\\.(c|cpp)")
    #message("After FILTER: ${${selectedmodelfiles}}")
    unset(selection)
    unset(selectednames)
  endif()

  #[[
  Configure the file "ListOfModels.h" or "ListOfModules.h"
  #]]
  set(NBOFSELECTEDMODELS)
  set(SELECTEDMODELMETHODS)
  set(SELECTEDMODELNAMES)
  list(LENGTH ${selectedmodelfiles} NBOFSELECTEDMODELS)
  foreach(i IN ITEMS ${${selectedmodelfiles}})
    get_filename_component(b ${i} NAME_WE)
    list(APPEND SELECTEDMODELMETHODS "${b}##m")
    list(APPEND SELECTEDMODELNAMES   "\"${b}\"")
  endforeach()
  string(REPLACE ";" "," SELECTEDMODELMETHODS "${SELECTEDMODELMETHODS}")
  string(REPLACE ";" "," SELECTEDMODELNAMES "${SELECTEDMODELNAMES}")
  #[[
  message("SELECTEDMODELMETHODS = ${SELECTEDMODELMETHODS}")
  message("SELECTEDMODELNAMES = ${SELECTEDMODELNAMES}")
  #]]
  string(REPLACE ".inc" ".h.in" file.h.in ${file})
  string(REPLACE ".inc" ".h" file.h ${file})
  #[[
  message("file.h.in = ${file.h.in}")
  message("file.h = ${file.h}")
  #]]
  configure_file(${file.h.in}  ${file.h})
  unset(NBOFSELECTEDMODELS)
  unset(SELECTEDMODELMETHODS)
  unset(SELECTEDMODELNAMES)
endmacro()
