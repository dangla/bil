
# Function to get the list of files of extension "ext" in a given directory
function(list_files dir ext return_var)
  #message("Entry in list_files!")
  #message("dir: ${dir}")
  #message("return_var: ${return_var}")
  unset(list)
  file(GLOB list LIST_DIRECTORIES false RELATIVE ${dir} ${dir}/*.${ext})
  #message("list: ${list}")
  set(${return_var} ${list} PARENT_SCOPE)
endfunction()
