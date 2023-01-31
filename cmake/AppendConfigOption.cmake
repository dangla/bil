# Append config_option with VARNAME
macro(append_config_option VARNAME STRING config_option)
  set(${VARNAME} TRUE)
  list(APPEND ${config_option} ${STRING})
  list(SORT ${config_option})
  message(STATUS "Found " ${STRING})
endmacro(append_config_option)
