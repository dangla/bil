
file(STRINGS OPTIONS opts REGEX "^[^# ]")

foreach(opt IN ITEMS ${opts})
  #message("opt = ${opt}")
  string(REGEX MATCH "^[^ ]*" opt_name "${opt}")
  string(REGEX REPLACE "^[ ]*[^ ]+[ ]+=[ ]+|[ ]+$" "" opt_bool "${opt}")

  #message("opt_name = ${opt_name}")
  #message("opt_bool = ${opt_bool}")
  
  if("${opt_bool}")
    
    string(REGEX REPLACE "ENABLE" "HAVE" opt_have "${opt_name}")
    
    #message("opt_have = ${opt_have}")
    
    set(${opt_have} TRUE)

  endif()
endforeach()


configure_file(${BIL_PATH}/BilConfig.h.in 
               ${BIL_PATH}/BilConfig.h)

