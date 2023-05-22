macro(create_bil_header headers bilheader)

  #message("headers: ${headers}")
  #message("bilheader: ${bilheader}")

  file(REMOVE ${bilheader})
  file(TOUCH  ${bilheader})
  
  #get_filename_component(bilheader_dir ${bilheader} DIRECTORY)
  #file(COPY ${headers} DESTINATION ${bilheader_dir})
  
  foreach(header IN ITEMS ${headers})
    get_filename_component(base ${header} NAME)

    #message("base: ${base}")
    
    file(APPEND ${bilheader} "#include \"${base}\"\n")
  endforeach()
endmacro()
