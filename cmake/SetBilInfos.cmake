


# Bil version
# -----------
file(STRINGS ${BIL_PATH}/VERSION BIL_VERSION)
#message("BIL_VERSION = ${BIL_VERSION}")


# General informations
# --------------------
set(BIL_SHORT_LICENSE "GNU General Public License")
string(TIMESTAMP BIL_DATE)
string(TIMESTAMP BIL_YEAR %Y)
cmake_host_system_information(RESULT BIL_HOST QUERY HOSTNAME)
#set(BIL_HOST    ${shell hostname}: ${shell hostname -I})
#cmake_host_system_information(RESULT BIL_DISTRIB QUERY DISTRIB_INFO)
set(BIL_OS ${CMAKE_HOST_SYSTEM_NAME})
set(BIL_URL       "http://bil.ifsttar.fr")
set(BIL_EMAIL     "patrick.dangla@univ-eiffel.fr")
set(BIL_COPYRIGHT "Copyright \(C\) 2002 Patrick Dangla")
set(BIL_PROGNAME  "Bil, a modeling platform based on FEM/FVM")
execute_process(COMMAND date "+%Y%m%d" OUTPUT_VARIABLE BIL_DATE 
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND hostname OUTPUT_VARIABLE BIL_HOSTNAME 
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND whoami OUTPUT_VARIABLE BIL_PACKAGER 
                OUTPUT_STRIP_TRAILING_WHITESPACE)
                
               

configure_file(${BIL_PATH}/BilInfo.h.in 
               ${BIL_PATH}/BilInfo.h)
configure_file(${BIL_PATH}/BilVersion.h.in
               ${BIL_PATH}/BilVersion.h)
configure_file(${BIL_PATH}/BilPath.h.in
               ${BIL_PATH}/BilPath.h)
               
               


#[[
if(APPLE)
  set(BIL_OS "MacOSX")
elseif(CYGWIN)
  set(BIL_OS "Windows")
else()
  set(BIL_OS "${CMAKE_SYSTEM_NAME}")
endif()
]]
