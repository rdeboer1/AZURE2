include_directories(../include ${AZURE2_SOURCE_DIR}/include)

set (AZUREAPI_SOURCES main.cpp AZUREAPI.cpp AZURESocket.cpp)

set (AZUREAPI_HEADERS ../include/AZUREAPI.h
../include/AZURESocket.h)

add_library(AZUREAPI ${AZUREAPI_SOURCES} ${AZUREAPI_HEADERS_MOC})

# SuperSubatomic - added this
#target_link_libraries (AZURESetup Qt5::Core Qt5::PrintSupport)
#add_executable(AZURE2 AZURESetup.qrc)

set (TEMP_VAR AZUREAPI )
set (AZUREAPI_LIBRARIES ${TEMP_VAR} PARENT_SCOPE)