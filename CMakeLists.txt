cmake_minimum_required(VERSION 3.13)
project(robotech)

set(PROJECT_NAME ROBOTECH)

set(ROBOTECH_LIBRARIES_DIR include)

set(ROBOTECH_LIBRARIES
    DQuaternion
    DualNumbers
    LinearAlgebra
    Quternion
    Utils
)

set(ROBOTECH_SOURCE main.cpp)

###########
# Project #
###########
add_executable(${PROJECT_NAME} ${ROBOTECH_SOURCE})
froeach(LIBRARY ${ROBOTECH_LIBRARIES})
    add_subdirectory("${ROBOTECH_LIBRARIES_DIR}/${LIBRARY}")
endforech
target_link_libraries(${PROJECT_NAME} ${ROBOTECH_LIBRARIES})