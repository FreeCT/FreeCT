set(TARGET_NAME fct_icd)

file( GLOB SOURCES src/*.cpp )
add_executable(${TARGET_NAME} ${SOURCES})
add_dependencies(${TARGET_NAME} fct_read)

target_include_directories(${TARGET_NAME} PRIVATE ${INCLUDE_DIRS} src/)
target_link_libraries(${TARGET_NAME} PRIVATE yaml-cpp fct_read m)
