file( GLOB SOURCES src/*.cpp )
add_library(fct_read ${SOURCES})

target_include_directories(fct_read PUBLIC include/ ${INCLUDE_DIRS})
#target_link_libraries(fct_read PUBLIC gdcmCommon)

file(GLOB HEADERS include/*)
install(FILES ${HEADERS} DESTINATION include/fct/)
install(TARGETS fct_read DESTINATION lib/)