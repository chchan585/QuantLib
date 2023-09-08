set (myApplicationFiles
    myApplication.hpp
    myApplication.cpp
)

add_executable(myApplication ${myApplicationFiles})

target_link_libraries(myApplication ql_library ${QL_THREAD_LIBRARIES})
if (QL_INSTALL_EXAMPLES)
    install(TARGETS CDS RUNTIME DESTINATION ${QL_INSTALL_EXAMPLESDIR})
endif()
