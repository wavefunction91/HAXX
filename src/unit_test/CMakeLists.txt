#set(Boost_USE_STATIC_LIBS ON)
find_package (Boost COMPONENTS unit_test_framework REQUIRED)

add_executable(haxx_unary_op haxx_unary_op.cxx)
add_executable(haxx_binary_op haxx_binary_op.cxx)
add_executable(haxx_algebra haxx_algebra.cxx)

target_link_librarieS(haxx_unary_op ${Boost_LIBRARIES})
target_link_librarieS(haxx_binary_op ${Boost_LIBRARIES})
target_link_librarieS(haxx_algebra ${Boost_LIBRARIES})


add_test(HAXX_UNARY_OP haxx_unary_op --report_level=detailed )  
add_test(HAXX_BINARY_OP haxx_binary_op --report_level=detailed )  
add_test(HAXX_ALGEBRA haxx_algebra --report_level=detailed )  
