project(BigInteger)
cmake_minimum_required(VERSION 2.8)

set(SOURCE
    gcd.cpp
    BigInteger.cpp
    BigIntegerSign.cpp
)

set(SOURCE_EXAMPLE
    test.cpp
    main.cpp
)

set(WARNING_FLAGS "-pedantic -Wall -Wextra -Wshadow -Wformat=2 -Wfloat-equal -Wlogical-op -Wcast-qual -Wcast-align")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_FLAGS} -std=c++14 -g")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D _GLIBCXX_DEBUG -fsanitize=address")

add_library(big_number_static STATIC ${SOURCE})
add_library(big_number_shared SHARED ${SOURCE})
set_target_properties(big_integer_static PROPERTIES OUTPUT big_integer)
set_target_properties(big_integer_shared PROPERTIES OUTPUT big_integer)

add_executable(big_int ${SOURCE_EXAMPLE})
target_link_libraries(big_int big_number_static)