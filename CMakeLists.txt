cmake_minimum_required(VERSION 3.0)

project(a4)
set(CMAKE_CXX_STANDARD 11)

find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)
find_package(SDL2 REQUIRED)

add_library(col781 animation/numerical_methods.cpp animation/animation.cpp animation/halfedge.cpp animation/shapes.cpp src/hw.cpp src/camera.cpp deps/src/gl.c)
target_include_directories(col781 PUBLIC deps/include)
target_link_libraries(col781 glm::glm OpenGL::GL SDL2::SDL2)

add_executable(example executables/example.cpp)
target_include_directories(example PUBLIC src)
target_link_libraries(example col781)

add_executable(example2 executables/example2.cpp)
target_include_directories(example2 PUBLIC src)
target_link_libraries(example2 col781)

add_executable(spring_example executables/spring_example.cpp)
target_include_directories(spring_example PUBLIC src)
target_link_libraries(spring_example col781)

add_executable(collision_example executables/collision_example.cpp)
target_include_directories(collision_example PUBLIC src)
target_link_libraries(collision_example col781)

add_executable(collision_example2 executables/collision_example2.cpp)
target_include_directories(collision_example2 PUBLIC src)
target_link_libraries(collision_example2 col781)
