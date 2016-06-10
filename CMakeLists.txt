<<<<<<< HEAD
cmake_minimum_required(VERSION 3.2)
project(sph_interface)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Find the QtWidgets library
find_package(Qt5Widgets)

=======
cmake_minimum_required(VERSION 3.3)
project(SPHSm6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

>>>>>>> 883cc5708c435c9115b864d6fdf26606dffc1703
set(
        SOURCE_FILES
        main.cpp
        Particile.cpp Particile.h
        Flow.cpp Flow.h
        Help.cpp Help.h
        SpaceParsing.cpp SpaceParsing.h
        Calculator.cpp Calculator.h
<<<<<<< HEAD
        MainWindow.cpp MainWindow.h
        MyQGraphicsView.cpp myqgraphicsview.h
)


# Tell CMake to create the helloworld executable
add_executable(sph_interface  ${SOURCE_FILES})


# Use the Widgets module from Qt 5.
target_link_libraries(sph_interface Qt5::Widgets)
=======
        )
add_executable(SPHSm6 ${SOURCE_FILES})
>>>>>>> 883cc5708c435c9115b864d6fdf26606dffc1703
