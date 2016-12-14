cmake_minimum_required(VERSION 3.2)
project(SPHSm6)



set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Find the QtWidgets library
set (CMAKE_PREFIX_PATH "/home/nikita/Qt5.7.0/5.7/gcc_64")

find_package(Qt5Widgets)

set(
        SOURCE_FILES
        main.cpp
        Particile.cpp Particile.h
        Flow.cpp Flow.h
        Help.cpp Help.h
        SpaceParsing.cpp SpaceParsing.h
        Calculator.cpp Calculator.h
        MainWindow.cpp MainWindow.h
        MyQGraphicsView.cpp myqgraphicsview.h
        Flow_drawer.cpp Flow_drawer.h
        Visualisator.cpp Visualisator.h W_functions.cpp W_functions.h)


# Tell CMake to create the helloworld executable
add_executable(SPHSm6  ${SOURCE_FILES})


# Use the Widgets module from Qt 5.
target_link_libraries(SPHSm6 Qt5::Widgets)
