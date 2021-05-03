cmake_minimum_required(VERSION 3.19)
project(Visual_Studio_CPPSerialLeArq)

set(CMAKE_CXX_STANDARD 14)

add_executable(Visual_Studio_CPPSerialLeArq
        initialize_serial.cpp
        read_u24h_serial_data.m
        u24h_serial.cpp)
