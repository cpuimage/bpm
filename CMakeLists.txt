cmake_minimum_required(VERSION 2.8)
project(bpm LANGUAGES C)
add_executable(bpm bpm.c)
target_link_libraries(bpm -lm)