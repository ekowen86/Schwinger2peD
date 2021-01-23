add_executable(wilson2peD wilson2peD.cpp)
target_link_libraries(wilson2peD schwinger2peD)
target_include_directories(wilson2peD PUBLIC ${PROJECT_SOURCE_DIR}/include)

