# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/CMake/3.14.3/bin/cmake

# The command to remove a file.
RM = /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/CMake/3.14.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build

# Include any dependencies generated for this target.
include src/TP4/CMakeFiles/img.exe.dir/depend.make

# Include the progress variables for this target.
include src/TP4/CMakeFiles/img.exe.dir/progress.make

# Include the compile flags for this target's objects.
include src/TP4/CMakeFiles/img.exe.dir/flags.make

src/TP4/CMakeFiles/img.exe.dir/img.cpp.o: src/TP4/CMakeFiles/img.exe.dir/flags.make
src/TP4/CMakeFiles/img.exe.dir/img.cpp.o: ../src/TP4/img.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/TP4/CMakeFiles/img.exe.dir/img.cpp.o"
	cd /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/src/TP4 && /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Core/GCCcore/7.3.0/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/img.exe.dir/img.cpp.o -c /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/src/TP4/img.cpp

src/TP4/CMakeFiles/img.exe.dir/img.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/img.exe.dir/img.cpp.i"
	cd /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/src/TP4 && /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Core/GCCcore/7.3.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/src/TP4/img.cpp > CMakeFiles/img.exe.dir/img.cpp.i

src/TP4/CMakeFiles/img.exe.dir/img.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/img.exe.dir/img.cpp.s"
	cd /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/src/TP4 && /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Core/GCCcore/7.3.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/src/TP4/img.cpp -o CMakeFiles/img.exe.dir/img.cpp.s

# Object files for target img.exe
img_exe_OBJECTS = \
"CMakeFiles/img.exe.dir/img.cpp.o"

# External object files for target img.exe
img_exe_EXTERNAL_OBJECTS =

src/TP4/img.exe: src/TP4/CMakeFiles/img.exe.dir/img.cpp.o
src/TP4/img.exe: src/TP4/CMakeFiles/img.exe.dir/build.make
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_mpi.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_thread.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_serialization.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_chrono.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_system.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_program_options.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_regex.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_filesystem.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/tbb/2018_U6/lib/libtbb.so
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCC/7.3.0-2.30/impi/2018.3.222/lib64/libmpi.so
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_mpi.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_thread.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_serialization.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_chrono.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_system.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_program_options.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_regex.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/Boost/1.67.0/lib/libboost_filesystem.a
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/tbb/2018_U6/lib/libtbb.so
src/TP4/img.exe: ../contrib/opencv/3.4.7/lib64/libopencv_highgui.so.3.4.7
src/TP4/img.exe: /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCC/7.3.0-2.30/impi/2018.3.222/lib64/libmpi.so
src/TP4/img.exe: ../contrib/opencv/3.4.7/lib64/libopencv_videoio.so.3.4.7
src/TP4/img.exe: ../contrib/opencv/3.4.7/lib64/libopencv_imgcodecs.so.3.4.7
src/TP4/img.exe: ../contrib/opencv/3.4.7/lib64/libopencv_imgproc.so.3.4.7
src/TP4/img.exe: ../contrib/opencv/3.4.7/lib64/libopencv_core.so.3.4.7
src/TP4/img.exe: src/TP4/CMakeFiles/img.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable img.exe"
	cd /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/src/TP4 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/img.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/TP4/CMakeFiles/img.exe.dir/build: src/TP4/img.exe

.PHONY : src/TP4/CMakeFiles/img.exe.dir/build

src/TP4/CMakeFiles/img.exe.dir/clean:
	cd /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/src/TP4 && $(CMAKE_COMMAND) -P CMakeFiles/img.exe.dir/cmake_clean.cmake
.PHONY : src/TP4/CMakeFiles/img.exe.dir/clean

src/TP4/CMakeFiles/img.exe.dir/depend:
	cd /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/src/TP4 /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/src/TP4 /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/src/TP4/CMakeFiles/img.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/TP4/CMakeFiles/img.exe.dir/depend

