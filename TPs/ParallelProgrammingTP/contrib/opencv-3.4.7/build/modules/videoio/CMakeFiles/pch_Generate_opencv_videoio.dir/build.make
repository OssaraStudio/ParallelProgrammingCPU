# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

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
CMAKE_COMMAND = /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/CMake/3.11.4/bin/cmake

# The command to remove a file.
RM = /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/CMake/3.11.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build

# Utility rule file for pch_Generate_opencv_videoio.

# Include the progress variables for this target.
include modules/videoio/CMakeFiles/pch_Generate_opencv_videoio.dir/progress.make

modules/videoio/CMakeFiles/pch_Generate_opencv_videoio: modules/videoio/precomp.hpp.gch/opencv_videoio_Release.gch


modules/videoio/precomp.hpp.gch/opencv_videoio_Release.gch: ../modules/videoio/src/precomp.hpp
modules/videoio/precomp.hpp.gch/opencv_videoio_Release.gch: modules/videoio/precomp.hpp.command.sh
modules/videoio/precomp.hpp.gch/opencv_videoio_Release.gch: modules/videoio/precomp.hpp
modules/videoio/precomp.hpp.gch/opencv_videoio_Release.gch: lib/libopencv_videoio_pch_dephelp.a
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating precomp.hpp.gch/opencv_videoio_Release.gch"
	cd /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/modules/videoio && /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/CMake/3.11.4/bin/cmake -E make_directory /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/modules/videoio/precomp.hpp.gch
	cd /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/modules/videoio && chmod +x /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/modules/videoio/precomp.hpp.command.sh
	cd /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/modules/videoio && ./precomp.hpp.command.sh

modules/videoio/precomp.hpp: ../modules/videoio/src/precomp.hpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating precomp.hpp"
	cd /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/modules/videoio && /work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/CMake/3.11.4/bin/cmake -E copy_if_different /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/videoio/src/precomp.hpp /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/modules/videoio/precomp.hpp

pch_Generate_opencv_videoio: modules/videoio/CMakeFiles/pch_Generate_opencv_videoio
pch_Generate_opencv_videoio: modules/videoio/precomp.hpp.gch/opencv_videoio_Release.gch
pch_Generate_opencv_videoio: modules/videoio/precomp.hpp
pch_Generate_opencv_videoio: modules/videoio/CMakeFiles/pch_Generate_opencv_videoio.dir/build.make

.PHONY : pch_Generate_opencv_videoio

# Rule to build all files generated by this target.
modules/videoio/CMakeFiles/pch_Generate_opencv_videoio.dir/build: pch_Generate_opencv_videoio

.PHONY : modules/videoio/CMakeFiles/pch_Generate_opencv_videoio.dir/build

modules/videoio/CMakeFiles/pch_Generate_opencv_videoio.dir/clean:
	cd /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/modules/videoio && $(CMAKE_COMMAND) -P CMakeFiles/pch_Generate_opencv_videoio.dir/cmake_clean.cmake
.PHONY : modules/videoio/CMakeFiles/pch_Generate_opencv_videoio.dir/clean

modules/videoio/CMakeFiles/pch_Generate_opencv_videoio.dir/depend:
	cd /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7 /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/videoio /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/modules/videoio /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/modules/videoio/CMakeFiles/pch_Generate_opencv_videoio.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : modules/videoio/CMakeFiles/pch_Generate_opencv_videoio.dir/depend

