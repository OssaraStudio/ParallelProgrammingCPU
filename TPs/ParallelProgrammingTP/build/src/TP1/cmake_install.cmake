# Install script for directory: /work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/src/TP1

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELEASE")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_stdthread.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_stdthread.exe")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_stdthread.exe"
         RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_stdthread.exe")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin" TYPE EXECUTABLE FILES "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/src/TP1/helloworld_stdthread.exe")
  if(EXISTS "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_stdthread.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_stdthread.exe")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_stdthread.exe"
         OLD_RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/contrib/opencv/3.4.7/lib64:"
         NEW_RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/binutils/2.30/bin/strip" "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_stdthread.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_openmp.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_openmp.exe")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_openmp.exe"
         RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_openmp.exe")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin" TYPE EXECUTABLE FILES "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/src/TP1/helloworld_openmp.exe")
  if(EXISTS "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_openmp.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_openmp.exe")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_openmp.exe"
         OLD_RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/contrib/opencv/3.4.7/lib64:"
         NEW_RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/binutils/2.30/bin/strip" "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_openmp.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_task.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_task.exe")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_task.exe"
         RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_task.exe")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin" TYPE EXECUTABLE FILES "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/src/TP1/helloworld_task.exe")
  if(EXISTS "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_task.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_task.exe")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_task.exe"
         OLD_RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/contrib/opencv/3.4.7/lib64:"
         NEW_RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/binutils/2.30/bin/strip" "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_task.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_tbb.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_tbb.exe")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_tbb.exe"
         RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_tbb.exe")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin" TYPE EXECUTABLE FILES "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/src/TP1/helloworld_tbb.exe")
  if(EXISTS "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_tbb.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_tbb.exe")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_tbb.exe"
         OLD_RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/contrib/opencv/3.4.7/lib64:"
         NEW_RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/binutils/2.30/bin/strip" "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_tbb.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_mpi.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_mpi.exe")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_mpi.exe"
         RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_mpi.exe")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin" TYPE EXECUTABLE FILES "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/build/src/TP1/helloworld_mpi.exe")
  if(EXISTS "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_mpi.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_mpi.exe")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_mpi.exe"
         OLD_RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/contrib/opencv/3.4.7/lib64:"
         NEW_RPATH "/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/binutils/2.30/bin/strip" "$ENV{DESTDIR}/work/irlin355_1/gratienj/ParallelProgrammingCourse/GIT/ParallelProgrammingCourse/TPs/ParallelProgrammingTP/bin/helloworld_mpi.exe")
    endif()
  endif()
endif()

