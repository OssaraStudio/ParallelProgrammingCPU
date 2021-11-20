#!/bin/sh
if [ -n "$VERBOSE" ]; then
  tail -n1 $0
fi
/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Core/GCCcore/7.3.0/bin/c++ -O3 -DNDEBUG -DNDEBUG -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build" -isystem"/work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/imkl/2018.3.222/mkl/include" -isystem"/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/GLib/2.54.3/include/gio-unix-2.0" -isystem"/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/GLib/2.54.3/include/glib-2.0" -isystem"/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/GLib/2.54.3/lib/glib-2.0/include" -isystem"/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/PCRE/8.41/include" -isystem"/usr/include/gtk-3.0" -isystem"/usr/include/atk-1.0" -isystem"/usr/include/at-spi2-atk/2.0" -isystem"/usr/include/pango-1.0" -isystem"/usr/include/cairo" -isystem"/usr/include/gdk-pixbuf-2.0" -isystem"/usr/include/at-spi-2.0" -isystem"/usr/include/dbus-1.0" -isystem"/usr/lib64/dbus-1.0/include" -isystem"/usr/include/harfbuzz" -isystem"/usr/include/freetype2" -isystem"/usr/include/pixman-1" -isystem"/usr/include/libpng15" -isystem"/usr/include/libdrm" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build" -isystem"/work/irlin355_1/gratienj/local/eb_new/centos_7/software/MPI/GCC/7.3.0-2.30/impi/2018.3.222/imkl/2018.3.222/mkl/include" -isystem"/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/GLib/2.54.3/include/gio-unix-2.0" -isystem"/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/GLib/2.54.3/include/glib-2.0" -isystem"/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/GLib/2.54.3/lib/glib-2.0/include" -isystem"/work/irlin355_1/gratienj/local/eb_new/centos_7/software/Compiler/GCCcore/7.3.0/PCRE/8.41/include" -isystem"/usr/include/gtk-3.0" -isystem"/usr/include/atk-1.0" -isystem"/usr/include/at-spi2-atk/2.0" -isystem"/usr/include/pango-1.0" -isystem"/usr/include/cairo" -isystem"/usr/include/gdk-pixbuf-2.0" -isystem"/usr/include/at-spi-2.0" -isystem"/usr/include/dbus-1.0" -isystem"/usr/lib64/dbus-1.0/include" -isystem"/usr/include/harfbuzz" -isystem"/usr/include/freetype2" -isystem"/usr/include/pixman-1" -isystem"/usr/include/libpng15" -isystem"/usr/include/libdrm" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/ts/include" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/highgui/include" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/imgcodecs/include" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/videoio/include" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/core/include" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/imgproc/include" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/imgcodecs/include" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/videoio/include" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/core/include" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/imgproc/include" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/imgcodecs/include" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/videoio/include" -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/highgui/include" -fsigned-char -W -Wall -Werror=return-type -Werror=non-virtual-dtor -Werror=address -Werror=sequence-point -Wformat -Werror=format-security -Wmissing-declarations -Wundef -Winit-self -Wpointer-arith -Wshadow -Wsign-promo -Wuninitialized -Winit-self -Wno-delete-non-virtual-dtor -Wno-comment -Wimplicit-fallthrough=3 -Wno-strict-overflow -fdiagnostics-show-option -Wno-long-long -pthread -fomit-frame-pointer -ffunction-sections -fdata-sections -msse -msse2 -msse3 -fvisibility=hidden -fvisibility-inlines-hidden -Wno-deprecated-declarations -I"/work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/modules/highgui/test" -fPIE -x c++-header -o /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/modules/highgui/test_precomp.hpp.gch/opencv_test_highgui_Release.gch /work/irlin355_1/gratienj/ParallelProgrammingCourse/ParallelProgrammingTP/contrib/opencv-3.4.7/build/modules/highgui/test_precomp.hpp '-D__OPENCV_TESTS=1' '-D__OPENCV_BUILD=1' '-D_USE_MATH_DEFINES' '-D__STDC_CONSTANT_MACROS' '-D__STDC_LIMIT_MACROS' '-D__STDC_FORMAT_MACROS' '-DHAVE_WEBP'
