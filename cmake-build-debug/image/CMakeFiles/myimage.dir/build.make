# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.6

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "D:\clion\CLion 2016.3.5\bin\cmake\bin\cmake.exe"

# The command to remove a file.
RM = "D:\clion\CLion 2016.3.5\bin\cmake\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\dllibrary

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\dllibrary\cmake-build-debug

# Include any dependencies generated for this target.
include image/CMakeFiles/myimage.dir/depend.make

# Include the progress variables for this target.
include image/CMakeFiles/myimage.dir/progress.make

# Include the compile flags for this target's objects.
include image/CMakeFiles/myimage.dir/flags.make

image/CMakeFiles/myimage.dir/image.c.obj: image/CMakeFiles/myimage.dir/flags.make
image/CMakeFiles/myimage.dir/image.c.obj: image/CMakeFiles/myimage.dir/includes_C.rsp
image/CMakeFiles/myimage.dir/image.c.obj: ../image/image.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object image/CMakeFiles/myimage.dir/image.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\image && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\myimage.dir\image.c.obj   -c D:\dllibrary\image\image.c

image/CMakeFiles/myimage.dir/image.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myimage.dir/image.c.i"
	cd /d D:\dllibrary\cmake-build-debug\image && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\image\image.c > CMakeFiles\myimage.dir\image.c.i

image/CMakeFiles/myimage.dir/image.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myimage.dir/image.c.s"
	cd /d D:\dllibrary\cmake-build-debug\image && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\image\image.c -o CMakeFiles\myimage.dir\image.c.s

image/CMakeFiles/myimage.dir/image.c.obj.requires:

.PHONY : image/CMakeFiles/myimage.dir/image.c.obj.requires

image/CMakeFiles/myimage.dir/image.c.obj.provides: image/CMakeFiles/myimage.dir/image.c.obj.requires
	$(MAKE) -f image\CMakeFiles\myimage.dir\build.make image/CMakeFiles/myimage.dir/image.c.obj.provides.build
.PHONY : image/CMakeFiles/myimage.dir/image.c.obj.provides

image/CMakeFiles/myimage.dir/image.c.obj.provides.build: image/CMakeFiles/myimage.dir/image.c.obj


image/CMakeFiles/myimage.dir/matrix.c.obj: image/CMakeFiles/myimage.dir/flags.make
image/CMakeFiles/myimage.dir/matrix.c.obj: image/CMakeFiles/myimage.dir/includes_C.rsp
image/CMakeFiles/myimage.dir/matrix.c.obj: ../image/matrix.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object image/CMakeFiles/myimage.dir/matrix.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\image && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\myimage.dir\matrix.c.obj   -c D:\dllibrary\image\matrix.c

image/CMakeFiles/myimage.dir/matrix.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myimage.dir/matrix.c.i"
	cd /d D:\dllibrary\cmake-build-debug\image && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\image\matrix.c > CMakeFiles\myimage.dir\matrix.c.i

image/CMakeFiles/myimage.dir/matrix.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myimage.dir/matrix.c.s"
	cd /d D:\dllibrary\cmake-build-debug\image && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\image\matrix.c -o CMakeFiles\myimage.dir\matrix.c.s

image/CMakeFiles/myimage.dir/matrix.c.obj.requires:

.PHONY : image/CMakeFiles/myimage.dir/matrix.c.obj.requires

image/CMakeFiles/myimage.dir/matrix.c.obj.provides: image/CMakeFiles/myimage.dir/matrix.c.obj.requires
	$(MAKE) -f image\CMakeFiles\myimage.dir\build.make image/CMakeFiles/myimage.dir/matrix.c.obj.provides.build
.PHONY : image/CMakeFiles/myimage.dir/matrix.c.obj.provides

image/CMakeFiles/myimage.dir/matrix.c.obj.provides.build: image/CMakeFiles/myimage.dir/matrix.c.obj


image/CMakeFiles/myimage.dir/panorama_image.c.obj: image/CMakeFiles/myimage.dir/flags.make
image/CMakeFiles/myimage.dir/panorama_image.c.obj: image/CMakeFiles/myimage.dir/includes_C.rsp
image/CMakeFiles/myimage.dir/panorama_image.c.obj: ../image/panorama_image.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object image/CMakeFiles/myimage.dir/panorama_image.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\image && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\myimage.dir\panorama_image.c.obj   -c D:\dllibrary\image\panorama_image.c

image/CMakeFiles/myimage.dir/panorama_image.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myimage.dir/panorama_image.c.i"
	cd /d D:\dllibrary\cmake-build-debug\image && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\image\panorama_image.c > CMakeFiles\myimage.dir\panorama_image.c.i

image/CMakeFiles/myimage.dir/panorama_image.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myimage.dir/panorama_image.c.s"
	cd /d D:\dllibrary\cmake-build-debug\image && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\image\panorama_image.c -o CMakeFiles\myimage.dir\panorama_image.c.s

image/CMakeFiles/myimage.dir/panorama_image.c.obj.requires:

.PHONY : image/CMakeFiles/myimage.dir/panorama_image.c.obj.requires

image/CMakeFiles/myimage.dir/panorama_image.c.obj.provides: image/CMakeFiles/myimage.dir/panorama_image.c.obj.requires
	$(MAKE) -f image\CMakeFiles\myimage.dir\build.make image/CMakeFiles/myimage.dir/panorama_image.c.obj.provides.build
.PHONY : image/CMakeFiles/myimage.dir/panorama_image.c.obj.provides

image/CMakeFiles/myimage.dir/panorama_image.c.obj.provides.build: image/CMakeFiles/myimage.dir/panorama_image.c.obj


image/CMakeFiles/myimage.dir/spatial_filter.c.obj: image/CMakeFiles/myimage.dir/flags.make
image/CMakeFiles/myimage.dir/spatial_filter.c.obj: image/CMakeFiles/myimage.dir/includes_C.rsp
image/CMakeFiles/myimage.dir/spatial_filter.c.obj: ../image/spatial_filter.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object image/CMakeFiles/myimage.dir/spatial_filter.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\image && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\myimage.dir\spatial_filter.c.obj   -c D:\dllibrary\image\spatial_filter.c

image/CMakeFiles/myimage.dir/spatial_filter.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myimage.dir/spatial_filter.c.i"
	cd /d D:\dllibrary\cmake-build-debug\image && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\image\spatial_filter.c > CMakeFiles\myimage.dir\spatial_filter.c.i

image/CMakeFiles/myimage.dir/spatial_filter.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myimage.dir/spatial_filter.c.s"
	cd /d D:\dllibrary\cmake-build-debug\image && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\image\spatial_filter.c -o CMakeFiles\myimage.dir\spatial_filter.c.s

image/CMakeFiles/myimage.dir/spatial_filter.c.obj.requires:

.PHONY : image/CMakeFiles/myimage.dir/spatial_filter.c.obj.requires

image/CMakeFiles/myimage.dir/spatial_filter.c.obj.provides: image/CMakeFiles/myimage.dir/spatial_filter.c.obj.requires
	$(MAKE) -f image\CMakeFiles\myimage.dir\build.make image/CMakeFiles/myimage.dir/spatial_filter.c.obj.provides.build
.PHONY : image/CMakeFiles/myimage.dir/spatial_filter.c.obj.provides

image/CMakeFiles/myimage.dir/spatial_filter.c.obj.provides.build: image/CMakeFiles/myimage.dir/spatial_filter.c.obj


# Object files for target myimage
myimage_OBJECTS = \
"CMakeFiles/myimage.dir/image.c.obj" \
"CMakeFiles/myimage.dir/matrix.c.obj" \
"CMakeFiles/myimage.dir/panorama_image.c.obj" \
"CMakeFiles/myimage.dir/spatial_filter.c.obj"

# External object files for target myimage
myimage_EXTERNAL_OBJECTS =

image/myimage.exe: image/CMakeFiles/myimage.dir/image.c.obj
image/myimage.exe: image/CMakeFiles/myimage.dir/matrix.c.obj
image/myimage.exe: image/CMakeFiles/myimage.dir/panorama_image.c.obj
image/myimage.exe: image/CMakeFiles/myimage.dir/spatial_filter.c.obj
image/myimage.exe: image/CMakeFiles/myimage.dir/build.make
image/myimage.exe: image/CMakeFiles/myimage.dir/linklibs.rsp
image/myimage.exe: image/CMakeFiles/myimage.dir/objects1.rsp
image/myimage.exe: image/CMakeFiles/myimage.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking C executable myimage.exe"
	cd /d D:\dllibrary\cmake-build-debug\image && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\myimage.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
image/CMakeFiles/myimage.dir/build: image/myimage.exe

.PHONY : image/CMakeFiles/myimage.dir/build

image/CMakeFiles/myimage.dir/requires: image/CMakeFiles/myimage.dir/image.c.obj.requires
image/CMakeFiles/myimage.dir/requires: image/CMakeFiles/myimage.dir/matrix.c.obj.requires
image/CMakeFiles/myimage.dir/requires: image/CMakeFiles/myimage.dir/panorama_image.c.obj.requires
image/CMakeFiles/myimage.dir/requires: image/CMakeFiles/myimage.dir/spatial_filter.c.obj.requires

.PHONY : image/CMakeFiles/myimage.dir/requires

image/CMakeFiles/myimage.dir/clean:
	cd /d D:\dllibrary\cmake-build-debug\image && $(CMAKE_COMMAND) -P CMakeFiles\myimage.dir\cmake_clean.cmake
.PHONY : image/CMakeFiles/myimage.dir/clean

image/CMakeFiles/myimage.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\dllibrary D:\dllibrary\image D:\dllibrary\cmake-build-debug D:\dllibrary\cmake-build-debug\image D:\dllibrary\cmake-build-debug\image\CMakeFiles\myimage.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : image/CMakeFiles/myimage.dir/depend
