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
include utils/CMakeFiles/utils.dir/depend.make

# Include the progress variables for this target.
include utils/CMakeFiles/utils.dir/progress.make

# Include the compile flags for this target's objects.
include utils/CMakeFiles/utils.dir/flags.make

utils/CMakeFiles/utils.dir/utils.c.obj: utils/CMakeFiles/utils.dir/flags.make
utils/CMakeFiles/utils.dir/utils.c.obj: utils/CMakeFiles/utils.dir/includes_C.rsp
utils/CMakeFiles/utils.dir/utils.c.obj: ../utils/utils.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object utils/CMakeFiles/utils.dir/utils.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\utils.dir\utils.c.obj   -c D:\dllibrary\utils\utils.c

utils/CMakeFiles/utils.dir/utils.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/utils.dir/utils.c.i"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\utils\utils.c > CMakeFiles\utils.dir\utils.c.i

utils/CMakeFiles/utils.dir/utils.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/utils.dir/utils.c.s"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\utils\utils.c -o CMakeFiles\utils.dir\utils.c.s

utils/CMakeFiles/utils.dir/utils.c.obj.requires:

.PHONY : utils/CMakeFiles/utils.dir/utils.c.obj.requires

utils/CMakeFiles/utils.dir/utils.c.obj.provides: utils/CMakeFiles/utils.dir/utils.c.obj.requires
	$(MAKE) -f utils\CMakeFiles\utils.dir\build.make utils/CMakeFiles/utils.dir/utils.c.obj.provides.build
.PHONY : utils/CMakeFiles/utils.dir/utils.c.obj.provides

utils/CMakeFiles/utils.dir/utils.c.obj.provides.build: utils/CMakeFiles/utils.dir/utils.c.obj


utils/CMakeFiles/utils.dir/blas.c.obj: utils/CMakeFiles/utils.dir/flags.make
utils/CMakeFiles/utils.dir/blas.c.obj: utils/CMakeFiles/utils.dir/includes_C.rsp
utils/CMakeFiles/utils.dir/blas.c.obj: ../utils/blas.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object utils/CMakeFiles/utils.dir/blas.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\utils.dir\blas.c.obj   -c D:\dllibrary\utils\blas.c

utils/CMakeFiles/utils.dir/blas.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/utils.dir/blas.c.i"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\utils\blas.c > CMakeFiles\utils.dir\blas.c.i

utils/CMakeFiles/utils.dir/blas.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/utils.dir/blas.c.s"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\utils\blas.c -o CMakeFiles\utils.dir\blas.c.s

utils/CMakeFiles/utils.dir/blas.c.obj.requires:

.PHONY : utils/CMakeFiles/utils.dir/blas.c.obj.requires

utils/CMakeFiles/utils.dir/blas.c.obj.provides: utils/CMakeFiles/utils.dir/blas.c.obj.requires
	$(MAKE) -f utils\CMakeFiles\utils.dir\build.make utils/CMakeFiles/utils.dir/blas.c.obj.provides.build
.PHONY : utils/CMakeFiles/utils.dir/blas.c.obj.provides

utils/CMakeFiles/utils.dir/blas.c.obj.provides.build: utils/CMakeFiles/utils.dir/blas.c.obj


utils/CMakeFiles/utils.dir/im2col.c.obj: utils/CMakeFiles/utils.dir/flags.make
utils/CMakeFiles/utils.dir/im2col.c.obj: utils/CMakeFiles/utils.dir/includes_C.rsp
utils/CMakeFiles/utils.dir/im2col.c.obj: ../utils/im2col.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object utils/CMakeFiles/utils.dir/im2col.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\utils.dir\im2col.c.obj   -c D:\dllibrary\utils\im2col.c

utils/CMakeFiles/utils.dir/im2col.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/utils.dir/im2col.c.i"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\utils\im2col.c > CMakeFiles\utils.dir\im2col.c.i

utils/CMakeFiles/utils.dir/im2col.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/utils.dir/im2col.c.s"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\utils\im2col.c -o CMakeFiles\utils.dir\im2col.c.s

utils/CMakeFiles/utils.dir/im2col.c.obj.requires:

.PHONY : utils/CMakeFiles/utils.dir/im2col.c.obj.requires

utils/CMakeFiles/utils.dir/im2col.c.obj.provides: utils/CMakeFiles/utils.dir/im2col.c.obj.requires
	$(MAKE) -f utils\CMakeFiles\utils.dir\build.make utils/CMakeFiles/utils.dir/im2col.c.obj.provides.build
.PHONY : utils/CMakeFiles/utils.dir/im2col.c.obj.provides

utils/CMakeFiles/utils.dir/im2col.c.obj.provides.build: utils/CMakeFiles/utils.dir/im2col.c.obj


utils/CMakeFiles/utils.dir/gemm.c.obj: utils/CMakeFiles/utils.dir/flags.make
utils/CMakeFiles/utils.dir/gemm.c.obj: utils/CMakeFiles/utils.dir/includes_C.rsp
utils/CMakeFiles/utils.dir/gemm.c.obj: ../utils/gemm.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object utils/CMakeFiles/utils.dir/gemm.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\utils.dir\gemm.c.obj   -c D:\dllibrary\utils\gemm.c

utils/CMakeFiles/utils.dir/gemm.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/utils.dir/gemm.c.i"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\utils\gemm.c > CMakeFiles\utils.dir\gemm.c.i

utils/CMakeFiles/utils.dir/gemm.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/utils.dir/gemm.c.s"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\utils\gemm.c -o CMakeFiles\utils.dir\gemm.c.s

utils/CMakeFiles/utils.dir/gemm.c.obj.requires:

.PHONY : utils/CMakeFiles/utils.dir/gemm.c.obj.requires

utils/CMakeFiles/utils.dir/gemm.c.obj.provides: utils/CMakeFiles/utils.dir/gemm.c.obj.requires
	$(MAKE) -f utils\CMakeFiles\utils.dir\build.make utils/CMakeFiles/utils.dir/gemm.c.obj.provides.build
.PHONY : utils/CMakeFiles/utils.dir/gemm.c.obj.provides

utils/CMakeFiles/utils.dir/gemm.c.obj.provides.build: utils/CMakeFiles/utils.dir/gemm.c.obj


utils/CMakeFiles/utils.dir/list.c.obj: utils/CMakeFiles/utils.dir/flags.make
utils/CMakeFiles/utils.dir/list.c.obj: utils/CMakeFiles/utils.dir/includes_C.rsp
utils/CMakeFiles/utils.dir/list.c.obj: ../utils/list.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object utils/CMakeFiles/utils.dir/list.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\utils.dir\list.c.obj   -c D:\dllibrary\utils\list.c

utils/CMakeFiles/utils.dir/list.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/utils.dir/list.c.i"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\utils\list.c > CMakeFiles\utils.dir\list.c.i

utils/CMakeFiles/utils.dir/list.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/utils.dir/list.c.s"
	cd /d D:\dllibrary\cmake-build-debug\utils && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\utils\list.c -o CMakeFiles\utils.dir\list.c.s

utils/CMakeFiles/utils.dir/list.c.obj.requires:

.PHONY : utils/CMakeFiles/utils.dir/list.c.obj.requires

utils/CMakeFiles/utils.dir/list.c.obj.provides: utils/CMakeFiles/utils.dir/list.c.obj.requires
	$(MAKE) -f utils\CMakeFiles\utils.dir\build.make utils/CMakeFiles/utils.dir/list.c.obj.provides.build
.PHONY : utils/CMakeFiles/utils.dir/list.c.obj.provides

utils/CMakeFiles/utils.dir/list.c.obj.provides.build: utils/CMakeFiles/utils.dir/list.c.obj


# Object files for target utils
utils_OBJECTS = \
"CMakeFiles/utils.dir/utils.c.obj" \
"CMakeFiles/utils.dir/blas.c.obj" \
"CMakeFiles/utils.dir/im2col.c.obj" \
"CMakeFiles/utils.dir/gemm.c.obj" \
"CMakeFiles/utils.dir/list.c.obj"

# External object files for target utils
utils_EXTERNAL_OBJECTS =

utils/libutils.a: utils/CMakeFiles/utils.dir/utils.c.obj
utils/libutils.a: utils/CMakeFiles/utils.dir/blas.c.obj
utils/libutils.a: utils/CMakeFiles/utils.dir/im2col.c.obj
utils/libutils.a: utils/CMakeFiles/utils.dir/gemm.c.obj
utils/libutils.a: utils/CMakeFiles/utils.dir/list.c.obj
utils/libutils.a: utils/CMakeFiles/utils.dir/build.make
utils/libutils.a: utils/CMakeFiles/utils.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking C static library libutils.a"
	cd /d D:\dllibrary\cmake-build-debug\utils && $(CMAKE_COMMAND) -P CMakeFiles\utils.dir\cmake_clean_target.cmake
	cd /d D:\dllibrary\cmake-build-debug\utils && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\utils.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
utils/CMakeFiles/utils.dir/build: utils/libutils.a

.PHONY : utils/CMakeFiles/utils.dir/build

utils/CMakeFiles/utils.dir/requires: utils/CMakeFiles/utils.dir/utils.c.obj.requires
utils/CMakeFiles/utils.dir/requires: utils/CMakeFiles/utils.dir/blas.c.obj.requires
utils/CMakeFiles/utils.dir/requires: utils/CMakeFiles/utils.dir/im2col.c.obj.requires
utils/CMakeFiles/utils.dir/requires: utils/CMakeFiles/utils.dir/gemm.c.obj.requires
utils/CMakeFiles/utils.dir/requires: utils/CMakeFiles/utils.dir/list.c.obj.requires

.PHONY : utils/CMakeFiles/utils.dir/requires

utils/CMakeFiles/utils.dir/clean:
	cd /d D:\dllibrary\cmake-build-debug\utils && $(CMAKE_COMMAND) -P CMakeFiles\utils.dir\cmake_clean.cmake
.PHONY : utils/CMakeFiles/utils.dir/clean

utils/CMakeFiles/utils.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\dllibrary D:\dllibrary\utils D:\dllibrary\cmake-build-debug D:\dllibrary\cmake-build-debug\utils D:\dllibrary\cmake-build-debug\utils\CMakeFiles\utils.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : utils/CMakeFiles/utils.dir/depend

