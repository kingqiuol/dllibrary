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
include src/CMakeFiles/src.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/src.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/src.dir/flags.make

src/CMakeFiles/src.dir/convolutional_layer.c.obj: src/CMakeFiles/src.dir/flags.make
src/CMakeFiles/src.dir/convolutional_layer.c.obj: src/CMakeFiles/src.dir/includes_C.rsp
src/CMakeFiles/src.dir/convolutional_layer.c.obj: ../src/convolutional_layer.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/CMakeFiles/src.dir/convolutional_layer.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\src.dir\convolutional_layer.c.obj   -c D:\dllibrary\src\convolutional_layer.c

src/CMakeFiles/src.dir/convolutional_layer.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/src.dir/convolutional_layer.c.i"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\src\convolutional_layer.c > CMakeFiles\src.dir\convolutional_layer.c.i

src/CMakeFiles/src.dir/convolutional_layer.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/src.dir/convolutional_layer.c.s"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\src\convolutional_layer.c -o CMakeFiles\src.dir\convolutional_layer.c.s

src/CMakeFiles/src.dir/convolutional_layer.c.obj.requires:

.PHONY : src/CMakeFiles/src.dir/convolutional_layer.c.obj.requires

src/CMakeFiles/src.dir/convolutional_layer.c.obj.provides: src/CMakeFiles/src.dir/convolutional_layer.c.obj.requires
	$(MAKE) -f src\CMakeFiles\src.dir\build.make src/CMakeFiles/src.dir/convolutional_layer.c.obj.provides.build
.PHONY : src/CMakeFiles/src.dir/convolutional_layer.c.obj.provides

src/CMakeFiles/src.dir/convolutional_layer.c.obj.provides.build: src/CMakeFiles/src.dir/convolutional_layer.c.obj


src/CMakeFiles/src.dir/layer.c.obj: src/CMakeFiles/src.dir/flags.make
src/CMakeFiles/src.dir/layer.c.obj: src/CMakeFiles/src.dir/includes_C.rsp
src/CMakeFiles/src.dir/layer.c.obj: ../src/layer.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/CMakeFiles/src.dir/layer.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\src.dir\layer.c.obj   -c D:\dllibrary\src\layer.c

src/CMakeFiles/src.dir/layer.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/src.dir/layer.c.i"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\src\layer.c > CMakeFiles\src.dir\layer.c.i

src/CMakeFiles/src.dir/layer.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/src.dir/layer.c.s"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\src\layer.c -o CMakeFiles\src.dir\layer.c.s

src/CMakeFiles/src.dir/layer.c.obj.requires:

.PHONY : src/CMakeFiles/src.dir/layer.c.obj.requires

src/CMakeFiles/src.dir/layer.c.obj.provides: src/CMakeFiles/src.dir/layer.c.obj.requires
	$(MAKE) -f src\CMakeFiles\src.dir\build.make src/CMakeFiles/src.dir/layer.c.obj.provides.build
.PHONY : src/CMakeFiles/src.dir/layer.c.obj.provides

src/CMakeFiles/src.dir/layer.c.obj.provides.build: src/CMakeFiles/src.dir/layer.c.obj


src/CMakeFiles/src.dir/batchnorm_layer.c.obj: src/CMakeFiles/src.dir/flags.make
src/CMakeFiles/src.dir/batchnorm_layer.c.obj: src/CMakeFiles/src.dir/includes_C.rsp
src/CMakeFiles/src.dir/batchnorm_layer.c.obj: ../src/batchnorm_layer.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object src/CMakeFiles/src.dir/batchnorm_layer.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\src.dir\batchnorm_layer.c.obj   -c D:\dllibrary\src\batchnorm_layer.c

src/CMakeFiles/src.dir/batchnorm_layer.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/src.dir/batchnorm_layer.c.i"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\src\batchnorm_layer.c > CMakeFiles\src.dir\batchnorm_layer.c.i

src/CMakeFiles/src.dir/batchnorm_layer.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/src.dir/batchnorm_layer.c.s"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\src\batchnorm_layer.c -o CMakeFiles\src.dir\batchnorm_layer.c.s

src/CMakeFiles/src.dir/batchnorm_layer.c.obj.requires:

.PHONY : src/CMakeFiles/src.dir/batchnorm_layer.c.obj.requires

src/CMakeFiles/src.dir/batchnorm_layer.c.obj.provides: src/CMakeFiles/src.dir/batchnorm_layer.c.obj.requires
	$(MAKE) -f src\CMakeFiles\src.dir\build.make src/CMakeFiles/src.dir/batchnorm_layer.c.obj.provides.build
.PHONY : src/CMakeFiles/src.dir/batchnorm_layer.c.obj.provides

src/CMakeFiles/src.dir/batchnorm_layer.c.obj.provides.build: src/CMakeFiles/src.dir/batchnorm_layer.c.obj


src/CMakeFiles/src.dir/network.c.obj: src/CMakeFiles/src.dir/flags.make
src/CMakeFiles/src.dir/network.c.obj: src/CMakeFiles/src.dir/includes_C.rsp
src/CMakeFiles/src.dir/network.c.obj: ../src/network.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object src/CMakeFiles/src.dir/network.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\src.dir\network.c.obj   -c D:\dllibrary\src\network.c

src/CMakeFiles/src.dir/network.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/src.dir/network.c.i"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\src\network.c > CMakeFiles\src.dir\network.c.i

src/CMakeFiles/src.dir/network.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/src.dir/network.c.s"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\src\network.c -o CMakeFiles\src.dir\network.c.s

src/CMakeFiles/src.dir/network.c.obj.requires:

.PHONY : src/CMakeFiles/src.dir/network.c.obj.requires

src/CMakeFiles/src.dir/network.c.obj.provides: src/CMakeFiles/src.dir/network.c.obj.requires
	$(MAKE) -f src\CMakeFiles\src.dir\build.make src/CMakeFiles/src.dir/network.c.obj.provides.build
.PHONY : src/CMakeFiles/src.dir/network.c.obj.provides

src/CMakeFiles/src.dir/network.c.obj.provides.build: src/CMakeFiles/src.dir/network.c.obj


src/CMakeFiles/src.dir/activations.c.obj: src/CMakeFiles/src.dir/flags.make
src/CMakeFiles/src.dir/activations.c.obj: src/CMakeFiles/src.dir/includes_C.rsp
src/CMakeFiles/src.dir/activations.c.obj: ../src/activations.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object src/CMakeFiles/src.dir/activations.c.obj"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\src.dir\activations.c.obj   -c D:\dllibrary\src\activations.c

src/CMakeFiles/src.dir/activations.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/src.dir/activations.c.i"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\dllibrary\src\activations.c > CMakeFiles\src.dir\activations.c.i

src/CMakeFiles/src.dir/activations.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/src.dir/activations.c.s"
	cd /d D:\dllibrary\cmake-build-debug\src && D:\vim\mingw\bin\gcc.exe  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\dllibrary\src\activations.c -o CMakeFiles\src.dir\activations.c.s

src/CMakeFiles/src.dir/activations.c.obj.requires:

.PHONY : src/CMakeFiles/src.dir/activations.c.obj.requires

src/CMakeFiles/src.dir/activations.c.obj.provides: src/CMakeFiles/src.dir/activations.c.obj.requires
	$(MAKE) -f src\CMakeFiles\src.dir\build.make src/CMakeFiles/src.dir/activations.c.obj.provides.build
.PHONY : src/CMakeFiles/src.dir/activations.c.obj.provides

src/CMakeFiles/src.dir/activations.c.obj.provides.build: src/CMakeFiles/src.dir/activations.c.obj


# Object files for target src
src_OBJECTS = \
"CMakeFiles/src.dir/convolutional_layer.c.obj" \
"CMakeFiles/src.dir/layer.c.obj" \
"CMakeFiles/src.dir/batchnorm_layer.c.obj" \
"CMakeFiles/src.dir/network.c.obj" \
"CMakeFiles/src.dir/activations.c.obj"

# External object files for target src
src_EXTERNAL_OBJECTS =

src/libsrc.a: src/CMakeFiles/src.dir/convolutional_layer.c.obj
src/libsrc.a: src/CMakeFiles/src.dir/layer.c.obj
src/libsrc.a: src/CMakeFiles/src.dir/batchnorm_layer.c.obj
src/libsrc.a: src/CMakeFiles/src.dir/network.c.obj
src/libsrc.a: src/CMakeFiles/src.dir/activations.c.obj
src/libsrc.a: src/CMakeFiles/src.dir/build.make
src/libsrc.a: src/CMakeFiles/src.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\dllibrary\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking C static library libsrc.a"
	cd /d D:\dllibrary\cmake-build-debug\src && $(CMAKE_COMMAND) -P CMakeFiles\src.dir\cmake_clean_target.cmake
	cd /d D:\dllibrary\cmake-build-debug\src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\src.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/src.dir/build: src/libsrc.a

.PHONY : src/CMakeFiles/src.dir/build

src/CMakeFiles/src.dir/requires: src/CMakeFiles/src.dir/convolutional_layer.c.obj.requires
src/CMakeFiles/src.dir/requires: src/CMakeFiles/src.dir/layer.c.obj.requires
src/CMakeFiles/src.dir/requires: src/CMakeFiles/src.dir/batchnorm_layer.c.obj.requires
src/CMakeFiles/src.dir/requires: src/CMakeFiles/src.dir/network.c.obj.requires
src/CMakeFiles/src.dir/requires: src/CMakeFiles/src.dir/activations.c.obj.requires

.PHONY : src/CMakeFiles/src.dir/requires

src/CMakeFiles/src.dir/clean:
	cd /d D:\dllibrary\cmake-build-debug\src && $(CMAKE_COMMAND) -P CMakeFiles\src.dir\cmake_clean.cmake
.PHONY : src/CMakeFiles/src.dir/clean

src/CMakeFiles/src.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\dllibrary D:\dllibrary\src D:\dllibrary\cmake-build-debug D:\dllibrary\cmake-build-debug\src D:\dllibrary\cmake-build-debug\src\CMakeFiles\src.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/src.dir/depend

