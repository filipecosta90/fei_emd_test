# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.8.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.8.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/filipeoliveira/Documents/fei_emd_test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/filipeoliveira/Documents/fei_emd_test/build

# Include any dependencies generated for this target.
include CMakeFiles/emd_test.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/emd_test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/emd_test.dir/flags.make

CMakeFiles/emd_test.dir/test_emd.cpp.o: CMakeFiles/emd_test.dir/flags.make
CMakeFiles/emd_test.dir/test_emd.cpp.o: ../test_emd.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/filipeoliveira/Documents/fei_emd_test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/emd_test.dir/test_emd.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/emd_test.dir/test_emd.cpp.o -c /Users/filipeoliveira/Documents/fei_emd_test/test_emd.cpp

CMakeFiles/emd_test.dir/test_emd.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/emd_test.dir/test_emd.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/filipeoliveira/Documents/fei_emd_test/test_emd.cpp > CMakeFiles/emd_test.dir/test_emd.cpp.i

CMakeFiles/emd_test.dir/test_emd.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/emd_test.dir/test_emd.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/filipeoliveira/Documents/fei_emd_test/test_emd.cpp -o CMakeFiles/emd_test.dir/test_emd.cpp.s

CMakeFiles/emd_test.dir/test_emd.cpp.o.requires:

.PHONY : CMakeFiles/emd_test.dir/test_emd.cpp.o.requires

CMakeFiles/emd_test.dir/test_emd.cpp.o.provides: CMakeFiles/emd_test.dir/test_emd.cpp.o.requires
	$(MAKE) -f CMakeFiles/emd_test.dir/build.make CMakeFiles/emd_test.dir/test_emd.cpp.o.provides.build
.PHONY : CMakeFiles/emd_test.dir/test_emd.cpp.o.provides

CMakeFiles/emd_test.dir/test_emd.cpp.o.provides.build: CMakeFiles/emd_test.dir/test_emd.cpp.o


# Object files for target emd_test
emd_test_OBJECTS = \
"CMakeFiles/emd_test.dir/test_emd.cpp.o"

# External object files for target emd_test
emd_test_EXTERNAL_OBJECTS =

emd_test: CMakeFiles/emd_test.dir/test_emd.cpp.o
emd_test: CMakeFiles/emd_test.dir/build.make
emd_test: /usr/local/Cellar/hdf5/1.10.1/lib/libhdf5_cpp.dylib
emd_test: /usr/local/Cellar/hdf5/1.10.1/lib/libhdf5.dylib
emd_test: /usr/local/opt/szip/lib/libsz.dylib
emd_test: /usr/lib/libz.dylib
emd_test: /usr/lib/libdl.dylib
emd_test: /usr/lib/libm.dylib
emd_test: CMakeFiles/emd_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/filipeoliveira/Documents/fei_emd_test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable emd_test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/emd_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/emd_test.dir/build: emd_test

.PHONY : CMakeFiles/emd_test.dir/build

CMakeFiles/emd_test.dir/requires: CMakeFiles/emd_test.dir/test_emd.cpp.o.requires

.PHONY : CMakeFiles/emd_test.dir/requires

CMakeFiles/emd_test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/emd_test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/emd_test.dir/clean

CMakeFiles/emd_test.dir/depend:
	cd /Users/filipeoliveira/Documents/fei_emd_test/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/filipeoliveira/Documents/fei_emd_test /Users/filipeoliveira/Documents/fei_emd_test /Users/filipeoliveira/Documents/fei_emd_test/build /Users/filipeoliveira/Documents/fei_emd_test/build /Users/filipeoliveira/Documents/fei_emd_test/build/CMakeFiles/emd_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/emd_test.dir/depend

