# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/moritz/computergraphics/02_Programming_Assignment_2/Solution

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build

# Include any dependencies generated for this target.
include external/tinyobjloader/CMakeFiles/tinyobjloader.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include external/tinyobjloader/CMakeFiles/tinyobjloader.dir/compiler_depend.make

# Include the progress variables for this target.
include external/tinyobjloader/CMakeFiles/tinyobjloader.dir/progress.make

# Include the compile flags for this target's objects.
include external/tinyobjloader/CMakeFiles/tinyobjloader.dir/flags.make

external/tinyobjloader/CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.o: external/tinyobjloader/CMakeFiles/tinyobjloader.dir/flags.make
external/tinyobjloader/CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.o: ../external/tinyobjloader/tiny_obj_loader.cc
external/tinyobjloader/CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.o: external/tinyobjloader/CMakeFiles/tinyobjloader.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/tinyobjloader/CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.o"
	cd /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build/external/tinyobjloader && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT external/tinyobjloader/CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.o -MF CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.o.d -o CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.o -c /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/external/tinyobjloader/tiny_obj_loader.cc

external/tinyobjloader/CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.i"
	cd /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build/external/tinyobjloader && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/external/tinyobjloader/tiny_obj_loader.cc > CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.i

external/tinyobjloader/CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.s"
	cd /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build/external/tinyobjloader && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/external/tinyobjloader/tiny_obj_loader.cc -o CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.s

# Object files for target tinyobjloader
tinyobjloader_OBJECTS = \
"CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.o"

# External object files for target tinyobjloader
tinyobjloader_EXTERNAL_OBJECTS =

lib/libtinyobjloader.a: external/tinyobjloader/CMakeFiles/tinyobjloader.dir/tiny_obj_loader.cc.o
lib/libtinyobjloader.a: external/tinyobjloader/CMakeFiles/tinyobjloader.dir/build.make
lib/libtinyobjloader.a: external/tinyobjloader/CMakeFiles/tinyobjloader.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library ../../lib/libtinyobjloader.a"
	cd /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build/external/tinyobjloader && $(CMAKE_COMMAND) -P CMakeFiles/tinyobjloader.dir/cmake_clean_target.cmake
	cd /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build/external/tinyobjloader && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tinyobjloader.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/tinyobjloader/CMakeFiles/tinyobjloader.dir/build: lib/libtinyobjloader.a
.PHONY : external/tinyobjloader/CMakeFiles/tinyobjloader.dir/build

external/tinyobjloader/CMakeFiles/tinyobjloader.dir/clean:
	cd /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build/external/tinyobjloader && $(CMAKE_COMMAND) -P CMakeFiles/tinyobjloader.dir/cmake_clean.cmake
.PHONY : external/tinyobjloader/CMakeFiles/tinyobjloader.dir/clean

external/tinyobjloader/CMakeFiles/tinyobjloader.dir/depend:
	cd /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/moritz/computergraphics/02_Programming_Assignment_2/Solution /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/external/tinyobjloader /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build/external/tinyobjloader /home/moritz/computergraphics/02_Programming_Assignment_2/Solution/build/external/tinyobjloader/CMakeFiles/tinyobjloader.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/tinyobjloader/CMakeFiles/tinyobjloader.dir/depend

