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
CMAKE_SOURCE_DIR = /home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution/build

# Include any dependencies generated for this target.
include CMakeFiles/PathTracing.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/PathTracing.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/PathTracing.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/PathTracing.dir/flags.make

CMakeFiles/PathTracing.dir/PathTracing.cpp.o: CMakeFiles/PathTracing.dir/flags.make
CMakeFiles/PathTracing.dir/PathTracing.cpp.o: ../PathTracing.cpp
CMakeFiles/PathTracing.dir/PathTracing.cpp.o: CMakeFiles/PathTracing.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/PathTracing.dir/PathTracing.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/PathTracing.dir/PathTracing.cpp.o -MF CMakeFiles/PathTracing.dir/PathTracing.cpp.o.d -o CMakeFiles/PathTracing.dir/PathTracing.cpp.o -c /home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution/PathTracing.cpp

CMakeFiles/PathTracing.dir/PathTracing.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PathTracing.dir/PathTracing.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution/PathTracing.cpp > CMakeFiles/PathTracing.dir/PathTracing.cpp.i

CMakeFiles/PathTracing.dir/PathTracing.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PathTracing.dir/PathTracing.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution/PathTracing.cpp -o CMakeFiles/PathTracing.dir/PathTracing.cpp.s

# Object files for target PathTracing
PathTracing_OBJECTS = \
"CMakeFiles/PathTracing.dir/PathTracing.cpp.o"

# External object files for target PathTracing
PathTracing_EXTERNAL_OBJECTS =

bin/PathTracing: CMakeFiles/PathTracing.dir/PathTracing.cpp.o
bin/PathTracing: CMakeFiles/PathTracing.dir/build.make
bin/PathTracing: lib/libtinyobjloader.a
bin/PathTracing: CMakeFiles/PathTracing.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bin/PathTracing"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PathTracing.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/PathTracing.dir/build: bin/PathTracing
.PHONY : CMakeFiles/PathTracing.dir/build

CMakeFiles/PathTracing.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/PathTracing.dir/cmake_clean.cmake
.PHONY : CMakeFiles/PathTracing.dir/clean

CMakeFiles/PathTracing.dir/depend:
	cd /home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution /home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution /home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution/build /home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution/build /home/moritz/computergraphics/02_Programming_Assignment_2/Submission/Solutions_PA2_Perschke_Mittermair/Solution/build/CMakeFiles/PathTracing.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/PathTracing.dir/depend
