# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.20.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.20.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9/build

# Include any dependencies generated for this target.
include CMakeFiles/TestProject.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/TestProject.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/TestProject.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TestProject.dir/flags.make

CMakeFiles/TestProject.dir/src/TestProject.cpp.o: CMakeFiles/TestProject.dir/flags.make
CMakeFiles/TestProject.dir/src/TestProject.cpp.o: ../src/TestProject.cpp
CMakeFiles/TestProject.dir/src/TestProject.cpp.o: CMakeFiles/TestProject.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/TestProject.dir/src/TestProject.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/TestProject.dir/src/TestProject.cpp.o -MF CMakeFiles/TestProject.dir/src/TestProject.cpp.o.d -o CMakeFiles/TestProject.dir/src/TestProject.cpp.o -c /Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9/src/TestProject.cpp

CMakeFiles/TestProject.dir/src/TestProject.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestProject.dir/src/TestProject.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9/src/TestProject.cpp > CMakeFiles/TestProject.dir/src/TestProject.cpp.i

CMakeFiles/TestProject.dir/src/TestProject.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestProject.dir/src/TestProject.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9/src/TestProject.cpp -o CMakeFiles/TestProject.dir/src/TestProject.cpp.s

# Object files for target TestProject
TestProject_OBJECTS = \
"CMakeFiles/TestProject.dir/src/TestProject.cpp.o"

# External object files for target TestProject
TestProject_EXTERNAL_OBJECTS =

TestProject: CMakeFiles/TestProject.dir/src/TestProject.cpp.o
TestProject: CMakeFiles/TestProject.dir/build.make
TestProject: CMakeFiles/TestProject.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable TestProject"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestProject.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TestProject.dir/build: TestProject
.PHONY : CMakeFiles/TestProject.dir/build

CMakeFiles/TestProject.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/TestProject.dir/cmake_clean.cmake
.PHONY : CMakeFiles/TestProject.dir/clean

CMakeFiles/TestProject.dir/depend:
	cd /Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9 /Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9 /Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9/build /Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9/build /Users/macbookpro/Desktop/ENSIMAG/3A/hedging-derivatives/Equipe_9/build/CMakeFiles/TestProject.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/TestProject.dir/depend

