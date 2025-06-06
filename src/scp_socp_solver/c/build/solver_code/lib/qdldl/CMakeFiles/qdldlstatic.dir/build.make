# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.31

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
CMAKE_COMMAND = /Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build

# Include any dependencies generated for this target.
include solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/compiler_depend.make

# Include the progress variables for this target.
include solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/progress.make

# Include the compile flags for this target's objects.
include solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/flags.make

solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/codegen:
.PHONY : solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/codegen

solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/src/qdldl.c.o: solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/flags.make
solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/src/qdldl.c.o: /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/solver_code/lib/qdldl/src/qdldl.c
solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/src/qdldl.c.o: solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/src/qdldl.c.o"
	cd /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build/solver_code/lib/qdldl && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/src/qdldl.c.o -MF CMakeFiles/qdldlstatic.dir/src/qdldl.c.o.d -o CMakeFiles/qdldlstatic.dir/src/qdldl.c.o -c /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/solver_code/lib/qdldl/src/qdldl.c

solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/src/qdldl.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/qdldlstatic.dir/src/qdldl.c.i"
	cd /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build/solver_code/lib/qdldl && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/solver_code/lib/qdldl/src/qdldl.c > CMakeFiles/qdldlstatic.dir/src/qdldl.c.i

solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/src/qdldl.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/qdldlstatic.dir/src/qdldl.c.s"
	cd /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build/solver_code/lib/qdldl && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/solver_code/lib/qdldl/src/qdldl.c -o CMakeFiles/qdldlstatic.dir/src/qdldl.c.s

# Object files for target qdldlstatic
qdldlstatic_OBJECTS = \
"CMakeFiles/qdldlstatic.dir/src/qdldl.c.o"

# External object files for target qdldlstatic
qdldlstatic_EXTERNAL_OBJECTS =

solver_code/lib/qdldl/out/libqdldl.a: solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/src/qdldl.c.o
solver_code/lib/qdldl/out/libqdldl.a: solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/build.make
solver_code/lib/qdldl/out/libqdldl.a: solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library out/libqdldl.a"
	cd /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build/solver_code/lib/qdldl && $(CMAKE_COMMAND) -P CMakeFiles/qdldlstatic.dir/cmake_clean_target.cmake
	cd /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build/solver_code/lib/qdldl && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/qdldlstatic.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/build: solver_code/lib/qdldl/out/libqdldl.a
.PHONY : solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/build

solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/clean:
	cd /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build/solver_code/lib/qdldl && $(CMAKE_COMMAND) -P CMakeFiles/qdldlstatic.dir/cmake_clean.cmake
.PHONY : solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/clean

solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/depend:
	cd /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/solver_code/lib/qdldl /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build/solver_code/lib/qdldl /Users/carlosm/Documents/guidance/scp_mpc/src/scp_socp_solver/c/build/solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : solver_code/lib/qdldl/CMakeFiles/qdldlstatic.dir/depend

