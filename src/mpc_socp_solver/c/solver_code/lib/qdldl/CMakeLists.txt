# Minimum version required
cmake_minimum_required (VERSION 3.13)

set(QDLDL_VERSION_MAJOR "0")
set(QDLDL_VERSION_MINOR "1")
set(QDLDL_VERSION_PATCH "7")
set(QDLDL_VERSION "${QDLDL_VERSION_MAJOR}.${QDLDL_VERSION_MINOR}.${QDLDL_VERSION_PATCH}")

# Project name
project(qdldl VERSION ${QDLDL_VERSION})

include( CMakeDependentOption )

option( QDLDL_BUILD_STATIC_LIB "Build the static library" ON )
option( QDLDL_BUILD_SHARED_LIB "Build the shared library" ON )

cmake_dependent_option( QDLDL_BUILD_DEMO_EXE
                        "Build the demo executable (requires the static library)"
                        ON    # Default to on
                        QDLDL_BUILD_STATIC_LIB OFF ) # Force off if the static library isn't built

cmake_dependent_option( QDLDL_UNITTESTS
                        "Build the unit testing suite"
                        OFF    # Default to off
                        QDLDL_BUILD_STATIC_LIB OFF ) # Force off if the static library isn't built

# Set the output folder where your program will be created
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/out)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/out)

# Some non-standard CMake modules
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/configure/cmake)

# Export compilation commands for IDEs and autocompletion
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
# ----------------------------------------------
# Use floats instead of doubles
option( QDLDL_FLOAT "Use float numbers instead of doubles" OFF )
if( ${QDLDL_FLOAT} )
    message( STATUS "Using single precision floats" )
else()
    message( STATUS "Using double precision floats" )
endif()

# Use long integers for indexing
option( QDLDL_LONG "Use long integers (64bit) for indexing" ON )

if( NOT (CMAKE_SIZEOF_VOID_P EQUAL 8) )
	message(STATUS "Disabling long integers (64bit) on 32bit machine")
	set(QDLDL_LONG OFF)
endif()
message(STATUS "Long integers (64bit) are ${QDLDL_LONG}")


# Set Compiler flags
# ----------------------------------------------
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # -fPIC


# Add compiler options if we are not on windows
if (NOT MSVC)

    if (COVERAGE)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")

	      if(FORTRAN)
		        set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} --coverage")
	      endif(FORTRAN)
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g")
endif (NOT MSVC)

# Generate header file with the global options
# ---------------------------------------------

# numeric types
if(QDLDL_FLOAT)
  set(QDLDL_FLOAT_TYPE "float")
  set(QDLDL_FLOAT 1)
else()
	set(QDLDL_FLOAT_TYPE "double")
endif()

if(QDLDL_LONG)
  set(QDLDL_INT_TYPE "long long")
	set(QDLDL_INT_TYPE_MAX "LLONG_MAX")
  set(QDLDL_LONG 1)
else()
	set(QDLDL_INT_TYPE "int")
	set(QDLDL_INT_TYPE_MAX "INT_MAX")
endif()

#boolean type is always unsigned char
#for now, since _Bool does not exist in
#C89 and we want to avoid interoperability
#problems when calling QDLDL from C++
set(QDLDL_BOOL_TYPE "unsigned char")

# Generate header file with the global options
# ---------------------------------------------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/configure/qdldl_types.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/include/qdldl_types.h
               NEWLINE_STYLE LF)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/configure/qdldl_version.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/include/qdldl_version.h
               NEWLINE_STYLE LF)


# Set sources
# ----------------------------------------------
set(
	qdldl_src
	src/qdldl.c
	)

set(
	qdldl_headers
	include/qdldl.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/qdldl_types.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/qdldl_version.h
	)

# Create object library
# ----------------------------------------------
add_library (qdldlobject OBJECT ${qdldl_src} ${qdldl_headers})
target_include_directories(qdldlobject PRIVATE
                           ${CMAKE_CURRENT_SOURCE_DIR}/include)


# Create Static Library
# ----------------------------------------------

include(GNUInstallDirs)

message( STATUS "Static library build is ${QDLDL_BUILD_STATIC_LIB}" )

if( QDLDL_BUILD_STATIC_LIB )
    # Static library
    add_library (qdldlstatic STATIC ${qdldl_src} ${qdldl_headers})
    # Give same name to static library output
    set_target_properties(qdldlstatic PROPERTIES OUTPUT_NAME qdldl)

    # Declare include directories for the cmake exported target
    target_include_directories(qdldlstatic
                               PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>"
                                      "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/qdldl>")

    # Install Static Library
    # ----------------------------------------------

    install(TARGETS qdldlstatic
            EXPORT  ${PROJECT_NAME}
            ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif()

# Install Headers
# ----------------------------------------------
install(FILES ${qdldl_headers} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/qdldl")


# Install Shared Library
# ----------------------------------------------
# Create qdldl shared library
message( STATUS "Shared library build is ${QDLDL_BUILD_SHARED_LIB}" )

if( QDLDL_BUILD_SHARED_LIB )
    add_library (qdldl SHARED ${qdldl_src} ${qdldl_headers})

    # Declare that we are building the shared library to get proper symbol exports.
    # Shared library consumers should also define QDLDL_SHARED_LIB to get the library
    # exports properly, so we do it for them in the CMake interface by defining it as
    # a PUBLIC compile definition.
    target_compile_definitions(qdldl PRIVATE BUILDING_QDLDL)
    target_compile_definitions(qdldl PUBLIC  QDLDL_SHARED_LIB)

    # Declare include directories for the cmake exported target
    target_include_directories(qdldl
        PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
               "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>"
               "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/qdldl>")

    # Install qdldl shared library
    install(TARGETS qdldl
        EXPORT  ${PROJECT_NAME}
        LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif()

# Create CMake packages for the build directory
# ----------------------------------------------
if( QDLDL_BUILD_SHARED_LIB OR QDLDL_BUILD_STATIC_LIB)
    include(CMakePackageConfigHelpers)

    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/qdldl-config-version.cmake"
      VERSION ${QDLDL_VERSION}
      COMPATIBILITY SameMajorVersion
    )

    export(EXPORT ${PROJECT_NAME}
      FILE "${CMAKE_CURRENT_BINARY_DIR}/qdldl-targets.cmake"
      NAMESPACE qdldl::)

    if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/qdldl-config.cmake)
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/qdldl-config.cmake "include(\"\${CMAKE_CURRENT_LIST_DIR}/qdldl-targets.cmake\")\n")
    endif()


    # Create CMake packages for the install directory
    # ----------------------------------------------

    set(ConfigPackageLocation ${CMAKE_INSTALL_LIBDIR}/cmake/qdldl)

    install(EXPORT ${PROJECT_NAME}
            FILE qdldl-targets.cmake
            NAMESPACE qdldl::
            DESTINATION ${ConfigPackageLocation})

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/qdldl-config.cmake
                  ${CMAKE_CURRENT_BINARY_DIR}/qdldl-config-version.cmake
            DESTINATION ${ConfigPackageLocation})

    # Add uninstall command
    # ----------------------------------------------
    if(NOT TARGET uninstall)
        configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/configure/cmake/cmake_uninstall.cmake.in"
            "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
            IMMEDIATE @ONLY)

        add_custom_target(uninstall
            COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
    endif()
endif()


# Add testing
# ----------------------------------------------
# Add custom command to generate tests
message( STATUS "Unit testing suite build is ${QDLDL_UNITTESTS}" )

if( QDLDL_UNITTESTS )
    # Add test_headers and codegen_test_headers
    add_subdirectory(tests)

    # Direct qdldl solver testing
    add_executable(qdldl_tester
                ${PROJECT_SOURCE_DIR}/tests/qdldl_tester.c ${PROJECT_SOURCE_DIR}/tests/minunit.h
                ${test_headers})
    target_link_libraries (qdldl_tester qdldlstatic)

    # Add testing
    include(CTest)
    enable_testing()
    add_test(NAME tester COMMAND $<TARGET_FILE:qdldl_tester>)
endif()
