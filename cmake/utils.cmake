MACRO(INITIALISE_PROJECT)

#    SET(CMAKE_VERBOSE_MAKEFILE ON)
    SET(CMAKE_INCLUDE_CURRENT_DIR ON)

    # Required packages
    find_package(Qt5Widgets REQUIRED)
    # Keep track of some information about Qt
    SET(QT_BINARY_DIR ${_qt5Widgets_install_prefix}/bin)
    SET(QT_LIBRARY_DIR ${_qt5Widgets_install_prefix}/lib)
    SET(QT_PLUGINS_DIR ${_qt5Widgets_install_prefix}/plugins)
    SET(QT_VERSION_MAJOR ${Qt5Widgets_VERSION_MAJOR})
    SET(QT_VERSION_MINOR ${Qt5Widgets_VERSION_MINOR})
    SET(QT_VERSION_PATCH ${Qt5Widgets_VERSION_PATCH})
    
    # Some settings which depend on whether we want a debug or release version
    # of stVi
    IF(WIN32)
        STRING(REPLACE "/W3" "/W3 /WX" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        SET(LINK_FLAGS_PROPERTIES "/STACK:10000000 /MACHINE:X86")
    ENDIF()

    IF(CMAKE_BUILD_TYPE MATCHES [Dd][Ee][Bb][Uu][Gg])
        MESSAGE("Building a debug version...")
        set(DEBUG_MODE ON)
        # Default compiler settings
        IF(WIN32)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_DEBUG /MDd /Zi /Ob0 /Od /RTC1")
            SET(LINK_FLAGS_PROPERTIES "${LINK_FLAGS_PROPERTIES} /DEBUG")
        ELSE()
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
        ENDIF()
        # Make sure that debugging is on for Qt
        ADD_DEFINITIONS(-DQT_DEBUG)
    ELSE()
        MESSAGE("Building a release version...")
        set(DEBUG_MODE OFF)
        # Default compiler and linker settings
        IF(WIN32)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DNDEBUG /MD /O2 /Ob2")
        ELSE()
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ffast-math")
        ENDIF()
        # Make sure that debugging is off for Qt
        ADD_DEFINITIONS(-DQT_NO_DEBUG_OUTPUT)
        ADD_DEFINITIONS(-DQT_NO_DEBUG)
    ENDIF()

    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG -DDEBUG")
    SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG -DDEBUG")

    # Reduce the number of warnings
    # Remove "warning: multi-character character constant"
    OPTION(TREAT_WARNINGS_AS_ERRORS "Treat warnings as errors" ON)

    IF(WIN32)
        #TODO
    ELSE()
        IF(TREAT_WARNINGS_AS_ERRORS)
            SET(WARNING_ERROR "-Werror")
        ENDIF(TREAT_WARNINGS_AS_ERRORS)
        SET(DISABLED_WARNINGS "-Wno-multichar -Wno-unused-variable -Wno-unused-function -Wno-return-type -Wno-switch")
        SET(DISABLED_WARNINGS_DEBUG "-Wno-unused-variable -Wno-unused-function")
    ENDIF()

    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DISABLED_WARNINGS} ${WARNING_ERROR}")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DISABLED_WARNINGS_DEBUG}")

    # Ask for Unicode to be used
    ADD_DEFINITIONS(-DUNICODE)

    IF(WIN32)
        ADD_DEFINITIONS(-D_UNICODE)
    ENDIF()

    # Set the RPATH information on Linux
    # Note: this prevent us from having to use the uncool LD_LIBRARY_PATH...
    IF(NOT WIN32 AND NOT APPLE)
        SET(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../plugins/${PROJECT_NAME}")
    ENDIF()
    
    #enable c++11
    check_for_cxx11_compiler(CXX11_COMPILER)

    # If a C++11 compiler is available, then set the appropriate flags
    if(CXX11_COMPILER)
        enable_cxx11()
    else()
        message(FATAL_ERROR "Your compiler does not support c++11, UPDATE!!")
    endif()

    if(APPLE)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmacosx-version-min=10.7 -stdlib=libc++")
        SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -mmacosx-version-min=10.7 -stdlib=libc++")
    endif()

ENDMACRO()

MACRO(use_qt5lib qt5lib)
    find_package(${qt5lib} REQUIRED)
    include_directories(${${qt5lib}_INCLUDE_DIRS})
    add_definitions(${${qt5lib}_DEFINITIONS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${${qt5lib}_EXECUTABLE_COMPILE_FLAGS}")
ENDMACRO()

MACRO(LINUX_DEPLOY_QT_PLUGIN PLUGIN_CATEGORY)
    FOREACH(PLUGIN_NAME ${ARGN})
        # Deploy the Qt plugin itself
        INSTALL(FILES ${QT_PLUGINS_DIR}/${PLUGIN_CATEGORY}/${CMAKE_SHARED_LIBRARY_PREFIX}${PLUGIN_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}
                DESTINATION plugins/${PLUGIN_CATEGORY})
    ENDFOREACH()
ENDMACRO()

MACRO(PROJECT_GROUP TARGET_NAME FOLDER_PATH)
    #Organize projects into folders
    SET_PROPERTY(TARGET ${TARGET_NAME} PROPERTY FOLDER ${FOLDER_PATH})
ENDMACRO()

# Determines whether or not the compiler supports C++11
macro(check_for_cxx11_compiler _VAR)
    message(STATUS "Checking for C++11 compiler")
    set(${_VAR})
    if((MSVC AND MSVC10) OR
       (CMAKE_COMPILER_IS_GNUCXX AND NOT ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.6) OR
       (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND NOT ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 3.1))
        set(${_VAR} 1)
        message(STATUS "Checking for C++11 compiler - available")
    else()
        message(STATUS "Checking for C++11 compiler - unavailable")
    endif()
endmacro()

# Sets the appropriate flag to enable C++11 support
macro(enable_cxx11)
    if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++0x")
    endif()
endmacro()

