# nest/CMakeLists.txt
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

set( nest_sources
    neststartup.cpp neststartup.h
    )

add_executable( nest main.cpp ${nest_sources} )

if ( NOT APPLE )
    set_target_properties( nest
        PROPERTIES
        LINK_FLAGS "-Wl,--no-as-needed"
        )
endif ()

add_library( nest_lib ${nest_sources} )

if ( NOT APPLE )
    set_target_properties( nest_lib
        PROPERTIES
        OUTPUT_NAME nest
        LINK_FLAGS "-Wl,--no-as-needed"
        )
else ()
    set_target_properties( nest_lib
        PROPERTIES
        OUTPUT_NAME nest

        # delay lookup of symbols from libpython when building with MPI4Py
        LINK_FLAGS "-Wl,-undefined -Wl,dynamic_lookup"
        )
endif ()

set_target_properties( nest_lib
    PROPERTIES
    VERSION ${NEST_VERSION}
    SOVERSION 3
    )

target_link_libraries( nest
    nestutil nestkernel sli_lib sli_readline models
    ${EXTERNAL_MODULE_LIBRARIES} OpenMP::OpenMP_CXX )

target_link_libraries( nest_lib
    nestutil nestkernel sli_lib models ${EXTERNAL_MODULE_LIBRARIES}
    OpenMP::OpenMP_CXX )

target_include_directories( nest PRIVATE
    ${PROJECT_SOURCE_DIR}/thirdparty
    ${PROJECT_BINARY_DIR}/nest
    ${PROJECT_BINARY_DIR}/libnestutil
    ${PROJECT_SOURCE_DIR}/libnestutil
    ${PROJECT_SOURCE_DIR}/sli
    ${PROJECT_SOURCE_DIR}/nestkernel
    ${PROJECT_SOURCE_DIR}/models
    )


target_compile_definitions( nest PRIVATE
    -D_BUILD_NEST_CLI
)

target_include_directories( nest_lib PRIVATE
    ${PROJECT_SOURCE_DIR}/thirdparty
    ${PROJECT_BINARY_DIR}/nest
    ${PROJECT_SOURCE_DIR}/libnestutil
    ${PROJECT_BINARY_DIR}/libnestutil
    ${PROJECT_SOURCE_DIR}/sli
    ${PROJECT_SOURCE_DIR}/nestkernel
    ${PROJECT_SOURCE_DIR}/models
    ${Python_INCLUDE_DIRS}
    )


if ( HAVE_PYTHON )
  target_compile_definitions( nest_lib PRIVATE
      -D_IS_PYNEST
      )
endif ()

install( TARGETS nest nest_lib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/nest
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/nest
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

FILTER_HEADERS("${nest_sources}" install_headers )
install( FILES ${install_headers}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nest)
