# pynest/CMakeLists.txt
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

if ( HAVE_PYTHON )

  if ( CYTHON_FOUND )
    include( UseCython )
    set_source_files_properties(
        pynestkernel.pyx
        PROPERTIES CYTHON_IS_CXX TRUE )

    cython_add_module( pynestkernel pynestkernel.pyx )
  else ()
    # require the source pynest/pynestkernel.cxx precompiled
    set( pynestkernel_generated_file "NOTFOUND" )
    if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/pynestkernel.cxx )
      set( pynestkernel_generated_file pynestkernel.cxx )
    elseif ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/pynestkernel.cpp )
      set( pynestkernel_generated_file pynestkernel.cpp )
    endif ()

    if ( pynestkernel_generated_file STREQUAL "NOTFOUND" )
      message( FATAL_ERROR "You have configured NEST to build the Python bindings PyNEST, but "
          "either there is no Cython version found on your system or you decided to "
          "not have NEST cythonize the bindings. Either install Cython (>=0.19.2) or cythonize the "
          "PyNEST bindings yourself on a system with a suitable version of Cython (>=0.19.2):\n"
          "    cd <NEST sources>/pynest\n"
          "    cythonize pynestkernel.pyx\n"
          "And copy the file `pynestkernel.cpp` into ${CMAKE_CURRENT_SOURCE_DIR}/pynest ." )
    endif ()
    python_add_module( pynestkernel ${pynestkernel_generated_file} pynestkernel.pxd )
    if ( APPLE )
      set_target_properties( pynestkernel PROPERTIES LINK_FLAGS "-undefined dynamic_lookup" )
    else ()
      set_target_properties( pynestkernel PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")
      target_link_libraries( pynestkernel ${Python_LIBRARIES} )
    endif ()
  endif ()

  target_link_libraries( pynestkernel
      nest_lib nestutil nestkernel sli_lib models ${EXTERNAL_MODULE_LIBRARIES}
      )

  target_include_directories( pynestkernel PRIVATE
      ${PROJECT_BINARY_DIR}/libnestutil
      ${PROJECT_SOURCE_DIR}/libnestutil
      ${PROJECT_SOURCE_DIR}/nest
      ${PROJECT_SOURCE_DIR}/nestkernel
      ${PROJECT_SOURCE_DIR}/nestkernel/spatial
      ${PROJECT_SOURCE_DIR}/sli
      ${PROJECT_SOURCE_DIR}/thirdparty
      ${PROJECT_SOURCE_DIR}/models
      ${Python_INCLUDE_DIRS}
      )

  target_compile_definitions( pynestkernel PRIVATE
      -D_IS_PYNEST
      )

  install(DIRECTORY  nest/ ${PROJECT_BINARY_DIR}/pynest/nest/
      DESTINATION ${CMAKE_INSTALL_PREFIX}/${PYEXECDIR}/nest
      PATTERN "versionchecker.py.in" EXCLUDE
  )
  install( TARGETS pynestkernel DESTINATION ${PYEXECDIR}/nest/ )

  install( DIRECTORY examples/
      DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples/pynest
  )

endif ()
