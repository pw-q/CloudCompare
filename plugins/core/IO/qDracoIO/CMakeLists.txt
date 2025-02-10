option( PLUGIN_IO_QDRACO "Install qDracoIO plugin to read & write DRC files" ON )

if( PLUGIN_IO_QDRACO )
	project( QDRACO_IO_PLUGIN )

	find_package( Draco REQUIRED )

	AddPlugin( NAME ${PROJECT_NAME} TYPE io )



	if( MSVC )

		target_include_directories( ${PROJECT_NAME} PRIVATE ${DRACO_INCLUDE_DIRS} )
		target_link_libraries( ${PROJECT_NAME} ${Draco_DIR}/../../lib/${DRACO_LIBRARIES}.lib )

	elseif( UNIX )

		target_include_directories( ${PROJECT_NAME} PRIVATE ${DRACO_INCLUDE_DIR} )
		set ( DRACO_INCLUDE_DIR "/usr/local/include/draco" CACHE PATH "Draco include dir" )
		set ( DRACO_LIB_DIR "/usr/local/lib" CACHE PATH "Draco library dir" )
		message("include:" ${DRACO_INCLUDE_DIR})
		
		if( NOT DRACO_INCLUDE_DIR )
			message( FATAL_ERROR "Draco include dir not specified (DRACO_INCLUDE_DIR)" )
		endif()
		if( NOT DRACO_LIB_DIR )
			message( FATAL_ERROR "Draco library dir not specified (DRACO_LIB_DIR)" )
		endif()
		
		target_link_libraries( ${PROJECT_NAME} ${DRACO_LIBRARIES} )

	elseif()

		target_link_libraries( ${PROJECT_NAME} ${DRACO_LIBRARIES} )

	endif()

	add_subdirectory( include )
	add_subdirectory( src )
	add_subdirectory( ui )
 
endif()
