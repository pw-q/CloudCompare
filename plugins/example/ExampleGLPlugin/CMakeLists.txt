option( PLUGIN_EXAMPLE_GL "Install example GL plugin" ON )

if( PLUGIN_EXAMPLE_GL )
    project( ExampleGLPlugin )
    
    AddPlugin(
        NAME ${PROJECT_NAME}
        TYPE gl
    )
        
    add_subdirectory( include )
    add_subdirectory( src )
    add_subdirectory( shaders )
endif()
