set(Geometry_Sources
    GeoGlobal.h
    Geometry2D.cpp
    Geometry2D.h
    Geometry3D.h
    List.h
    MathDefs.h
#    Operation3D.cpp
    Operation3D.h
    Triangulate.cpp
    Triangulate.h
    Vector2D.h
    Vector2D_float.h
    Vector3D.h
    Vector3D_float.h)

#set(Geometry_Required_Libs
#   debug MSVCRTD.LIB
#   debug MSVCPRTD.LIB
#   optimized MSVCRT.LIB
#   optimized MSVCPRT.LIB)


add_library(Geometry ${Geometry_Sources})
# libraries that needs to be linked to the present library
#target_link_libraries(Geometry ${Geometry_Required_Libs})



