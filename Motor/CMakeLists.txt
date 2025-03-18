# CMakeList.txt : CMake project for Motor, include source and define
# project specific logic here.
#

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Réglage de la politique CMP0002 pour éviter les conflits de cibles
if(POLICY CMP0002)
    cmake_policy(SET CMP0002 NEW)
endif()

# Réglage de la politique CMP0079 pour autoriser la liaison de bibliothèques non construites dans le même répertoire
if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)
endif()


set(SOURCE_FILES
    "src/pch.cpp"
    "src/Motor.cpp"

    "src/General/Vector.cpp"

    "src/Motor/2D/MassiveObject2D.cpp"
    "src/Motor/2D/Planet2D.cpp"
    "src/Motor/Entity.cpp"

    "src/Motor/Core/ObjectsHandler.cpp"
    "src/Motor/Core/MassiveObject.cpp"
)

# Create executable
add_executable(Motor ${SOURCE_FILES})

# Include directories
target_include_directories(Motor PUBLIC
    "include"
    "${CMAKE_CURRENT_SOURCE_DIR}/libs/PGL/include"  # Inclure les headers de PGL
    "${CMAKE_CURRENT_SOURCE_DIR}/libs/PGL/include/pgl"
)



# Ajouter PGL comme sous-projet
add_subdirectory(libs/PGL)


# preprocessor macro
target_compile_definitions(PGL PRIVATE 
    $<$<CONFIG:Debug>:Mt_Debug>
    $<$<CONFIG:Release>:Mt_Release>
)



# Utiliser un header précompilé
target_precompile_headers (Motor PRIVATE "include/pch.h")

# Lier la bibliothèque PGL à l'exécutable Motor
target_link_libraries(Motor PRIVATE PGL)


if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Motor PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
