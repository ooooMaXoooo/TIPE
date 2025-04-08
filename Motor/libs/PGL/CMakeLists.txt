cmake_minimum_required(VERSION 3.10)

# Définir le projet
project(PGL VERSION 1.1)

# Réglage de la politique CMP0002 pour éviter les conflits de cibles
if(POLICY CMP0002)
    cmake_policy(SET CMP0002 NEW)
endif()

# Réglage de la politique CMP0079 pour autoriser la liaison de bibliothèques non construites dans le même répertoire
if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)
endif()

# Définir la version C++
set(CMAKE_CXX_STANDARD 20)

# Définir le répertoire de sortie des exécutables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Définir les fichiers sources pour ImGui
set(ImGUI_SOURCES
    "libs/ImGUI/imgui.cpp"
    "libs/ImGUI/imgui_demo.cpp"
    "libs/ImGUI/imgui_draw.cpp"
    "libs/ImGUI/imgui_impl_glfw.cpp"
    "libs/ImGUI/imgui_impl_opengl3.cpp"
    "libs/ImGUI/imgui_stdlib.cpp"
    "libs/ImGUI/imgui_tables.cpp"
    "libs/ImGUI/imgui_widgets.cpp"
)

# Définir les fichiers sources pour gérer OpenGL
set(OPENGL_SOURCES
    "src/openGL/IndexBuffer.cpp"
    "src/openGL/Renderer.cpp"
    "src/openGL/Shader.cpp"
    "src/openGL/VertexArray.cpp"
    "src/openGL/VertexBuffer.cpp"
    "src/openGL/VertexBufferLayout.cpp"
)

set (GLAD_SOURCES
    "libs/glad/glad.c"
)

set (BASIC_SOURCES
    "src/basic/Circle.cpp"
)

set (CORE_SOURCE
    "src/Core/Window.cpp"
    #"src/Core/ObjectsHandler.cpp"
)

set (PGL_SOURCES
    "src/Logger.cpp"
    #"src/pgl.cpp"
    "src/Drawable.cpp"
    "src/Camera.cpp"
    "${BASIC_SOURCES}"
    "${CORE_SOURCE}"
)

# Définir les fichiers sources principaux
set(SOURCES
    "${PGL_SOURCES}"
    "${ImGUI_SOURCES}"
    "${OPENGL_SOURCES}"
    "${GLAD_SOURCES}"
)


message(STATUS "Configurations disponibles : ${CMAKE_CONFIGURATION_TYPES}")
message(STATUS "Configuration sélectionnée : ${CMAKE_BUILD_TYPE}")

# Ajout des sous-répertoires pour GLFW et GLM
add_subdirectory("libs/glfw")
add_subdirectory("libs/glm")

# Création de la bibliothèque statique
add_library(PGL STATIC "${SOURCES}")


# Ajouter les répertoires d'inclusion
target_include_directories(PGL PUBLIC
    "include/"
    
    "include/ImGui"

    "include/pgl/openGL"
    "include/pgl"

    "libs/glfw/include"

    "libs/glm/"
)

# Lier les bibliothèques
target_link_libraries(PGL PRIVATE
    glfw
    glm
    OpenGL32  # OpenGL standard sur Windows
)


install(TARGETS PGL
    CONFIGURATIONS Debug Release
    ARCHIVE DESTINATION lib/$<CONFIG>
    LIBRARY DESTINATION lib/$<CONFIG>
    RUNTIME DESTINATION bin/$<CONFIG>
)

install(DIRECTORY include/ DESTINATION include)