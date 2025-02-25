#include "basic/Circle.h"
#include "Logger.h"

namespace pgl {

    Circle::Circle(glm::vec2 center, float radius)
        : m_Center(center), m_Radius(radius)
    {
        float vertices[] = {
            -0.5f, -0.5f,  0.0f,   // vertex 1
            -0.5f,  0.5f,  0.0f,   // vertex 2
             0.5f,  0.5f,  0.0f,   // vertex 3
             0.5f, -0.5f,  0.0f    // vertex 4
        };

        uint indices[] = {
            0, 1, 2,
            2, 3, 0
        };

        // vertex buffer 
        m_VBO = std::make_unique<VertexBuffer>(vertices, sizeof(vertices), GL_STATIC_DRAW);

        // vertex layout
        VertexBufferLayout layout;
        layout.Push<float>(3);

        // vertex array
        m_VAO = std::make_unique<VertexArray>();
        m_VAO->AddBuffer(*m_VBO, layout);


        // index array
        m_IBO = std::make_unique<IndexBuffer>(indices, 6);

        m_renderer = Renderer();
        
        m_MVP = m_Proj * m_View * m_Model;
        
        // Shaders

        const char* vert_source =   "#version 330 core\n"
                                    "layout (location = 0) in vec3 pos;\n"
                                    "out vec4 vertexColor;\n"
                                    "void main()\n"
                                    "{\n"
                                    "    gl_Position = vec4(pos, 1.0);\n"
                                    "    vertexColor = vec4(1, 0.0, 0.0, 1.0);\n"
                                    "}\0";

        const char* frag_source =
            "#version 330 core\n"
            "out vec4 FragColor;\n"
            "in vec4 vertexColor;\n"
            "void main()\n"
            "{\n"
            "    vec2 coord = 2. * gl_FragCoord.xy / vec2(500.0, 500.0) - vec2(1.0);\n"
            "    //if(gl_FragCoord.x * gl_FragCoord.x + gl_FragCoord.y * gl_FragCoord.y < 100000)\n"
            "    //if(gl_FragCoord.x < 470)\n"
            "    if (length(coord) < 0.5)\n"
            "    {\n"
            "        FragColor = vec4(1.0, 0.0, 0.0, 1.0);"
            "    }\n"
            "    else {\n"
            "        FragColor = vec4(0.0);\n"
            "    }\n"
            "}\0";


        m_Shader = std::make_unique<Shader>(vert_source, frag_source, false);
        m_Shader->Bind();
    }

    
    void Circle::OnRender() const
    {
        m_VAO->Bind();
        m_IBO->Bind();
        m_Shader->Bind();
        m_renderer.Draw(*m_VAO, *m_IBO, *m_Shader);
    }

    void Circle::OnImGuiRender()
    {
    }
};

