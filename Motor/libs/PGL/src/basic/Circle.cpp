#include "basic/Circle.h"
#include "Logger.h"

namespace pgl {

    Circle::Circle(glm::vec2 center, float radius)
        : m_Center(center), m_Radius(radius)
    {
        float vertices[] = {
            -1.0f, -1.0f,  0.0f,   // vertex 1
            -1.0f,  1.0f,  0.0f,   // vertex 2
             1.0f,  1.0f,  0.0f,   // vertex 3
             1.0f, -1.0f,  0.0f    // vertex 4
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
        
        
        // Shaders

        const char* vert_source =   
            "#version 330 core\n"
            "layout (location = 0) in vec3 pos;\n"
            "out vec4 vertexColor;\n"
            "out vec3 Epos;\n"
            "uniform mat4 iMVP;\n"
            "uniform float iTime;\n"
            "void main()\n"
            "{\n"
            "    gl_Position = iMVP * vec4(pos, 1.0);\n"
            "    vertexColor = gl_Position;\n"
            "    Epos = gl_Position.xyz;"
            "}\0";

        const char* frag_source =
            "#version 330 core\n"
            "out vec4 FragColor;\n"
            "in vec4 vertexColor;\n"
            "in vec3 Epos;\n"
            "uniform float iTime;\n"
            "void main()\n"
            "{\n"
            "    vec3 coord = 2. * Epos / vec3(500.0, 500.0, 10.0f) - vec3(1.0);\n"
            "    //if(gl_FragCoord.x * gl_FragCoord.x + gl_FragCoord.y * gl_FragCoord.y < 100000)\n"
            "    //if(gl_FragCoord.x < 470)\n"
            "    if (length(coord) < 1.5)\n"
            "    {\n"
            "        FragColor = vec4(vertexColor.xyz, 1.0f) * abs(sin(iTime));"
            "    }\n"
            "    else {\n"
            "        FragColor = vec4(0.0);\n"
            "    }\n"
            "    //FragColor = vec4(vertexColor.xyz, 1.0f) * abs(sin(iTime));\n"
            "}\0";


        m_Shader = std::make_unique<Shader>(vert_source, frag_source, false);
        m_Shader->Bind();


        m_Model = glm::mat4(1.0f);
        m_Model = glm::rotate(m_Model, glm::radians(-55.0f), glm::vec3(1.0f, 0.0f, 0.0f));

        m_View = glm::mat4(1.0f);
        // note that we're translating the scene in the reverse direction of where we want to move
        m_View = glm::translate(m_View, glm::vec3(0.0f, 0.0f, -6.0f));

        m_Proj = glm::perspective(glm::radians(45.0f), 500.0f / 500.0f, 0.1f, 100.0f);

        m_MVP = m_Proj * m_View * m_Model;
        m_Shader->SetUniformMat4f("iMVP", m_MVP);

        m_Rotations[0] = -55.0f;
        m_Rotations[1] = 0;
        m_Rotations[2] = 0;
    }

    
    void Circle::OnRender()
    {
        m_VAO->Bind();
        m_IBO->Bind();
        m_Shader->Bind();


        m_MVP = m_Proj * m_View * m_Model;

        m_Shader->SetUniformMat4f("iMVP", m_MVP);
        m_Shader->SetUniform1f("iTime", glfwGetTime());

        m_renderer.Draw(*m_VAO, *m_IBO, *m_Shader);
    }

    void Circle::OnImGuiRender()
    {
        ImGui::Begin("Circle");

        float temp[3];
        for (int i = 0; i < 3; i++) {
            temp[i] = m_Rotations[i];
        }
        
        ImGui::SliderFloat3("rotation (x,y,z)", m_Rotations, 0, 360);

        if (m_Rotations[0] != temp[0])
        {
            m_Model = glm::rotate(m_Model, glm::radians(glm::radians(m_Rotations[0])), glm::vec3(1.0f, 0.0f, 0.0f));
        }
        if (m_Rotations[1] != temp[1])
        {
            m_Model = glm::rotate(m_Model, glm::radians(glm::radians(m_Rotations[1])), glm::vec3(0.0f, 1.0f, 0.0f));
        }
        if (m_Rotations[2] != temp[2])
        {
            m_Model = glm::rotate(m_Model, glm::radians(glm::radians(m_Rotations[2])), glm::vec3(0.0f, 0.0f, 1.0f));
        }

        ImGui::End();
    }
};

