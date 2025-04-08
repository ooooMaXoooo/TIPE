#include "basic/Circle.h"
#include "Logger.h"

#include <cmath>


namespace pgl {

    Circle::Circle(glm::vec3 center, float radius, GLFWwindow* parent_window)
        : m_Center(center), m_Radius(radius), Drawable(parent_window)
    {
        Init();


        m_renderer = Renderer();


        UpdateTransform();
        m_Shader->SetUniformMat4f("iMVP", m_MVP);

           

        m_Rotations[0] = 55.0f;
        m_Rotations[1] = 0.0f;
        m_Rotations[2] = 0.0f;

        /*
        m_Translations[0] = 0.0f;
        m_Translations[1] = 0.0f;
        m_Translations[2] = -10.0f;
        */

        m_Scale[0] = 1.0f;
        m_Scale[1] = 1.0f;
        m_Scale[2] = 1.0f;
    }




    void Circle::Init()
    {
        // find all points on the circle
        constexpr uint segment_count = 100;

        std::vector<float> vertices;
        std::vector<uint> indices;

        // we reserve the place for our vertices
        // 3 coordinates for each point on the circle
        // +3 for the center
        vertices.reserve(3 * segment_count + 3);

        // we add the center first
        vertices.emplace_back(m_Center.x);
        vertices.emplace_back(m_Center.y);
        vertices.emplace_back(0.0f);

        // we fill the vertices with trigonometry
        for (uint i = 0; i < segment_count; i++)
        {
            float theta = 2 * PI * i / segment_count;

            float x = m_Radius * cos(theta);
            float y = m_Radius * sin(theta);

            vertices.emplace_back(x);
            vertices.emplace_back(y);
            vertices.emplace_back(0.0f);
        }


        // we reserve the place for the indices
        indices.reserve(3 * segment_count);

        // we set indices vector exept the last one
        for (uint i = 1; i < segment_count; i++)
        {
            indices.emplace_back(0);
            indices.emplace_back(i);
            indices.emplace_back(i + 1);
        }

        // we put the last one in :

        indices.emplace_back(0);
        indices.emplace_back(segment_count);
        indices.emplace_back(1);




        /*
        float vertices[] = {
            -1.0f, -1.0f,  0.0f,   // vertex 1
            -1.0f,  1.0f,  0.0f,   // vertex 2
             1.0f,  1.0f,  0.0f,   // vertex 3
             1.0f, -1.0f,  0.0f    // vertex 4
        };

        uint indices[] = {
            0, 1, 2,
            2, 3, 0
        };*/

        // vertex buffer
        m_VBO.reset();
        m_VBO = std::make_unique<VertexBuffer>(vertices.data(), vertices.size() * sizeof(float), GL_STATIC_DRAW);

        // vertex layout
        VertexBufferLayout layout;
        layout.Push<float>(3); // x, y, z

        // vertex array
        m_VAO.reset();
        m_VAO = std::make_unique<VertexArray>();
        m_VAO->AddBuffer(*m_VBO, layout);


        // index array
        m_IBO.reset();
        m_IBO = std::make_unique<IndexBuffer>(indices.data(), indices.size());



        // Shaders

        const char* vert_source =
            "#version 330 core\n"
            "layout (location = 0) in vec3 pos;\n"
            "out vec4 vertexColor;\n"
            "uniform mat4 iMVP;\n"
            "uniform float iTime;\n"
            "void main()\n"
            "{\n"
            "    gl_Position = iMVP * vec4(pos, 1.0);\n"
            "    vertexColor = gl_Position;\n"
            "}\0";

        const char* frag_source =
            "#version 330 core\n"
            "out vec4 FragColor;\n"
            "in vec4 vertexColor;\n"
            "uniform float iTime;\n"
            "uniform mat4 iMVP;\n"
            "void main()\n"
            "{\n"
            "    /*vec3 coord = (2. * (vertexColor).xyz / vec3(500.0, 500.0, 5000000.0f)) - vec3(1.0);\n"
            "    if (length(coord.xyz) < 2.0f)\n"
            "    {\n"
            "        FragColor = vec4(2. * gl_FragCoord.xyz / vec3(500.0, 500.0, 500.0f) - vec3(1.0), 1.0f);//vec4(vertexColor.xyz, 1.0f);\n"
            "    }\n"
            "    else {\n"
            "        FragColor = vec4(0.0);\n"
            "    }*/\n"
            "    FragColor = vec4(vertexColor.xyz, 1.0f) * abs(sin(iTime));\n"
            "}\0";

        m_Shader.reset();
        m_Shader = std::make_unique<Shader>(vert_source, frag_source, false);
        m_Shader->Bind();
    }


    void Circle::UpdateTransform()
    {
        // handle transform
        m_Model = glm::mat4(1.0f);

        // Apply the translation to take the center into account
        m_Model = glm::translate(m_Model, m_Center);

        m_Model = glm::rotate(m_Model, glm::radians(m_Rotations[0]), glm::vec3(1.0f, 0.0f, 0.0f));
        m_Model = glm::rotate(m_Model, glm::radians(m_Rotations[1]), glm::vec3(0.0f, 1.0f, 0.0f));
        m_Model = glm::rotate(m_Model, glm::radians(m_Rotations[2]), glm::vec3(0.0f, 0.0f, 1.0f));


        // Scale the object using the radius
        m_Model = glm::scale(m_Model, glm::vec3(m_Radius));


        // handle scale
        m_Model = glm::scale(m_Model, glm::vec3(m_Scale[0], m_Scale[1], m_Scale[2]));


        // handle translations
        m_View = glm::mat4(1.0f);
        //m_View = glm::translate(m_View, glm::vec3(m_Translations[0], m_Translations[1], m_Translations[2]));



        m_Proj = glm::perspective(glm::radians(45.0f), 500.0f / 500.0f, 0.1f, 100.0f);

        m_MVP = m_Proj * m_View * m_Model;
    }



    void Circle::SetRadius(float radius)
    {
        m_Radius = radius;
    }

    void Circle::SetPosition(glm::vec3 center)
    {
        m_Center = center;
    }

    void Circle::SetScale(glm::vec3 scale)
    {
        m_Scale[0] = scale[0];
        m_Scale[1] = scale[1];
        m_Scale[2] = scale[2];
    }

    void Circle::SetRotation(glm::vec3 rotation)
    {
        m_Rotations[0] = rotation.x;
        m_Rotations[1] = rotation.y;
        m_Rotations[2] = rotation.z;
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

        ImGui::SliderFloat3("position (x,y,z)", &m_Center[0], -20, 20);
        ImGui::SliderFloat3("rotation (x,y,z)", m_Rotations, 0, 360);
        ImGui::SliderFloat3("scaling (x,y,z)", m_Scale, 0, 5);
        ImGui::Separator();
        ImGui::Text("In the shader : length(coord) < %.3f", 2.0f * abs(glm::sin(glfwGetTime())));

        /*// we change the position only if it change
        if (ImGui::IsItemDeactivatedAfterEdit()) {
            SetPosition(m_Center);
        }*/

        UpdateTransform();


        ImGui::End();
    }
};

