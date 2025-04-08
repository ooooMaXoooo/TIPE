#include "Motor/3D/Planet3D.h"

namespace Motor {

    Planet3D::Planet3D(ldouble mass, glm::vec3 center, double radius, GLFWwindow* parent_window, Renderer& renderer, glm::vec3 v0, glm::vec3 a0)
        : MassiveObject3D{ mass, center, parent_window, v0, a0 }, m_Radius(radius)
    {
        m_renderer = renderer;
        InitPlanet();
    }

    Planet3D::Planet3D(const Planet3D& p)
        : m_Radius(p.m_Radius), MassiveObject3D{ p.m_Mass, p.m_Pos, p.m_Parent_window, p.m_Velocity, p.m_Acceleration }
    {
        m_renderer = p.m_renderer;
        InitPlanet();
    }

    Planet3D::~Planet3D()
    {
        // need to handle everything
    }

    void Planet3D::Update(float ts)
    {
        // look for forces  --> job done by the objectHandler
        // apply the resultant to the acceleration
        // ( WARNING : do not ADD the forces to the acceleration / velocity, Newton's laws is   acceleration IS EQUAL TO the sum of forces)
        m_Acceleration = m_Forces / glm::vec3(m_Mass);

        // check for potential collisions

        // if there is no collisions, let the physics make its job
        // 
        // if there is collisions, delete both planets
        //  ---->  the desctructor need to clear everything potentially including rendering stuff
        //  ----> it's handle by the objectHandler

        // apply half the acceleration to velocity
        const glm::vec3 half_acc = m_Acceleration * ts * 0.5f;
        m_Velocity += half_acc;

        // apply velocity to position
        m_Pos += m_Velocity * ts;

        // apply the remaining half of acceleration to velocity in order to correct the mean
        m_Velocity += half_acc;


        // update the view
        UpdateTransform();

        // done ?
    }


    void Planet3D::OnRender()
    {
        m_VAO->Bind();
        m_IBO->Bind();
        m_Shader->Bind();


        m_MVP = m_Proj * m_View * m_Model;

        m_Shader->SetUniformMat4f("iMVP", m_MVP);
        m_Shader->SetUniform1f("iTime", glfwGetTime());

        m_renderer.Draw(*m_VAO, *m_IBO, *m_Shader);
    }

    void Planet3D::OnImGuiRender()
    {
    }





    void Planet3D::InitPlanet()
    {
        float vertices[] =
        {   -0.5f, -0.5f, 0.0f,
            -0.5f,  0.5f, 0.0f,
             0.5f,  0.5f, 0.0f,
             0.5f, -0.5f, 0.0f
        };

        uint indices[] =
        {   0, 1, 2,
            2, 3, 0
        };

        m_VBO = std::make_unique<VertexBuffer>(vertices, sizeof(vertices));

        // vertex buffer layout
        VertexBufferLayout layout;
        layout.Push<float>(3); // x, y, z

        // vertex array
        m_VAO = std::make_unique<VertexArray>();
        m_VAO->AddBuffer(*m_VBO, layout);


        // index buffer
        m_IBO = std::make_unique<IndexBuffer>(indices, 6);


        /// Shaders
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
            "//uniform mat4 iMVP;\n"
            "uniform float iRadius;\n"
            "void main()\n"
            "{\n"
            "    vec3 coord = (2. * (vertexColor).xyz / vec3(500.0, 500.0, 1.0f)) - vec3(1.0);\n"
            "    if (length(vertexColor.xy) < iRadius)\n"
            "    {\n"
            "        FragColor = vec4(vertexColor.xyz, 1.0f);\n"
            "    }\n"
            "    else {\n"
            "        FragColor = vec4(0.0);\n"
            "    }\n"
            "    //FragColor = vec4(vertexColor.xyz, 1.0f);\n"
            "}\0";

        m_Shader = std::make_unique<Shader>(vert_source, frag_source, false);
        m_Shader->Bind();

        m_Shader->SetUniform1f("iRadius", m_Radius);

        UpdateTransform();
        m_Shader->SetUniformMat4f("iMVP", m_MVP);
    }

    void Planet3D::UpdateTransform()
    {
        // handle transform
        m_Model = glm::mat4(1.0f);

        // Apply the translation to take the center into account
        m_Model = glm::translate(m_Model, m_Pos);

        /*m_Model = glm::rotate(m_Model, glm::radians(m_Rotations[0]), glm::vec3(1.0f, 0.0f, 0.0f));
        m_Model = glm::rotate(m_Model, glm::radians(m_Rotations[1]), glm::vec3(0.0f, 1.0f, 0.0f));
        m_Model = glm::rotate(m_Model, glm::radians(m_Rotations[2]), glm::vec3(0.0f, 0.0f, 1.0f));*/


        // Scale the object using the radius
        m_Model = glm::scale(m_Model, glm::vec3(m_Radius));


        // handle scale
        //m_Model = glm::scale(m_Model, glm::vec3(m_Scale[0], m_Scale[1], m_Scale[2]));


        // handle translations
        m_View = glm::mat4(1.0f);
        //m_View = glm::translate(m_View, glm::vec3(m_Translations[0], m_Translations[1], m_Translations[2]));

        int width = 0, height = 0;
        glfwGetWindowSize(m_Parent_window, &width, &height);
        m_Proj = glm::perspective(glm::radians(45.0f), (float)width / (float)height, 0.1f, 100.0f);

        m_MVP = m_Proj * m_View * m_Model;
    }

};