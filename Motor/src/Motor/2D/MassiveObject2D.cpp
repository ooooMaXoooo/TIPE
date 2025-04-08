#include "Motor\2D\MassiveObject2D.h"


namespace Motor {

    MassiveObject2D::MassiveObject2D(ldouble mass, glm::vec2 pos, GLFWwindow* parent_window, glm::vec2 initial_velocity, glm::vec2 initial_acceleration)
        : Core::MassiveObject(mass, parent_window),
          m_Pos(pos), m_Velocity(initial_velocity), m_Acceleration(initial_acceleration),
          m_Forces(glm::vec2(0))
    {

    }
    

    MassiveObject2D::MassiveObject2D(const MassiveObject2D& other)
        : Core::MassiveObject(other.m_Mass, other.m_Parent_window), m_Pos(other.m_Pos), m_Velocity(other.m_Velocity), m_Acceleration(other.m_Acceleration),
            m_Forces(other.m_Forces)
    {
    }

    MassiveObject2D::~MassiveObject2D()
    {
    }

    void MassiveObject2D::Update(float ts) {

    }

    void MassiveObject2D::OnRender() {

    }

    void MassiveObject2D::OnImGuiRender() {

    }

    void MassiveObject2D::ApplyForce(const glm::vec2& force)
    {
        m_Forces += force;
    }
}