#include "Motor\2D\MassiveObject2D.h"


namespace Motor {

    MassiveObject2D::MassiveObject2D(ldouble mass, glm::vec2 pos, glm::vec2 initial_velocity, glm::vec2 initial_acceleration)
        : Core::MassiveObject(mass),
          m_Pos(pos), m_Velocity(initial_velocity), m_Acceleration(initial_acceleration)
    {

    }
    

    MassiveObject2D::MassiveObject2D(const MassiveObject2D& other)
        : Core::MassiveObject(other.m_Mass), m_Pos(other.m_Pos), m_Velocity(other.m_Velocity), m_Acceleration(other.m_Acceleration)
    {
    }

    MassiveObject2D::~MassiveObject2D()
    {
    }

    void MassiveObject2D::ApplyForce(const glm::vec2& force)
    {
        m_Forces += force;
    }
}