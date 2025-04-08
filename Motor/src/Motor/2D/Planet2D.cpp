#include "Motor\2D\Planet2D.h"

namespace Motor {
    Planet2D::Planet2D(ldouble mass, glm::vec2 center, double radius, GLFWwindow* parent_window, glm::vec2 v0, glm::vec2 a0)
        : MassiveObject2D{ mass, center, parent_window, v0, a0 }, m_Radius(radius)
    {
    }

    Planet2D::Planet2D(const Planet2D& p)
        : m_Radius(p.m_Radius), MassiveObject2D{ p.m_Mass, p.m_Pos, p.m_Parent_window, p.m_Velocity, p.m_Acceleration }
    {
    }

    Planet2D::~Planet2D()
    {
        // need to handle everything
    }

    void Planet2D::Update(float ts)
    {
        // look for forces  --> job done by the objectHandler
        // apply the resultant to the acceleration
        // ( WARNING : do not ADD the forces to the acceleration / velocity, Newton's laws is   acceleration IS EQUAL TO the sum of forces)
        m_Acceleration = m_Forces;

        // check for potential collisions

        // if there is no collisions, let the physics make its job
        // 
        // if there is collisions, delete both planets
        //  ---->  the desctructor need to clear everything potentially including rendering stuff
        //  ----> it's handle by the objectHandler

        // apply half the acceleration to velocity
        const glm::vec2 half_acc = m_Acceleration * ts * 0.5f;
        m_Velocity += half_acc;

        // apply velocity to position
        m_Pos += m_Velocity * ts;

        // apply the remaining half of acceleration to velocity in order to correct the mean
        m_Velocity += half_acc;

        // done ?
    }

    void Planet2D::OnRender()
    {
        // rendering stuff
    }

    void Planet2D::OnImGuiRender()
    {

    }
};