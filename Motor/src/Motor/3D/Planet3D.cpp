#include "Motor/3D/Planet3D.h"

namespace Motor {

    Planet3D::Planet3D(ldouble mass, glm::vec3 center, double radius, glm::vec3 v0, glm::vec3 a0)
        : MassiveObject3D{ mass, center, v0, a0 }, m_Radius(radius)
    {
    }

    Planet3D::Planet3D(const Planet3D& p)
        : m_Radius(p.m_Radius), MassiveObject3D{ p.m_Mass, p.m_Pos, p.m_Velocity, p.m_Acceleration }
    {
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
        m_Acceleration = m_Forces;

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

        // done ?
    }

    void Planet3D::OnRender()
    {
        // rendering stuff
    }

    void Planet3D::OnImGuiRender()
    {

    }

};