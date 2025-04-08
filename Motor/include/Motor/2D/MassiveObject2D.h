#pragma once

#include "Motor\Core\MassiveObject.h"

/*
#include "Motor/Entity.h"
#include "Vector.h"
*/

namespace Motor {
    class MassiveObject2D : public Core::MassiveObject
    {
    protected:
        glm::vec2 m_Pos;
        glm::vec2 m_Velocity;
        glm::vec2 m_Acceleration;

        glm::vec2 m_Forces;

    public:

        MassiveObject2D(ldouble mass, glm::vec2 pos, GLFWwindow* parent_window, glm::vec2 initial_velocity = glm::vec2(0, 0), glm::vec2 initial_acceleration = glm::vec2(0, 0));
        MassiveObject2D(const MassiveObject2D& massiveObject_2D);

        virtual ~MassiveObject2D();

        // virtual method of entity
        virtual void Update(float ts) override;

        // virtual method of drawable
        virtual void OnRender() override;
        virtual void OnImGuiRender() override;


        void ApplyForce(const glm::vec2& force);

    };

}

