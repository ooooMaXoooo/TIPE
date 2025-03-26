#pragma once

#include "Motor/Core/MassiveObject.h"

namespace Motor
{
	class MassiveObject3D : public Core::MassiveObject
	{
    protected:
        glm::vec3 m_Pos;
        glm::vec3 m_Velocity;
        glm::vec3 m_Acceleration;
                
        glm::vec3 m_Forces;

    public:

        MassiveObject3D(ldouble mass, glm::vec3 pos, glm::vec3 initial_velocity = glm::vec3(0, 0, 0), glm::vec3 initial_acceleration = glm::vec3(0, 0, 0));
        MassiveObject3D(const MassiveObject3D& massiveObject_3D);

        virtual ~MassiveObject3D();

        // virtual method of entity
        virtual void Update(float ts) override;

        // virtual method of drawable
        virtual void OnRender() override;
        virtual void OnImGuiRender() override;


        void ApplyForce(const glm::vec3& force);
	};
};