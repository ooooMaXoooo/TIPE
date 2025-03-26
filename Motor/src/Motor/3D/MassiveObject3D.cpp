#include "Motor/3D/MassiveObject3D.h"

namespace Motor {
	MassiveObject3D::MassiveObject3D(ldouble mass, glm::vec3 pos, glm::vec3 v0, glm::vec3 a0)
		: Core::MassiveObject(mass), m_Pos(pos), m_Velocity(v0), m_Acceleration(a0)
	{

	}

	MassiveObject3D::MassiveObject3D(const MassiveObject3D& massiveObject_3D)
		:   Core::MassiveObject(massiveObject_3D.m_Mass),
			m_Pos(massiveObject_3D.m_Pos),
			m_Velocity(massiveObject_3D.m_Velocity),
			m_Acceleration(massiveObject_3D.m_Acceleration)
	{
	}

	// virtual method of entity
	void MassiveObject3D::Update(float ts)
	{

	}

	// virtual method of drawable
	void MassiveObject3D::OnRender()
	{

	}

	void MassiveObject3D::OnImGuiRender()
	{

	}


	void MassiveObject3D::ApplyForce(const glm::vec3& force)
	{
		m_Forces += force;
	}

};