#include "Motor/3D/MassiveObject3D.h"

namespace Motor {
	MassiveObject3D::MassiveObject3D(ldouble mass, glm::vec3 pos, GLFWwindow* parent_window, Renderer& renderer, glm::vec3 v0, glm::vec3 a0)
		: Core::MassiveObject(mass, parent_window, renderer), m_Pos(pos), m_Velocity(v0), m_Acceleration(a0),
		m_Forces(0)
	{
	}

	MassiveObject3D::MassiveObject3D(const MassiveObject3D& massiveObject_3D)
		:   Core::MassiveObject(massiveObject_3D.m_Mass, massiveObject_3D.m_Parent_window, massiveObject_3D.m_renderer),
			m_Pos(massiveObject_3D.m_Pos),
			m_Velocity(massiveObject_3D.m_Velocity),
			m_Acceleration(massiveObject_3D.m_Acceleration),
			m_Forces(0)
	{
	}

	MassiveObject3D::~MassiveObject3D() {

	}

	// virtual method of entity
	void MassiveObject3D::Update(float ts)
	{
		// look for forces  --> job done by the objectHandler
		// apply the resultant to the acceleration
		// ( WARNING : do not ADD the forces to the acceleration / velocity, Newton's laws is   acceleration IS EQUAL TO the sum of forces)
		m_Acceleration = static_cast<float>(1 / m_Mass) * m_Forces;

		// check for potential collisions

		// if there is no collisions, let the physics make its job
		// 
		// if there is collisions, delete both planets
		//  ---->  the desctructor need to clear everything potentially including rendering stuff
		//  ----> it's handle by the objectHandler

		// apply half the acceleration to velocity
		const glm::vec3 half_acc = static_cast<float>(ts * 0.5) * m_Acceleration;
		m_Velocity += half_acc;

		// apply velocity to position
		m_Pos += m_Velocity * ts;

		// apply the remaining half of acceleration to velocity in order to correct the mean
		m_Velocity += half_acc;

		// reset forces to avoid old forces to interact at the next frame
		m_Forces = glm::vec3(0.0);
	}


	void MassiveObject3D::UpdateFirstPart(float ts) {
		m_Acceleration = static_cast<float>(1/m_Mass) * m_Forces;
		m_Velocity += static_cast<float>(ts*0.5) * m_Acceleration;
		m_Pos += ts * m_Velocity;

		m_Forces = glm::vec3(0.0);
	}

	void MassiveObject3D::UpdateSecondPart(float ts) {
		m_Acceleration = static_cast<float>(1 / m_Mass) * m_Forces;
		m_Velocity += static_cast<float>(ts * 0.5) * m_Acceleration;
		m_Forces = glm::vec3(0.0);
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