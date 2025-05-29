#pragma once

#include "Motor/3D/MassiveObject3D.h"

namespace Motor {
	class Planet3D : public MassiveObject3D
	{
	protected:
		const double m_Radius;

	public:

		Planet3D(ldouble mass, glm::vec3 center, double radius, GLFWwindow* parent_window, Renderer& renderer, glm::vec3 initial_speed = glm::vec3(0, 0, 0), glm::vec3 initial_acceleration = glm::vec3(0, 0, 0));
		Planet3D(const Planet3D& planet);
		//Planet2D(Planet2D&& planet);

		~Planet3D();

		// virtual method of entity
		virtual void Update(float ts) override;
		virtual void UpdateFirstPart(float ts) override;
		virtual void UpdateSecondPart(float ts) override;


		// virtual method of drawable
		virtual void OnRender() override;
		virtual void OnImGuiRender() override;

		double GetRadius() const { return m_Radius; }

	private:
		// private methods
		void UpdateTransform();

		void InitPlanet();
	};
};