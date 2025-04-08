#pragma once

#include "MassiveObject2D.h"
#include "pch.h"

namespace Motor
{
	class Planet2D : public MassiveObject2D
	{
	protected:
		const double m_Radius;

	public:
		Planet2D(ldouble mass, glm::vec2 center, double radius, GLFWwindow* parent_window, glm::vec2 initial_speed = glm::vec2(0, 0), glm::vec2 initial_acceleration = glm::vec2(0, 0));
		Planet2D(const Planet2D& planet);
		//Planet2D(Planet2D&& planet);

		~Planet2D();

		// virtual method of entity
		virtual void Update(float ts) override;

		// virtual method of drawable
		virtual void OnRender() override;
		virtual void OnImGuiRender() override;

	private:
		// private methods
	};
}
