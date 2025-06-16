#include "Drawable.h"

namespace pgl {

	Drawable::Drawable(GLFWwindow* parent_window, Renderer& renderer) : 
		m_Parent_window(parent_window),
		m_Model(1.0f), m_View (1.0f), m_Proj(1.0f),
		m_renderer(renderer)
	{
		m_MVP = m_Proj * m_View * m_Model;
	}

	void Drawable::OnRender() {

	}

	void Drawable::OnImGuiRender() {

	}

	void Drawable::UpdateViewMatrix(const glm::mat4& view)
	{
		m_View = view;
	}
};

