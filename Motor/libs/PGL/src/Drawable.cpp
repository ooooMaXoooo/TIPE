#include "Drawable.h"

namespace pgl {

	void Drawable::OnRender() {

	}

	void Drawable::OnImGuiRender() {

	}

	void Drawable::UpdateViewMatrix(glm::mat4 view)
	{
		m_View = view;
	}
};

