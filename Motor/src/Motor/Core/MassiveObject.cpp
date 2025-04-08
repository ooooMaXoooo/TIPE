#include "Motor/Core/MassiveObject.h"

namespace Motor {
	namespace Core {
		MassiveObject::MassiveObject(ldouble mass, GLFWwindow* parent_window)
			: m_Mass(mass), Entity(), pgl::Drawable(parent_window)
		{

		}

	};
};