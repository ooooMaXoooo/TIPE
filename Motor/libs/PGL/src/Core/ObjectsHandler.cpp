#include "Core/ObjectsHandler.h"

namespace pgl {
	namespace Core {

		ObjectsHandler::ObjectsHandler(Camera& activeCam)
		{
			m_ActiveCamera = std::make_unique<Camera>(activeCam);
		}
	};
};