#include "Motor/Core/ObjectsHandler.h"

namespace Motor {
	namespace Core {

		ObjectsHandler::ObjectsHandler(Camera& activeCam)
		{
			m_ActiveCamera = std::make_unique<Camera>(activeCam);
		}
		

		void ObjectsHandler::AddObject(const Drawable& object)
		{
			std::unique_ptr<Drawable> obj_ptr = std::make_unique<Drawable>(object);

			m_Objects.push_back(obj_ptr);
		}

		void ObjectsHandler::UpdateViewMatrix()
		{
			const auto camera =  m_ActiveCamera->GetViewMatrix();

			for (auto& object : m_Objects)
			{
				object->UpdateViewMatrix(camera);
			}
		}
		void ObjectsHandler::UpdateObjects()
		{
		}
	};
};