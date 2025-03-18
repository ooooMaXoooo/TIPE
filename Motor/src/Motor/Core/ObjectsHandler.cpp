#include "Motor/Core/ObjectsHandler.h"

namespace Motor {
	namespace Core {

		ObjectsHandler::ObjectsHandler(pgl::Camera& activeCam)
		{
			m_ActiveCamera = std::make_unique<pgl::Camera>(activeCam);
		}
		

		void ObjectsHandler::AddObject(const pgl::Drawable& object)
		{
			std::unique_ptr<pgl::Drawable> obj_ptr = std::make_unique<pgl::Drawable>(object);

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