#include "Motor/Core/ObjectsHandler.h"

namespace Motor {
	namespace Core {

		ObjectsHandler::ObjectsHandler(pgl::Camera& activeCam)
			: m_ActiveCamera(activeCam)
		{
			//m_ActiveCamera = std::make_unique<pgl::Camera>(activeCam);
		}
		

		void ObjectsHandler::AddObject(const Object& object)
		{
			std::unique_ptr<Object> obj_ptr = std::make_unique<Object>(object);

			m_Objects.push_back(obj_ptr);
		}

		void ObjectsHandler::UpdateViewMatrix()
		{
			const auto camera =  m_ActiveCamera.GetViewMatrix();

			for (auto& object : m_Objects)
			{
				object->UpdateViewMatrix(camera);
			}
		}


		void ObjectsHandler::UpdateObjects(float ts)
		{
			// detect collisions
			// destroys objects

			// handle forces


			for (auto& object : m_Objects)
			{
				object->Update(ts);
			}
		}
	};
};