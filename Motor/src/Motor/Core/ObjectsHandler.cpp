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

			m_Objects.push_back(std::move(obj_ptr));
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

			// handle forces and update every object at the same time
			const size_t NB_OBJECTS = m_Objects.size();
			for (int i = 0; i < NB_OBJECTS; i++)
			{
				std::unique_ptr<Object>& obj1 = m_Objects[i];
				for (int j = i + 1; j < NB_OBJECTS; j++)
				{
					//calculate the force of i on j
					std::unique_ptr<Object>& obj2 = m_Objects[j];
					const glm::vec3 direction_vector = obj1->GetPosition() - obj2->GetPosition();
					const glm::vec3 unit_vector = glm::normalize(direction_vector);
					glm::vec3 test = unit_vector;
					glm::vec3 force = static_cast<float>(_G * obj1->GetMass() * obj2->GetMass() / direction_vector.length()) * unit_vector;

					// apply this force on j and the opposite on i
					obj1->ApplyForce(force);
					obj2->ApplyForce(-force);
				}
				// we finished looking for all forces on i.
				// update i
				obj1->Update(ts);
			}

		}
	};
};