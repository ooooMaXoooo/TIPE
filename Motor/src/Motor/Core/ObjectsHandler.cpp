#include "Motor/Core/ObjectsHandler.h"

namespace Motor {
	namespace Core {

		ObjectsHandler::ObjectsHandler(GLFWwindow* parent_window) : m_Parent_window(parent_window)
		{
			m_ActiveCamera = pgl::Camera();
			m_Renderer = Renderer();
		}
		

		void ObjectsHandler::AddObject(const Object& object)
		{
			std::unique_ptr<Object> obj_ptr = std::make_unique<Object>(object);

			m_Objects.push_back(std::move(obj_ptr));
		}

		void ObjectsHandler::UpdateViewMatrix()
		{
			const auto& camera =  m_ActiveCamera.GetViewMatrix();

			for (auto& object : m_Objects)
			{
				object->UpdateViewMatrix(camera);
			}
		}


		void ObjectsHandler::UpdateObjects(float ts)
		{
			// collisions
			HandleCollisions();

			// handle forces and update every object at the same time
			const size_t NB_OBJECTS = m_Objects.size();
			for (int i = 0; i < NB_OBJECTS; i++)
			{
				std::unique_ptr<Object>& obj1 = m_Objects[i];
				for (int j = i + 1; j < NB_OBJECTS; j++)
				{
					//calculate the force of i on j
					std::unique_ptr<Object>& obj2 = m_Objects[j];
					const glm::vec3 direction_vector = obj1->GetPosition() - obj2->GetPosition(); // i - j
					const glm::vec3 unit_vector = glm::normalize(direction_vector);
					const ldouble distance = glm::length(direction_vector);

					glm::vec3 force(0);

					if(distance != 0)
						force = (static_cast<float>(_G * obj1->GetMass() * obj2->GetMass() / distance) * unit_vector);

					// apply this force on j and the opposite on i
					obj1->ApplyForce(-force); // i
					obj2->ApplyForce(force);  // j
				}
				// we have finish looking for all forces on i.
				// update i
				obj1->Update(ts);
			}
		}

		/// <summary>
		/// check for collisions with an octree
		/// if there is collision between 2 objects, we destroy both of them
		/// </summary>
		void ObjectsHandler::HandleCollisions()
		{

		}

		void ObjectsHandler::RenderObjects()
		{
			const glm::mat4& viewMatrix = m_ActiveCamera.GetViewMatrix();
			for (auto& object : m_Objects)
			{
				object->UpdateViewMatrix(viewMatrix);
				object->OnRender();
			}
		}

		void ObjectsHandler::RenderImGuiObjects()
		{
			for (auto& object : m_Objects)
			{
				object->OnImGuiRender();
			}
		}

		void ObjectsHandler::ImGuiRender()
		{
			ImGui::Begin("Object Handler");
			if (ImGui::Button("Add Object") && !m_isAddingObject)
			{
				m_isAddingObject = true;
			}

			if (m_isAddingObject)
			{
				ImGui::Begin("Attributes");



				if (ImGui::TreeNodeEx("Transform", ImGuiTreeNodeFlags_DefaultOpen))
				{

					ImGui::SliderFloat3("Position (x, y, z)", &pos_to_add[0], -20, 20);
					ImGui::SliderFloat3("Rotation (x, y, z)", &rot_to_add[0], -20, 20);
					ImGui::SliderFloat3("Scale (x, y, z)", &sca_to_add[0], 0, 5);

					ImGui::TreePop();
				}

				if (ImGui::TreeNodeEx("Physics"))
				{
					ImGui::InputDouble("Mass (kg)", &mass_to_add, 10, 1000);
					ImGui::InputDouble("Radius (km)", &radius_to_add, 1, 500);

					ImGui::InputFloat3("initial velocity (km/s)", &v0_to_add[0]);
					ImGui::InputFloat3("initial acceleration (km/s2)", &a0_to_add[0]);

					ImGui::TreePop();
				}



				if (ImGui::Button("Add"))
				{

					Object obj(mass_to_add, pos_to_add, radius_to_add, m_Parent_window, m_Renderer, v0_to_add, a0_to_add);
					AddObject(obj);
					m_isAddingObject = false;
				}

				if (ImGui::Button("Cancel"))
				{
					m_isAddingObject = false;
				}

				ImGui::End();
			}

			ImGui::SeparatorText("Objects");

			// render all object in a list
			constexpr uint nbColumns = 1;
			if (ImGui::BeginTable("table1", nbColumns))
			{

				const size_t nbRows = m_Objects.size();
				for (size_t row = 0; row < nbRows; row++)
				{
					ImGui::TableNextRow();


					ImGui::TableSetColumnIndex(0);
					ImGui::Text("Object %d\t(x, y, z) : (%.6f, %.6f, %.6f)", row, m_Objects[row]->GetPosition()[0], m_Objects[row]->GetPosition()[1], m_Objects[row]->GetPosition()[2]);
					//ImGui::Text("(x, y, z) : (%.6f, %.6f, %.6f)", m_Objects[row]->GetPosition()[0], m_Objects[row]->GetPosition()[1], m_Objects[row]->GetPosition()[2]);
				}
				ImGui::EndTable();
			}

			ImGui::End();
		}


		void ObjectsHandler::processInputs(float deltaTime)
		{
			// handle camera movements
			glm::vec3 direction = glm::vec3(0);
			if (ImGui::IsKeyDown(ImGuiKey_Z))
			{
				direction += glm::vec3(0.0f, 1.0f, 0.0f);
			} 
			if (ImGui::IsKeyDown(ImGuiKey_S))
			{
				direction += glm::vec3(0.0f, -1.0f, 0.0f);
			}
			if (ImGui::IsKeyDown(ImGuiKey_Q))
			{
				direction += glm::vec3(-1.0f, 0.0f, 0.0f);
			}
			if (ImGui::IsKeyDown(ImGuiKey_D))
			{
				direction += glm::vec3(1.0f, 0.0f, 0.0f);
			}
			if (ImGui::IsKeyDown(ImGuiKey_A))
			{
				direction += (glm::vec3(0.0f, 0.0f, -1.0f));
			}
			if (ImGui::IsKeyDown(ImGuiKey_E))
			{
				direction += (glm::vec3(0.0f, 0.0f, 1.0f));
			}

			const float DIRECTION_LENGTH = glm::length(direction);

			if (DIRECTION_LENGTH)
				direction = glm::normalize(direction);

			m_ActiveCamera.ProcessKeyboard(direction, deltaTime);
		}
	};
};