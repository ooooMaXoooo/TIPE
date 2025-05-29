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
			const ldouble softening = 1e-3; // valeur ajustable selon les unit�s

			for (int i = 0; i < NB_OBJECTS; i++)
			{
				std::unique_ptr<Object>& obj1 = m_Objects[i];

				for (int j = i + 1; j < NB_OBJECTS; j++)
				{
					std::unique_ptr<Object>& obj2 = m_Objects[j];

					// Vecteur direction de i vers j
					const glm::vec3 direction_vector = obj1->GetPosition() - obj2->GetPosition();

					// Distance au carr� avec softening
					const ldouble distance_squared = glm::dot(direction_vector, direction_vector) + softening * softening;
					const ldouble distance = std::sqrt(distance_squared);

					// Vecteur unitaire
					const glm::vec3 unit_vector = direction_vector / static_cast<float>(distance);

					// Force gravitationnelle
					glm::vec3 force = static_cast<float>(_G * obj1->GetMass() * obj2->GetMass() / distance_squared) * unit_vector;

					// Application de la force : action = -r�action
					obj1->ApplyForce(-force);
					obj2->ApplyForce(force);
				}

				// Une fois toutes les forces appliqu�es sur obj1, on met � jour sa position
				obj1->UpdateFirstPart(ts);
			}

			// tous les objets ont �t� mis � jours, il faut maintenant faire la seconde partie de la m�thode de Leapfrog
			for (int i = 0; i < NB_OBJECTS; i++)
			{
				std::unique_ptr<Object>& obj1 = m_Objects[i];

				for (int j = i + 1; j < NB_OBJECTS; j++)
				{
					std::unique_ptr<Object>& obj2 = m_Objects[j];

					// Vecteur direction de i vers j
					const glm::vec3 direction_vector = obj1->GetPosition() - obj2->GetPosition();

					// Distance au carr� avec softening
					const ldouble distance_squared = glm::dot(direction_vector, direction_vector) + softening * softening;
					const ldouble distance = std::sqrt(distance_squared);

					// Vecteur unitaire
					const glm::vec3 unit_vector = direction_vector / static_cast<float>(distance);

					// Force gravitationnelle
					glm::vec3 force = static_cast<float>(_G * obj1->GetMass() * obj2->GetMass() / distance_squared) * unit_vector;

					// Application de la force : action = -r�action
					obj1->ApplyForce(-force);
					obj2->ApplyForce(force);
				}

				// Une fois toutes les forces appliqu�es sur obj1, on met � jour sa position
				obj1->UpdateSecondPart(ts);
			}
		}

		void ObjectsHandler::HandleCollisions()
		{
			// Le facteur de masse pour d�terminer quel objet est "tr�s peu massique"
			const double massFactor = 1e4;

			// On utilise une taille dynamique pour �viter des probl�mes avec les indices apr�s suppression
			size_t NB_OBJECTS = m_Objects.size();
			for (size_t i = 0; i < NB_OBJECTS; ++i) {
				std::unique_ptr<Object>& obj1 = m_Objects[i];

				for (size_t j = i + 1; j < NB_OBJECTS; ++j) {
					std::unique_ptr<Object>& obj2 = m_Objects[j];

					// V�rifier si une collision a eu lieu (avec ta m�thode de d�tection de collision, par exemple)
					if (DetectCollision(*obj1, *obj2)) {
						// Calcul des masses relatives pour d�terminer quel objet est dominant
						double massRatio = obj1->GetMass() / obj2->GetMass();

						// Traiter la collision selon le rapport de masse
						if (massRatio >= massFactor) {
							// obj2 est tr�s l�ger, on applique un impact l�ger sur obj1 et on d�truit obj2
							obj1->ApplyForce((obj2->GetVelocity() - obj1->GetVelocity()) * 0.01f);
							m_Objects.erase(m_Objects.begin() + j);  // Supprimer obj2 du tableau
							--j; // Ajuster l'index j apr�s la suppression
							NB_OBJECTS = m_Objects.size(); // Mettre � jour la taille apr�s suppression
							continue; // Passer � la prochaine paire
						}
						else if (massRatio <= 1.0 / massFactor) {
							// obj1 est tr�s l�ger, on applique un impact l�ger sur obj2 et on d�truit obj1
							obj2->ApplyForce((obj1->GetVelocity() - obj2->GetVelocity()) * 0.01f);
							m_Objects.erase(m_Objects.begin() + i);  // Supprimer obj1 du tableau
							--i; // Ajuster l'index i apr�s la suppression
							NB_OBJECTS = m_Objects.size(); // Mettre � jour la taille apr�s suppression
							break; // Passer � l'objet suivant
						}
						else {
							// Les deux objets ont des masses similaires, on les d�truit tous les deux
							m_Objects.erase(m_Objects.begin() + j);  // Supprimer obj2
							m_Objects.erase(m_Objects.begin() + i);  // Supprimer obj1
							--j; // Ajuster l'index j apr�s la suppression de obj2
							NB_OBJECTS = m_Objects.size(); // Mettre � jour la taille apr�s suppression
							break; // Passer � l'objet suivant
						}
					}
				}
			}
		}

		bool ObjectsHandler::DetectCollision(const Object& obj1, const Object& obj2) const
		{
			// R�cup�rer les positions des objets
			const glm::vec3 pos1 = obj1.GetPosition();
			const glm::vec3 pos2 = obj2.GetPosition();

			// Calculer la distance entre les deux objets
			float distance = glm::length(pos1 - pos2);

			// R�cup�rer les rayons des objets (suppos�s �tre des sph�res)
			float radius1 = obj1.GetRadius();
			float radius2 = obj2.GetRadius();

			// V�rifier si la distance est inf�rieure � la somme des rayons
			if (distance < (radius1 + radius2)) {
				return true; // Collision d�tect�e
			}
			return false; // Pas de collision
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
					ImGui::InputFloat3("Position (x, y, z)", &pos_to_add[0]);
					ImGui::InputFloat3("Rotation (x, y, z)", &rot_to_add[0]);
					ImGui::InputFloat3("Scale (x, y, z)", &sca_to_add[0]);

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

			// camera component
			ImGui::SeparatorText("Camera");
			m_ActiveCamera.OnImGuiRender();

			// renderer component
			m_Renderer.OnImGuiRender();

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