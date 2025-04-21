#pragma once

#include "basic/basic.h"
// drawable is already included since Circle.h include it
// same goes for <memory>

#include "Camera.h"

#include "pch.h"


#include "Motor/3D/Planet3D.h"


namespace Motor {
	namespace Core {

		/// <summary>
		/// A class that abstract the method to handle many objects in a single scene
		/// </summary>
		class ObjectsHandler
		{
		private:
			using Object = Planet3D;

			// the vector of objects to handle
			std::vector<std::unique_ptr<Object>> m_Objects;

			// the active camera
			pgl::Camera m_ActiveCamera;

			// the renderer of the block
			Renderer m_Renderer;


			////// Dear ImGui variables
			
			bool m_isAddingObject = false;

			glm::vec3 pos_to_add = glm::vec3(0);
			glm::vec3 rot_to_add = glm::vec3(0);
			glm::vec3 sca_to_add = glm::vec3(1);

			double mass_to_add = 1e5;
			double radius_to_add = 50;

			glm::vec3 v0_to_add = glm::vec3(0);
			glm::vec3 a0_to_add = glm::vec3(0);


		//// other
			GLFWwindow* m_Parent_window;

		public:
			ObjectsHandler(GLFWwindow* parent_window);


			void AddObject(const Object& object); // maybe return an id ?  --> 
			//    need to update the same id if objects are removed so its a bad approach  -->
			//    maybe not sa bad if we use nullptr when looking for an object


			void processInputs(float deltaTime);
			void UpdateParentWindow(GLFWwindow* parent_window) { m_Parent_window = parent_window; };

			////////////// object specific

			void UpdateViewMatrix();
			void UpdateObjects(float ts);
			void RenderObjects();
			void RenderImGuiObjects();

			void ImGuiRender();

			///////////// camera specific

		private :

			/*
			*	check for collisions --> octree ?
			*   if there is collision between 2 objects, we destroy both of them
			*/
			void HandleCollisions();

			bool DetectCollision(const Object& obj1, const Object& obj2) const;

		};
	};
};