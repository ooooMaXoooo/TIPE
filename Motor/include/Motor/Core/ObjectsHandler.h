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
			const pgl::Camera& m_ActiveCamera;

		public:
			ObjectsHandler(pgl::Camera& activeCam);


			void AddObject(const Object& object); // maybe return an id ?  --> 
			//    need to update the same id if objects are removed so its a bad approach  -->
			//    maybe not sa bad if we use nullptr when looking for an object

			void UpdateViewMatrix();
			void UpdateObjects(float ts);
		};
	};
};