#pragma once

#include "pch.h"
#include "Motor/Entity.h"

using ldouble = long double;

namespace Motor
{
	namespace Core
	{
		class MassiveObject : Entity, pgl::Drawable {
			ldouble m_Mass;

		public:

			MassiveObject();

			// virtual method of entity
			virtual void Update(float ts) override;

			// virtual method of drawable
			virtual void OnRender() override;
			virtual void OnImGuiRender() override;
		};
	};
};