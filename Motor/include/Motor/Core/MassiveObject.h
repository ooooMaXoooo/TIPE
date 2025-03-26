#pragma once

#include "pch.h"
#include "Motor/Entity.h"

#define _G 6.67e-11

using ldouble = long double;

namespace Motor
{
	namespace Core
	{
		class MassiveObject : public Entity, public pgl::Drawable {
		protected:
			const ldouble m_Mass;

		public:

			MassiveObject(ldouble mass);

			// virtual method of entity
			virtual void Update(float ts) override;

			// virtual method of drawable
			virtual void OnRender() override;
			virtual void OnImGuiRender() override;
		};
	};
};