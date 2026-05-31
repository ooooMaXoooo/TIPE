#pragma once

#include <pch.h>
#include <SimuCore/structures/Rocket.h>

#include <DataExport/RocketData.h>

namespace SimuCore {
	namespace Statistics {

		class Patate {
		public:
			Patate(size_t n_rockets);
			~Patate() = default;

			friend double dissimilarity(Patate const& p1, Patate const& p2);

			void AddRocket(const RocketData* rocket);

			void Merge(Patate& other);

			size_t size() const noexcept { return m_rockets.size(); }

		private:

			std::vector<const RocketData*> m_rockets; // faire attention à la durée de vie des pointeurs

		public :

			std::vector<const RocketData*> getRockets() const noexcept { return m_rockets; }
		};


		double dissimilarity(Patate const& p1, Patate const& p2);
	};
};