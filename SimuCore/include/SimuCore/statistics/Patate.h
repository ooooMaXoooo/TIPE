#pragma once

#include <pch.h>
#include <SimuCore/structures/Rocket.h>

namespace SimuCore {
	namespace Statistics {

		class Patate {
		public:
			Patate(size_t n_rockets);
			~Patate() = default;

			friend double dissimilarity(Patate const& p1, Patate const& p2);

			void AddRocket(const SimuCore::Structures::Rocket* rocket, size_t indice);

			void Merge(Patate& other);

			size_t size() const noexcept { return m_rockets.size(); }

		private:

			std::vector<std::pair<const SimuCore::Structures::Rocket*, size_t>> m_rockets; // faire attention à la durée de vie des pointeurs

		public :

			std::vector<std::pair<const SimuCore::Structures::Rocket*, size_t>> getRockets() const noexcept { return m_rockets; }
		};


		double dissimilarity(Patate const& p1, Patate const& p2);
	};
};