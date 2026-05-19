#include <pch.h>
#include <SimuCore/statistics/Patate.h>

#include <SimuCore/structures/Rocket.h>

namespace SimuCore {
	namespace Statistics {
		Patate::Patate(size_t n_rockets) {
			m_rockets.reserve(n_rockets);
		}

		void Patate::AddRocket(const SimuCore::Structures::Rocket* rocket, size_t indice) {
			m_rockets.push_back({ rocket, indice });
		}

		void Patate::Merge(Patate& other) {
			m_rockets.insert(m_rockets.end(), other.m_rockets.begin(), other.m_rockets.end());
			other.m_rockets.clear();
		}

		double dissimilarity(Patate const& p1, Patate const& p2) {
			/* on utilise la dissimilarité Dmin : 
			* 
			* Pour tout couple de fusées (f1, f2) avec f1 dans p1 et f2 dans p2, on calcule la distance d(f1, f2) entre les deux fusées.
			* En notant Dmin la plus petite de ces distances, on définit la dissimilarité entre p1 et p2 comme étant Dmin.
			* 
			* */

			double Dmin = std::numeric_limits<double>::max();

			for (size_t i = 0; i < p1.m_rockets.size(); ++i) {
				for (size_t j = 0; j < p2.m_rockets.size(); ++j) {
					double d = SimuCore::Structures::distance(*(p1.m_rockets[i].first), *(p2.m_rockets[j].first));

					if (d < Dmin) {
						Dmin = d;
					}
				}
			}

			return Dmin;
		}
	}; // namespace Statistics
}; // namespace SimuCore