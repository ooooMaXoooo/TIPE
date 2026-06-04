#include <pch.h>
#include <SimuCore/statistics/Patate.h>

#include <SimuCore/structures/Rocket.h>

namespace SimuCore {
	namespace Statistics {
		Patate::Patate(size_t n_rockets) {
			m_rockets.reserve(n_rockets);
		}

		void Patate::AddRocket(const RocketData* rocket) {
			m_rockets.push_back(rocket);
		}

		void Patate::Merge(Patate& other) {
			m_rockets.insert(m_rockets.end(), other.m_rockets.begin(), other.m_rockets.end());
			other.m_rockets.clear();
		}

		RocketDataVectorSpace Patate::GetCentroid() const { // TODO : à check
			RocketDataVectorSpace barycentre;

			size_t size = m_rockets.size();

			for (int i = 0; i < size; i++) {
				barycentre += m_rockets[i]->GetAssociatedVector();
			}

			return barycentre / size;
		}

		double dissimilarity_mean(Patate const& p1, Patate const& p2) {
			double Dmean = 0;
			size_t total = 0;


			for (size_t i = 0; i < p1.m_rockets.size(); ++i) {
				for (size_t j = 0; j < p2.m_rockets.size(); ++j) {
					Dmean += distance(*(p1.m_rockets[i]), *(p2.m_rockets[j]));
					total++;
				}
			}

			return Dmean / total;
		}

		double dissimilarity_centroids(Patate const& p1, Patate const& p2) {
			RocketDataVectorSpace c1 = p1.GetCentroid();
			RocketDataVectorSpace c2 = p1.GetCentroid();

			return (c1 - c2).Length();
		}

		double dissimilarity(Patate const& p1, Patate const& p2) {
			return dissimilarity_mean(p1, p2);
		}
	}; // namespace Statistics
}; // namespace SimuCore