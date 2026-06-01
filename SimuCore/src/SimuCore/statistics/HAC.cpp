#include <pch.h>
#include <SimuCore/statistics/HAC.h>
#include <SimuCore/statistics/Patate.h>

namespace SimuCore {
	namespace Statistics {

			std::vector<Cluster> SimuCore::Statistics::HAC(const std::vector<std::unique_ptr<RocketData>>& rockets)
			{
				constexpr double dissimilarity_threshold = 100; // à ajuster en fonction des résultats

				// on initialise les patates avec une fusée chacune
				std::vector<Patate> patates;
				patates.reserve(rockets.size());

				for (const auto& rocket : rockets) {
					Patate patate(5); // TODO : à modifier
					patate.AddRocket(rocket.get());

					patates.emplace_back(patate);
				}

				// tant que la dissimilitude entre les patates est inférieure au seuil, on fusionne les deux patates les plus proches
				double dissim = std::numeric_limits<double>::lowest();

				while (dissim < dissimilarity_threshold && patates.size() > 1) {

					size_t idx1 = 0, idx2 = 0;
					dissim = std::numeric_limits<double>::max();

					// on parcourt toutes les paires de patates pour trouver les deux patates les plus proches
					for (size_t i = 0; i < patates.size(); i++) {
						for (size_t j = i + 1; j < patates.size(); j++) {
							double d = dissimilarity(patates[i], patates[j]);
							if (d < dissim) {
								dissim = d;
								idx1 = i;
								idx2 = j;
							}
						}
					}

					if (dissim < dissimilarity_threshold) {
						patates[idx1].Merge(patates[idx2]);
						patates.erase(patates.begin() + idx2);
					}
				}

				// on retourne les clusters composés des pointeurs vers les individus de la population
				std::vector<Cluster> clusters;

				for (const auto& patate : patates) {
					Cluster cluster;
					cluster.reserve(patate.size());
					
					for (const auto& rocket_ptr : patate.getRockets()) {
						cluster.push_back(rocket_ptr->GetId());
					}
					clusters.emplace_back(cluster);
				}
				return clusters;
			}

	} // namespace Statistics
}; // namespace SimuCore
