#include <pch.h>
#include <SimuCore\structures\Rocket.h>


namespace SimuCore {
	namespace Structures {

		Impulsion operator-(const Impulsion& I1, const Impulsion& I2) {
			return Impulsion(I1.velocity - I2.velocity);
		}

		Rocket operator-(const Rocket& r1, const Rocket& r2) {

			std::vector<std::pair<Impulsion, double>> impulsions;
			impulsions.reserve(r1.m_Impulsions.size());

			for (int i = 0; i < r1.m_Impulsions.size(); i++) {
				auto& [imp1, date1] = r1.m_Impulsions[i];
				auto& [imp2, date2] = r2.m_Impulsions[i];

				impulsions.emplace_back(imp1 - imp2, date1 - date2);
			}


			return Rocket(r1.lifetime, impulsions, r1.mass, r1.m_Vitesse_ejection_gaz, r1.position - r2.position, r1.velocity - r2.velocity);
		}

		double distance(const Rocket& r1, const Rocket& r2) {
			return (r1 - r2).Norme();
		}

	}; // namespace Structures
}; // namespace SimuCore