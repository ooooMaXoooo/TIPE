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
			// On va définir une "distance" entre deux fusées sans utiliser de normes.
			// 
			// Afin de coller à la version python, nous avons besoin de connaître le tof.
			// Il nous est donc nécessaire de faire une simulation de la trajectoire pour en déduire le tof.
			// Cependant, si nous faisons une simulation complète, autant faire une distance de Fréchet, et cela serait très coûteux en temps de calcul.
		}

	}; // namespace Structures
}; // namespace SimuCore