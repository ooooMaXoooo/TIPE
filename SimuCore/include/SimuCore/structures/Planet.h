#pragma once

#include <pch.h>
#include "SimuCore/structures/Entity.h"
#include <SimuCore/constants.h>

namespace SimuCore {
	namespace Structures {

        class Planet : public Entity {
            const char* m_Name;
            double m_radius;                        // Rayon physique de la planète (m)
			double m_muPlanet;                      // Paramètre gravitationnel de la planète (G * mass) (m^3/s²)
			double m_mass;                          // Masse de la planète (kg)

			double m_exobase;                       // altitude du début de l'exosphère (m)
			double m_maxAltitude;                   // altitude maximale de l'anneau dans lequel la fusée dot être capturée (m)

			//double m_constanteDesAires;             // Constante des aires (m²/s)

        public:
            /**
             * @brief Constructeur de Planet.
             * @param name Nom de la planète.
			 * @param m Masse de la planète. (kg)
             * @param radius Rayon physique de la planète. (km)
             * @param exobase Altitude du début de l'exosphère. (km)
             * @param maxAltitude Altitude maximale de l'anneau dans lequel la fusée dot être capturée. (km)
			 * @param p0 Position initiale (glm::dvec3). (m)
			 * @param v0 Vitesse initiale (glm::dvec3). (m/s)
             */
			Planet(const char* name, double mass, double radius, double exobase, double maxAltitude,
                glm::dvec3 p0 = glm::dvec3(0, 0, 0), glm::dvec3 v0 = glm::dvec3(0, 0, 0))
                :
                Entity(mass, p0, v0),
                m_Name(name),
                m_radius(radius * 1e3),
				m_exobase(exobase * 1e3),
				m_maxAltitude(maxAltitude * 1e3),
                m_mass(mass)
            {
                SimuCore::constants::G;

                // Calcule le paramètre gravitationnel de la planète elle-même
                m_muPlanet = SimuCore::constants::G * mass;
            }

            std::string Name() const { return std::string(m_Name); }

            // --- FONCTIONS D'ORBITE AUTOMATIQUE (pour la phase finale) ---

            /** Rayon orbital cible. */
            double orbitRadius() const {
                return m_radius + ((m_exobase + m_maxAltitude) / 2.0);
            }

            /** Rayon minimum toléré pour l'orbite. */
            double minOrbitRadius() const {
                return m_exobase + m_radius;
            }

            /** Rayon maximum toléré pour l'orbite. */
            double maxOrbitRadius() const {
				return m_radius + m_maxAltitude;
            }

            /** Vitesse nécessaire pour une orbite circulaire stable (Vitesse relative cible). */
            double targetCaptureVelocity() const {
                // Formule de la vitesse orbitale circulaire : V = sqrt(mu_planete / rayon_orbital)
                return std::sqrt(m_muPlanet / orbitRadius());
            }

            double targetRadius() const {
                return m_radius + (0.5 * (m_exobase + m_maxAltitude));
            }

            /** Vitesse relative cible attendue (égale à la vitesse de capture). */
            double orbitVelocity() const {
                return targetCaptureVelocity();
            }

            /** Vitesse relative minimum tolérée. */
            double minOrbitVelocity() const {
				// On ne veut pas que la fusée soit trop lente pour éviter une chute vers la planète.
				// On calcule donc l'énergie potentielle effective avec un rayon égal au rayon minimal toléré.

                // Impossible ?
            }

            /** Vitesse relative maximum tolérée. */
            double maxOrbitVelocity() const {
            }

            /** Vitesse minimal à l'extraction de la planète à partir d'un rayon initial en m */
            double extractionVelocity(double initial_radius) {
                return std::sqrt(2 * m_muPlanet / initial_radius);
            }

            double getMu() const { return m_muPlanet; }
            double getRadius() const { return m_radius; }
        };
	}
}