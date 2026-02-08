#pragma once

#include <pch.h>
#include "SimuCore/structures/Entity.h"
#include <SimuCore/constants.h>
#include <SimuCore/Units/all.h>

namespace SimuCore {
    namespace Structures {

        class Planet : public Entity {
            const char* m_Name;
            double m_radius;                        // Rayon physique de la planète (km)
            double m_muPlanet;                      // Paramètre gravitationnel de la planète (G * mass) (m^3/s²)
            double m_mass;                          // Masse de la planète (kg)

            double m_exobase;                       // altitude du début de l'exosphère (km)
            double m_maxAltitude;                   // altitude maximale de l'anneau dans lequel la fusée dot être capturée (km)

        public:
            /**
             * @brief Constructeur de Planet.
             * @param name Nom de la planète.
             * @param m Masse de la planète. (kg)
             * @param radius Rayon physique de la planète. (km)
             * @param exobase Altitude du début de l'exosphère. (km)
             * @param maxAltitude Altitude maximale de l'anneau dans lequel la fusée doit être capturée. (km)
             * @param p0 Position initiale (glm::dvec3). (AU)
             * @param v0 Vitesse initiale (glm::dvec3). (km/s)
             */
            Planet(const char* name, double mass, double radius, double exobase, double maxAltitude,
                glm::dvec3 p0 = glm::dvec3(0, 0, 0), glm::dvec3 v0 = glm::dvec3(0, 0, 0))
                :
                Entity(mass, p0, v0),
                m_Name(name),
                m_mass(mass),
                m_radius(radius),
                m_exobase(exobase),
                m_maxAltitude(maxAltitude),
				m_muPlanet(constants::G * mass)
            {
            }

            Planet& operator=(const Planet& other) {
                if (this != &other) {
                    Entity::operator=(other);
                    this->m_Name = other.m_Name;
                    this->m_radius = other.m_radius;
                    this->m_muPlanet = other.m_muPlanet;
                    this->m_mass = other.m_mass;
                    this->m_exobase = other.m_exobase;
                    this->m_maxAltitude = other.m_maxAltitude;
                }
                return *this;
			}
            

            std::string Name() const { return std::string(m_Name); }

            // --- FONCTIONS D'ORBITE AUTOMATIQUE (pour la phase finale) ---

            /// <summary>
            /// Renvoi le rayon cible pour l'orbite de capture, qui est {le rayon de la planète} + {la moyenne entre l'altitude de l'exobase et l'altitude maximale}.
            /// </summary>
            /// <returns>km</returns>
            double orbitRadius() const {
                return m_radius + ((m_exobase + m_maxAltitude) * 0.5);
            }

            /// <summary>
            /// Rayon minimum toléré pour l'orbite.
            /// </summary>
            /// <returns>km</returns>
            double minOrbitRadius() const {
                return m_exobase + m_radius;
            }

            
            /// <summary>
            /// Rayon maximum toléré pour l'orbite
            /// </summary>
            /// <returns>km</returns>
            double maxOrbitRadius() const {
                return m_radius + m_maxAltitude;
            }

            /** Vitesse nécessaire pour une orbite circulaire stable (Vitesse relative cible). */

            /// <summary>
            /// 
            /// </summary>
            /// <returns> km/s </returns>
            double targetCaptureVelocity() const {
                // Formule de la vitesse orbitale circulaire : V = sqrt(mu_planete / rayon_orbital)
                return meters_per_seconds_to_kilometers_per_seconds(
                    std::sqrt(m_muPlanet / kilometers_to_meters(orbitRadius()))
                    );
            }


            /// <summary>
			/// Vitesse relative cible attendue (égale à la vitesse de capture) (en km/h)
            /// </summary>
            /// <returns>km/s</returns>
            double orbitVelocity() const {
                return targetCaptureVelocity();
            }

            /// <summary>
            /// Vitesse minimal à l'extraction de la planète à partir d'un rayon initial en km
            /// </summary>
            /// <param name="initial_radius">km</param>
            /// <returns> vitesse d'extraction en km/s </returns>
            double extractionVelocity(double initial_radius) const {
                initial_radius = kilometers_to_meters(initial_radius); // passage en USI

                return meters_per_seconds_to_kilometers_per_seconds(
                    std::sqrt(2 * m_muPlanet / initial_radius)
                    );
            }


            /// <summary>
            /// Renvoi la masse de la planète en kg
            /// </summary>
            /// <returns>kg</returns>
            double getMass() const { return m_mass; }

            /// <summary>
            /// Renvoi le paramètre mu de la planète = masse * G (en m^3/s²)
            /// </summary>
            /// <returns>m^3/s² (USI) </returns>
            double getMu() const { return m_muPlanet; }

            /// <summary>
            /// Renvoi le rayon de la planète en km
            /// </summary>
            /// <returns>km</returns>
            double getRadius() const { return m_radius; }

            /// <summary>
            /// Renvoi l'altitude de l'exobase de la planète en km
            /// </summary>
            /// <returns>km</returns>
            double getExobase() const { return m_exobase; }

            /// <summary>
            /// Renvoi l'altitude maximal autorisé pour que la planète garde son influence par rapport à un astre central (en km)
            /// </summary>
            /// <returns>km</returns>
            double getMaxAltitude() const { return m_maxAltitude; }

            /// <summary>
            /// Calcul la vitesse angulaire de la planète autour d'un astre central
            /// </summary>
            /// <param name="distance_to_central_star">AU</param>
            /// <param name="mu_central_star">m^3/s²</param>
            /// <returns>vitesse angulaire rad/h </returns>
            double getAngularVelocity(double distance_to_central_star, double mu_central_star) const {
                distance_to_central_star = AU_to_meters(distance_to_central_star);
			    double omega = std::sqrt(
                    mu_central_star /(distance_to_central_star * distance_to_central_star * distance_to_central_star)
                ); // deuxième loi de Kepler

				// conversion de rad/s en rad/h
				return omega * 3600;
            };
		}; // class Planet
	} // namespace Structures
} // namespace SimuCore