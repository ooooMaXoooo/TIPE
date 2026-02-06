#pragma once

#include <pch.h>
#include "GenerationStats.h"
#include <DataExport/Snapshot.h>

#include <GaLib\genetic.hpp>

#include <SimuCore.h>

class ScalarAccumulator;
struct VectorAccumulator;
class GenerationAccumulator;

inline dataExport::GenerationValidityStats calculate_validity_stats(
	const std::vector<dataExport::IndividualValidity>& ind_valid,
	const std::vector<dataExport::TrajectoryValidity>& traj_valid
);

// calcule les stats d'un vecteur de doubles
inline dataExport::ScalarStats calculate_scalar_stats(const std::vector<double>& data);

template<genetic::ConfigConcept ConfigType>
void calculate_genes_stats(
	const std::vector<genetic::Individu<ConfigType>>& population,
	GenerationAccumulator& accumulator
);


template<genetic::ConfigConcept ConfigType>
dataExport::TrajectorySoA calculate_trajectory(
	const genetic::Individu<ConfigType>& individual,
	SimuCore::Systems::AdaptedSystem system
);

template<genetic::ConfigConcept ConfigType>
void push_population_to_accumulator(
	GenerationAccumulator& accumulator,
	const std::vector<genetic::Individu<ConfigType>>& population,
	const std::vector<dataExport::IndividualValidity>& ind_valid,
	const std::vector<dataExport::TrajectoryValidity>& traj_valid
);

template<genetic::ConfigConcept ConfigType>
void calculate_population_validity(
	const std::vector<genetic::Individu<ConfigType>>& population,
	const SimuCore::Systems::AdaptedSystem& system,
	std::vector<dataExport::IndividualValidity>& ind_valid,
	std::vector<dataExport::TrajectoryValidity>& traj_valid
);



































class ScalarAccumulator {
public:
    void reserve(size_t n) { m_values.reserve(n); }

    void push(double v) {
        m_values.push_back(v);

        m_min = std::min(m_min, v);
        m_max = std::max(m_max, v);

        // Welford online variance
        m_count++;
        double delta = v - m_mean;
        m_mean += delta / m_count;
        double delta2 = v - m_mean;
        m_M2 += delta * delta2;
    }

    void finalize(int bin_count);

    const dataExport::ScalarStats& stats() const { return m_stats; }
    const dataExport::Histogram& histogram() const { return m_hist; }

private:
    std::vector<double> m_values;

    size_t m_count = 0;
    double m_mean = 0.0;
    double m_M2 = 0.0;
    double m_min = std::numeric_limits<double>::infinity();
    double m_max = -std::numeric_limits<double>::infinity();

    dataExport::ScalarStats m_stats;
    dataExport::Histogram m_hist;
};

struct VectorAccumulator {
    ScalarAccumulator x, y, z;

    void reserve(size_t n) {
        x.reserve(n);
        y.reserve(n);
        z.reserve(n);
    }

    // remplacement de glm::dvec3 par std::array<double,3>
    void push(const std::array<double, 3>& v) {
        x.push(v[0]);
        y.push(v[1]);
        z.push(v[2]);
    }

    void finalize(int bins) {
        x.finalize(bins);
        y.finalize(bins);
        z.finalize(bins);
    }

	const dataExport::VectorStats stats() const {
		dataExport::VectorStats stats;
		stats.x = x.stats();
		stats.y = y.stats();
		stats.z = z.stats();
		return stats;
	}
};


class GenerationAccumulator {
public:
    explicit GenerationAccumulator(uint8_t impulse_count)
        : impulses(impulse_count),
        delta_times(impulse_count) {
    }

    void reserve(size_t pop) {
        initial_position.reserve(pop);
        scores.reserve(pop);
        for (auto& i : impulses) i.reserve(pop);
        for (auto& d : delta_times) d.reserve(pop);
        individual_validity.reserve(pop);
        trajectory_validity.reserve(pop);
    }

    void push_individual(
        const std::array<double, 3>& pos,
        std::span<const std::array<double, 3>> imp,
        std::span<const double> dt,
        double fitness,
        dataExport::IndividualValidity ind_valid,
        dataExport::TrajectoryValidity traj_valid
    ) {
        initial_position.push(pos);
        scores.push(fitness);

        for (size_t i = 0; i < imp.size(); ++i) {
            impulses[i].push(imp[i]);
            delta_times[i].push(dt[i]);
        }

        individual_validity.push_back(ind_valid);
        trajectory_validity.push_back(traj_valid);
    }

    dataExport::Snapshot finalize(size_t gen, const dataExport::BestIndividualUpdate& best_update);

	template<genetic::ConfigConcept ConfigType>
	void push_population(
		const std::vector<genetic::Individu<ConfigType>>& population,
		const SimuCore::Systems::AdaptedSystem& system
	) {
		reserve(population.size());

		std::vector<dataExport::IndividualValidity> ind_valid(population.size());
		std::vector<dataExport::TrajectoryValidity> traj_valid(population.size());

		calculate_population_validity(population, system, ind_valid, traj_valid);

		push_population_to_accumulator(*this, population, ind_valid, traj_valid);
	}


private:
    VectorAccumulator initial_position;
    std::vector<VectorAccumulator> impulses;
    std::vector<ScalarAccumulator> delta_times;
    ScalarAccumulator scores;

    std::vector<dataExport::IndividualValidity> individual_validity;
    std::vector<dataExport::TrajectoryValidity> trajectory_validity;
};




inline dataExport::GenerationValidityStats calculate_validity_stats(
    const std::vector<dataExport::IndividualValidity>& ind_valid,
    const std::vector<dataExport::TrajectoryValidity>& traj_valid
) {
    dataExport::GenerationValidityStats stats;
    stats.population_size = ind_valid.size();

    for (auto v : ind_valid) {
        if (v == dataExport::IndividualValidity::Valid)
            stats.valid_individuals++;
        else {
            stats.invalid_individuals++;
            stats.individual_failure_counts[v]++;
        }
    }

    for (auto v : traj_valid) {
        if (v == dataExport::TrajectoryValidity::Valid)
            stats.valid_trajectories++;
        else {
            stats.invalid_trajectories++;
            stats.trajectory_failure_counts[v]++;
        }
    }

    return stats;
}



// calcule les stats d'un vecteur de doubles
inline dataExport::ScalarStats calculate_scalar_stats(const std::vector<double>& data) {
	dataExport::ScalarStats stats;
	if (data.empty()) return stats;

	// min/max
	auto [min_it, max_it] = std::minmax_element(data.begin(), data.end());
	stats.min = *min_it;
	stats.max = *max_it;

	// mean
	double sum = std::accumulate(data.begin(), data.end(), 0.0);
	stats.mean = sum / data.size();

	// variance / stddev
	double sq_sum = std::accumulate(data.begin(), data.end(), 0.0,
		[&](double acc, double v) { return acc + (v - stats.mean) * (v - stats.mean); });
	stats.variance = sq_sum / data.size();
	stats.standard_deviation = std::sqrt(stats.variance);

	// median
	std::vector<double> copy = data;
	std::sort(copy.begin(), copy.end());
	size_t mid = copy.size() / 2;
	stats.median = (copy.size() % 2 == 0) ? 0.5 * (copy[mid - 1] + copy[mid]) : copy[mid];

	return stats;
}

template<genetic::ConfigConcept ConfigType>
void calculate_genes_stats(
	const std::vector<genetic::Individu<ConfigType>>& population,
	GenerationAccumulator& accumulator
) {
	const size_t NB_INDIVIDUS = population.size();
	if (NB_INDIVIDUS == 0) return;

	// On suppose que chaque individu peut se transformer en vecteurs réels
	// la première "gène" est la position initiale, puis impulsions et delta_times
	const auto& genome0 = population[0].to_real_vectors();
	const size_t NB_GENES = genome0.size();
	const size_t NB_IMPULSIONS = (NB_GENES - 1) / 2;

	accumulator.reserve(NB_INDIVIDUS);

	for (const auto& ind : population) {
		const auto& genome = ind.to_real_vectors();

		// 1) Position initiale
		std::array<double, 3> pos = {
			genome[0][0],
			genome[0][1],
			genome[0][2]
		};

		// 2) Impulsions
		std::vector<std::array<double, 3>> imp(NB_IMPULSIONS);
		for (size_t i = 0; i < NB_IMPULSIONS; ++i) {
			const auto& vec = genome[2 * i + 1]; // indices 1,3,5,...
			imp[i] = { vec[0], vec[1], vec[2] };
		}

		// 3) Delta times
		std::vector<double> dt(NB_IMPULSIONS);
		for (size_t i = 0; i < NB_IMPULSIONS; ++i) {
			dt[i] = genome[2 * i + 2][0]; // le premier composant stocke le dt
		}

		// 4) Fitness
		double fitness = ind.get_fitness();

		// push vers l'accumulateur
		accumulator.push_individual(pos, imp, dt, fitness);
	}
}


template<genetic::ConfigConcept ConfigType>
dataExport::TrajectorySoA calculate_trajectory(const genetic::Individu<ConfigType>& individual, SimuCore::Systems::AdaptedSystem system) {
	dataExport::TrajectorySoA traj;
	traj.is_valid = true;

	const size_t MAX_ITERATIONS = system.getMaxIterations();
	double current_time = 0.0;
	const double dt = system.getDeltaTime();
	SimuCore::Systems::AdaptedSystem::RocketState state = SimuCore::Systems::AdaptedSystem::RocketState::NEUTRAL;


	traj.t.reserve(MAX_ITERATIONS);
	traj.px.reserve(MAX_ITERATIONS);
	traj.py.reserve(MAX_ITERATIONS);
	traj.pz.reserve(MAX_ITERATIONS);
	traj.vx.reserve(MAX_ITERATIONS);
	traj.vy.reserve(MAX_ITERATIONS);
	traj.vz.reserve(MAX_ITERATIONS);


	// Transformation de l'individu en fusée
	auto [rocket, gen_state] = SimuCore::IndividualToRocket(individual.to_real_vectors(), system);

	if (gen_state != SimuCore::GenerationState::VALID) {
		traj.is_valid = false;
		return traj;
	}

	// On prépare entities une seule fois
	std::vector<SimuCore::Structures::Entity*> entities;
	entities.push_back(&rocket); 
	entities.push_back(const_cast<SimuCore::Structures::Planet*>(&system.getSun()));
	entities.push_back(const_cast<SimuCore::Structures::Planet*>(&system.getStartPlanet()));
	entities.push_back(const_cast<SimuCore::Structures::Planet*>(&system.getFinalPlanet()));

	for (size_t iter = 0; iter < MAX_ITERATIONS; ++iter) {
		// Intégration
		SimuCore::Integrator::IntegrateStep(entities, dt, current_time);

		// Stockage des positions / vitesses
		traj.t.push_back(current_time);
		traj.px.push_back(entities[0]->position[0]);
		traj.py.push_back(entities[0]->position[1]);
		traj.pz.push_back(entities[0]->position[2]);
		traj.vx.push_back(entities[0]->velocity[0]);
		traj.vy.push_back(entities[0]->velocity[1]);
		traj.vz.push_back(entities[0]->velocity[2]);

		// Vérification de l'état
		state = SimuCore::Systems::GetRocketState(
			*static_cast<SimuCore::Structures::Rocket*>(entities[0]), system);

		if (state != SimuCore::Systems::AdaptedSystem::RocketState::NEUTRAL) {
			break;
		}

		current_time += dt;
	}

	//// redimensionner les vecteurs pour qu'ils correspondent à la taille réelle
	//traj.t.shrink_to_fit();
	//traj.px.shrink_to_fit();
	//traj.py.shrink_to_fit();
	//traj.pz.shrink_to_fit();
	//traj.vx.shrink_to_fit();
	//traj.vy.shrink_to_fit();
	//traj.vz.shrink_to_fit();

	traj.is_valid = state == SimuCore::Systems::AdaptedSystem::RocketState::VALID || state == SimuCore::Systems::AdaptedSystem::RocketState::NEUTRAL;
	return traj;
}


template<genetic::ConfigConcept ConfigType>
void push_population_to_accumulator(
	GenerationAccumulator& accumulator,
	const std::vector<genetic::Individu<ConfigType>>& population,
	const std::vector<dataExport::IndividualValidity>& ind_valid,
	const std::vector<dataExport::TrajectoryValidity>& traj_valid
)
{
	assert(population.size() == ind_valid.size());
	assert(population.size() == traj_valid.size());

	for (size_t idx = 0; idx < population.size(); ++idx) {
		const auto& ind = population[idx];

		// On transforme l'individu en vecteurs réels
		const auto real_vectors = ind.to_real_vectors(); // std::vector<std::vector<double>>

		// Position initiale : on prend le premier vecteur
		std::array<double, 3> init_pos = { real_vectors[0][0], real_vectors[0][1], real_vectors[0][2] };

		// Impulsions et delta_times
		std::vector<std::array<double, 3>> impulses;
		std::vector<double> delta_times;

		impulses.reserve((real_vectors.size() - 1) / 2);
		delta_times.reserve((real_vectors.size() - 1) / 2);

		for (size_t i = 1; i < real_vectors.size(); i += 2) {
			impulses.push_back({ real_vectors[i][0], real_vectors[i][1], real_vectors[i][2] });
			delta_times.push_back(real_vectors[i + 1][0]); // premier composant = dt_i
		}

		// On passe aussi les validités individuelles et trajectoires
		accumulator.push_individual(
			init_pos,
			impulses,
			delta_times,
			ind.get_fitness(),
			ind_valid[idx],
			traj_valid[idx]
		);
	}
}


template<genetic::ConfigConcept ConfigType>
void calculate_population_validity(
	const std::vector<genetic::Individu<ConfigType>>& population,
	const SimuCore::Systems::AdaptedSystem& system,
	std::vector<dataExport::IndividualValidity>& ind_valid,
	std::vector<dataExport::TrajectoryValidity>& traj_valid
) {
	ind_valid.clear();
	traj_valid.clear();
	ind_valid.reserve(population.size());
	traj_valid.reserve(population.size());

	SimuCore::Systems::AdaptedSystem local_system = system;

	for (const auto& ind : population) {
		auto [rocket, gen_state] = SimuCore::IndividualToRocket(ind.to_real_vectors(), local_system);

		if (gen_state != SimuCore::GenerationState::VALID) {
			ind_valid.push_back(dataExport::IndividualValidity::InvalidGenome);
			traj_valid.push_back(dataExport::TrajectoryValidity::EarlyAbort);
			continue;
		}
		else {
			ind_valid.push_back(dataExport::IndividualValidity::Valid);
		}

		auto traj = calculate_trajectory(ind, local_system);
		if (!traj.is_valid) {
			traj_valid.push_back(dataExport::TrajectoryValidity::NumericalDivergence);
		}
		else {
			traj_valid.push_back(dataExport::TrajectoryValidity::Valid);
		}
	}
}
