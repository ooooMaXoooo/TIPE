#pragma once

#include "HDF5AsyncWriter.h"
#include <DataExport/Writables/GenerationSnapshot.h>

#include <pch.h>
#include <GaLib/genetic_config.hpp>
#include <SimuCore/structures/System.h>


/**
* Cette classe a pour but d'écrire des données générales à propos d'une simulation.
*/

class SimulationAsyncWriter : public HDF5AsyncWriter {
public :

	/// <summary>
	/// 
	/// </summary>
	/// <typeparam name="ConfigType">Un type de configuration spécifique utilisé par la simulation</typeparam>
	/// <param name="directory">l'endroit où vont être écrit les fichiers</param>
	/// <param name="config"></param>
	/// <param name="start_planet_name"></param>
	/// <param name="final_planet_name"></param>
	/// <param name="simulation_duration"></param>
	/// <param name="timestep"></param>
	/// <param name="start_planet_trajectory"></param>
	/// <param name="final_planet_trajectory"></param>
	template <genetic::ConfigConcept ConfigType>
	SimulationAsyncWriter(
		const std::string& directory,
		const ConfigType& config,
		unsigned int seed,
		const SimuCore::Systems::PlanetsName start_planet_name,
		const SimuCore::Systems::PlanetsName final_planet_name,
		double simulation_duration,
		double timestep,
		const std::vector<glm::dvec3>& start_planet_trajectory,
		const std::vector<glm::dvec3>& final_planet_trajectory)
		: 
		HDF5AsyncWriter(directory)
	{
		{ // On demande un lock pour écrire la première partie des données
			std::unique_lock lock(m_mutex);
			H5::H5File simulation_info_file{ (m_directory / GenerateFilename()).c_str(), H5F_ACC_TRUNC};

			// Groupe physique
			{
				H5::Group physics_group = simulation_info_file.createGroup("Physics info");
				WritePhysicsInfo(physics_group,
					start_planet_name, final_planet_name,
					simulation_duration, timestep,
					start_planet_trajectory, final_planet_trajectory);
				physics_group.close();
			}

			// Groupe génétique
			{
				H5::Group genetics_group = simulation_info_file.createGroup("Genetics info");
				WriteGeneticInfo(genetics_group, config, seed);
				genetics_group.close();
			}

			simulation_info_file.close();
		}
	}

	void enqueue(std::shared_ptr<Writable> data, const std::string filename);

private :
	/// <summary>
	/// Generate the filename of the simulation information file.
	/// </summary>
	std::string GenerateFilename();

	void WritePhysicsInfo(H5::Group& physics_group,
		const SimuCore::Systems::PlanetsName start_planet_name,
		const SimuCore::Systems::PlanetsName final_planet_name,
		double simulation_duration,
		double timestep,
		const std::vector<glm::dvec3>& start_planet_trajectory,
		const std::vector<glm::dvec3>& final_planet_trajectory
	);

	template <genetic::ConfigConcept ConfigType>
	void WriteGeneticInfo(H5::Group& genetics_group, const ConfigType& config, unsigned int seed);

	void write_data(H5::H5File& file, const GenerationSnapshot& data);
};





template<genetic::ConfigConcept ConfigType>
inline void SimulationAsyncWriter::WriteGeneticInfo(H5::Group& genetics_group, const ConfigType& config, unsigned int seed) {
	/// Seed
	{
		H5::DataSpace ds(1, std::array<hsize_t, 1>{1}.data());
		genetics_group.createDataSet("Seed", H5::PredType::NATIVE_UINT32, ds).write(&seed, H5::PredType::NATIVE_UINT32);
	}

	/// Numbers representation
	{
		double interval[2] = { config.min_real, config.max_real };
		size_t integer_info[2] = { config.integer_bits, config.integer_max };

		H5::DataSpace ds_reals(1, std::array<hsize_t, 1>{2}.data());
		H5::DataSpace ds_integers(1, std::array<hsize_t, 1>{2}.data());

		genetics_group.createDataSet("Reals", H5::PredType::NATIVE_DOUBLE, ds_reals).write(interval, H5::PredType::NATIVE_DOUBLE);
		genetics_group.createDataSet("Integers", H5::PredType::NATIVE_HSIZE, ds_integers).write(integer_info, H5::PredType::NATIVE_HSIZE);
	}

	/// Genetic algo info
	{
		H5::DataSpace ds_pop_size(1, std::array<hsize_t, 1>{1}.data());
		genetics_group.createDataSet("Population size", H5::PredType::NATIVE_HSIZE, ds_pop_size).write(&config.population_size, H5::PredType::NATIVE_HSIZE);

		H5::DataSpace ds_max_gen(1, std::array<hsize_t, 1>{1}.data());
		genetics_group.createDataSet("Max generation", H5::PredType::NATIVE_HSIZE, ds_max_gen).write(&config.max_generations, H5::PredType::NATIVE_HSIZE);

		H5::DataSpace ds_tournament_size(1, std::array<hsize_t, 1>{1}.data());
		genetics_group.createDataSet("Tournament size", H5::PredType::NATIVE_UINT16, ds_tournament_size).write(&seed, H5::PredType::NATIVE_UINT16);

		uint16_t espace_vectoriel_info[2] = { config.number_of_vectors, config.dimension };
		H5::DataSpace ds_genes(1, std::array<hsize_t, 1>{2}.data());
		genetics_group.createDataSet("Espace vectoriel", H5::PredType::NATIVE_UINT16, ds_genes).write(
			espace_vectoriel_info, H5::PredType::NATIVE_UINT16
		);

		bool booleans[3] = {
			config.enable_elitism,
			config.enable_auto_adaptation,
			config.crossover_method == genetic::CrossoverType::SINGLE_POINT_BIT_LEVEL
		};
		H5::DataSpace ds_booleans(1, std::array<hsize_t, 1>{3}.data());
		genetics_group.createDataSet("Algorithm used", H5::PredType::NATIVE_HBOOL, ds_booleans).write(booleans, H5::PredType::NATIVE_HBOOL);



		if (config.enable_auto_adaptation) {
			double probas[2] = {
				config.initial_mutation_probability,
				config.initial_self_adaptation_probability
			};
			H5::DataSpace ds_proba(1, std::array<hsize_t, 1>{2}.data());
			genetics_group.createDataSet("Mutations probabilities", H5::PredType::NATIVE_DOUBLE, ds_proba).write(
				probas, H5::PredType::NATIVE_DOUBLE
			);
		}
		else {
			double proba = config.get_fixed_proba();
			H5::DataSpace ds_proba(1, std::array<hsize_t, 1>{1}.data());
			genetics_group.createDataSet("Mutations probabilities", H5::PredType::NATIVE_DOUBLE, ds_proba).write(
				&proba, H5::PredType::NATIVE_DOUBLE
			);
		}


		if (!booleans[2]) {
			double proba = config.uniform_crossover_probability;
			H5::DataSpace ds_proba_uniform_crossover(1, std::array<hsize_t, 1>{1}.data());
			genetics_group.createDataSet("Uniform crossover probability", H5::PredType::NATIVE_DOUBLE, ds_proba_uniform_crossover).write(
				&proba, H5::PredType::NATIVE_DOUBLE
			);
		}
	}
}
