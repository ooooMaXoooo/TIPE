#include <pch.h>
#include <SimuCore/Optimization/optimization.h>

#include <DataExport/TrajectoryData.h>
#include <DataExport/PhysicsData.h>
#include <DataExport/StatisticsData.h>

namespace SimuCore {
namespace Optimization {
	void writeTrajectory(std::ofstream& filestream, const std::vector<glm::dvec3>& trajectory)
	{
		if (!filestream.is_open()) {
			std::cerr << "\nErreur lors de l'ouverture d'un fichier" << std::endl;
			std::abort();
		}

		for (const glm::dvec3& vec : trajectory) {
			filestream << vec.x << ';' << vec.y << '\n';
		}
	}

	std::string generate_snapshot_directory()
	{
		// Récupérer le temps actuel
		auto now = std::chrono::system_clock::now();
		std::time_t t_now = std::chrono::system_clock::to_time_t(now);
		std::tm local_tm;

#if defined(_WIN32) || defined(_WIN64)
		localtime_s(&local_tm, &t_now);  // Windows
#else
		localtime_r(&t_now, &local_tm);  // Linux / macOS
#endif

		// Construire le nom sous la forme simu_jour_mois_annee_heure_minute_seconde
		std::ostringstream oss;
		oss << "simu_"
			<< std::setfill('0') << std::setw(2) << local_tm.tm_mday << "_"
			<< std::setfill('0') << std::setw(2) << local_tm.tm_mon + 1 << "_"
			<< local_tm.tm_year + 1900 << "_"
			<< std::setfill('0') << std::setw(2) << local_tm.tm_hour << "_"
			<< std::setfill('0') << std::setw(2) << local_tm.tm_min << "_"
			<< std::setfill('0') << std::setw(2) << local_tm.tm_sec;

		return oss.str();
	}

	std::string generate_snapshot_filename(const char* prefix, size_t id, const char* suffix, const std::string& extension)
	{
		std::ostringstream oss;
		oss << prefix << id << suffix << "." << extension;

		return oss.str();
	}

	void sendPlanetsTrajectory(const std::filesystem::path& filepath_start, const std::filesystem::path& filepath_final, SimuCore::Systems::AdaptedSystem system) {
		AsyncDataExporter planetsDataExporter;
		TrajectoryData traj_data_start(system.getStartPlanetPositions(), 2);
		TrajectoryData traj_data_final(system.getFinalPlanetPositions(), 2);

		try
		{
			planetsDataExporter.enqueue(std::make_shared<TrajectoryData>(traj_data_start), filepath_start);
			planetsDataExporter.enqueue(std::make_shared<TrajectoryData>(traj_data_final), filepath_final);
		}
		catch (const std::exception& e)
		{
			std::cerr << "\n\n" << e.what() << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	void sendIndividualTrajectory(const SimuCore::Structures::Rocket& rocket, GenerationState gen_state, const std::filesystem::path& filepath, SimuCore::Systems::AdaptedSystem* system, AsyncDataExporter& exporter) {
		/*
		system->SetRocket(rocket);
		std::vector<glm::dvec3> trajectory;
		system->GetRocketTrajectory(trajectory);
		*/

		std::vector<glm::dvec3> trajectory;
		const size_t trajectory_size = static_cast<const size_t>(daysInSeconds(system->getMaxTime()) / system->getDeltaTime());
		trajectory.reserve(trajectory_size + 1);

		auto position_callback = [&trajectory](const glm::dvec3& pos) {
			trajectory.push_back(pos);
		};

		system->Score(rocket, gen_state, position_callback);

		TrajectoryData trajectoryData{ trajectory, 2 };
		try {
			exporter.enqueue(std::make_shared<TrajectoryData>(trajectoryData), filepath);
		}
		catch (const std::exception& e) {
			std::cerr << "\n\nTrajectory writing error ||\t" << e.what() << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	void sendRocketPhysics(const SimuCore::Structures::Rocket& rocket, const std::filesystem::path& filepath, SimuCore::Systems::AdaptedSystem* system, AsyncDataExporter& exporter) {
		// en jours
		double max_time = system->getMaxTime();
		// en secondes
		double dt = system->getDeltaTime();

		uint8_t startPlanetIndex = SimuCore::Systems::AdaptedSystem::GetStartPlanetID();
		uint8_t finalPlanetIndex = SimuCore::Systems::AdaptedSystem::GetFinalPlanetID();

		glm::dvec3 startPlanet_initialPosition = system->GetStartPlanet_StartPosition();
		glm::dvec3 startPlanet_finalPosition = system->GetStartPlanet_CurrentPosition();
		glm::dvec3 finalPlanet_initialPosition = system->GetFinalPlanet_StartPosition();
		glm::dvec3 finalPlanet_finalPosition = system->GetFinalPlanet_CurrentPosition();

		// en km/s
		glm::dvec3 finalVelocity_rocket = system->GetRocketVelocity();


		std::vector<std::pair<SimuCore::Structures::Impulsion, double>> impulsions = rocket.getImpulsions();

		std::vector<double> impulsions_times;
		std::vector<glm::dvec3> impulsions_vectors;

		impulsions_times.reserve(impulsions.size());
		impulsions_vectors.reserve(impulsions.size());

		for (auto& impuls : impulsions) {
			auto& [vec, instant] = impuls;

			impulsions_times.emplace_back(instant);
			impulsions_vectors.emplace_back(vec.GetDeltaV_vec());
		}

		bool etat_lie = false;

		// en km
		auto [r_min, r_max] = system->GetApsidesAroundFinalPlanet(&etat_lie);


		r_min = meters_to_kilometers(r_min);
		r_max = meters_to_kilometers(r_max);

		constexpr uint8_t dimension = 2;


		double tof = system->getCurrentTime(); // en jours
		double delta_v = rocket.getDeltaV(tof);

		std::cout << "\tTime of flight : " << tof << " (jours)\n";
		std::cout << "\tDelta V : " << delta_v << " (km/s)\n";

		PhysicsData phys_data{ max_time, dt, startPlanetIndex, finalPlanetIndex,
			startPlanet_initialPosition,
			finalPlanet_initialPosition,
			startPlanet_finalPosition,
			finalPlanet_finalPosition,
			finalVelocity_rocket,
			impulsions_times,
			impulsions_vectors,
			tof, delta_v,
			etat_lie, r_min, r_max,
			dimension
		};

		try {
			exporter.enqueue(std::make_shared<PhysicsData>(phys_data), filepath);
		}
		catch (const std::exception& e) {
			std::cerr << "\n\nTrajectory writing error ||\t" << e.what() << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	void sendStatistics(double best_fit, double worst_fit, double mean_score, int kinds[12], const std::vector<SimuCore::Statistics::Cluster>& clusters, const std::filesystem::path& filepath, AsyncDataExporter& exporter) {		
		std::array<int, 12> kindsCount;
		for (int i = 0; i < 12; i++) {
			kindsCount[i] = kinds[i];
		}
		
		StatisticsData stats_data{ clusters.size(), best_fit, worst_fit, mean_score, kindsCount };

		try {
			exporter.enqueue(std::make_shared<StatisticsData>(stats_data), filepath);
		}
		catch (const std::exception& e) {
			std::cerr << "\n\nTrajectory writing error ||\t" << e.what() << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}; // namespace Optimization
}; // namespace SimuCore