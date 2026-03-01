#include <pch.h>
#include <DataExport/Writers/SimulationAsyncWriter.h>

std::string SimulationAsyncWriter::GenerateFilename()
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

	// Construire le nom sous la forme simulation_jour_mois_annee_heure_minute_seconde.h5
	std::ostringstream oss;
	oss << "simulation_"
		<< std::setfill('0') << std::setw(2) << local_tm.tm_mday << "_"
		<< std::setfill('0') << std::setw(2) << local_tm.tm_mon + 1 << "_"
		<< local_tm.tm_year + 1900 << "_"
		<< std::setfill('0') << std::setw(2) << local_tm.tm_hour << "_"
		<< std::setfill('0') << std::setw(2) << local_tm.tm_min << "_"
		<< std::setfill('0') << std::setw(2) << local_tm.tm_sec
		<< "." << "h5";

	return oss.str();
}

void SimulationAsyncWriter::WritePhysicsInfo(H5::Group& physics_group,
	const SimuCore::Systems::PlanetsName start_planet_name,
	const SimuCore::Systems::PlanetsName final_planet_name,
	double simulation_duration,
	double timestep,
	const std::vector<glm::dvec3>& start_planet_trajectory,
	const std::vector<glm::dvec3>& final_planet_trajectory
	) {

	/// Pas de temps et durée de simulation
	{
		double times[2] = {simulation_duration, timestep};
		H5::DataSpace ds(1, std::array<hsize_t, 1>{2}.data());
		physics_group.createDataSet("Times", H5::PredType::NATIVE_DOUBLE, ds).write(times, H5::PredType::NATIVE_DOUBLE);
	}

	/// Start and final planet's name
	{
		uint8_t id_planets[2] = { static_cast<uint8_t>(start_planet_name), static_cast<uint8_t>(final_planet_name) };
		H5::DataSpace ds(1, std::array<hsize_t, 1>{2}.data());
		physics_group.createDataSet("Planet's name", H5::PredType::NATIVE_UINT8, ds).write(id_planets, H5::PredType::NATIVE_UINT8);
	}

	/// Start and final planet's trajectory
	{
		auto create_vector_dataset = [&](const std::string& name, const std::vector<double>& data) {
			if (data.empty()) return;
			H5::DataSpace ds(1, std::array<hsize_t, 1>{data.size()}.data());
			physics_group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds).write(data.data(), H5::PredType::NATIVE_DOUBLE);
			};

		auto create_trajectory_dataset = [create_vector_dataset](const std::string& base_name, const std::vector<glm::dvec3>&trajectory) {
			std::vector<double> x_coord, y_coord, z_coord;
			x_coord.reserve(trajectory.size());
			y_coord.reserve(trajectory.size());
			z_coord.reserve(trajectory.size());

			for (const glm::dvec3& position : trajectory) {
				x_coord.emplace_back(position.x);
				y_coord.emplace_back(position.y);
				z_coord.emplace_back(position.z);
			}

			create_vector_dataset(base_name + " x", x_coord);
			create_vector_dataset(base_name + " y", y_coord);
			create_vector_dataset(base_name + " z", z_coord);
		};
	}
}
