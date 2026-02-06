#pragma once
#include <pch.h>

namespace dataExport {

    enum class IndividualValidity : uint8_t {
        Valid,
        InvalidGenome,
        InvalidRocketConfig,
        ConstraintViolation
    };


    enum class TrajectoryValidity : uint8_t {
        Valid,
        NumericalDivergence,
        SimulationFailure,
        EarlyAbort
    }; 


    struct GenerationValidityStats {
        std::size_t population_size;

        std::size_t valid_individuals;
        std::size_t invalid_individuals;

        std::size_t valid_trajectories;
        std::size_t invalid_trajectories;

        // optionnel mais puissant
        std::map<IndividualValidity, std::size_t> individual_failure_counts;
        std::map<TrajectoryValidity, std::size_t> trajectory_failure_counts;
    }; // struct GenerationValidityStats



    struct Histogram {
        double min;
        double max;
        uint32_t bin_count;
        std::vector<uint32_t> bins; // taille = bin_count
    };

    struct ScalarStats {
        double mean;
        double variance;
        double standard_deviation;
        double median;
        double min;
        double max;

        Histogram histogram;
    };

    struct VectorStats {
        ScalarStats x;
        ScalarStats y;
        ScalarStats z;
    };

	struct GenerationStats {
		VectorStats initial_position;
		std::vector<VectorStats> impulsions;
		std::vector<ScalarStats> delta_times; // les écarts de temps entre deux impulsions

		double best_score;
		double worst_score;
		double mean_score;
		double median_score;
        double variance_score;
        double standard_deviation_score;
	}; // struct GenerationStats

}; // namespace dataExport