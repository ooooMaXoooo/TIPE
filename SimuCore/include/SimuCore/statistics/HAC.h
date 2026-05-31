#pragma once

#include <pch.h>
#include <DataExport/RocketData.h>

namespace SimuCore {
	namespace Statistics {
		using Cluster = std::vector<int>;
		std::vector<Cluster> HAC(const std::vector<std::unique_ptr<RocketData>>& rockets);
	}; // namespace Statistics
}; // namespace SimuCore