#pragma once
#include <DataExport/GenerationStats.h>
#include <DataExport/TrajectorySoA.h>

namespace dataExport {

    struct Snapshot {
        size_t generation = 0;

        // Statistiques de la génération, y compris validité
        GenerationStats stats;

        // Informations sur le meilleur individu
        BestIndividualUpdate best_update;

        // Validité globale de la génération
        GenerationValidityStats validity;
	}; // struct Snapshot

}; // namespace dataExport