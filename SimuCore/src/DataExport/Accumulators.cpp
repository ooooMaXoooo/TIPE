#include <pch.h>
#include <DataExport\Accumulators.h>


void ScalarAccumulator::finalize(int bin_count) {
    if (m_count == 0) {
        m_stats = {};
        m_hist = {};
        return;
    }

    // stats scalaires
    m_stats.mean = m_mean;
    m_stats.variance = (m_count > 1) ? (m_M2 / m_count) : 0.0;
    m_stats.standard_deviation = std::sqrt(m_stats.variance);
    m_stats.min = m_min;
    m_stats.max = m_max;

    // médiane
    std::vector<double> sorted = m_values;
    std::sort(sorted.begin(), sorted.end());
    size_t mid = sorted.size() / 2;
    m_stats.median = (sorted.size() % 2 == 0)
        ? 0.5 * (sorted[mid - 1] + sorted[mid])
        : sorted[mid];

    // histogramme
    if (bin_count <= 0) {
        bin_count = static_cast<int>(std::sqrt(m_values.size()));
        bin_count = std::max(bin_count, 1);
    }

    m_hist.min = m_min;
    m_hist.max = m_max;
    m_hist.bin_count = static_cast<uint32_t>(bin_count);
    m_hist.bins.assign(bin_count, 0);

    if (m_min < m_max) {
        double inv_width = bin_count / (m_max - m_min);
        for (double v : m_values) {
            size_t idx = std::min(
                static_cast<size_t>((v - m_min) * inv_width),
                static_cast<size_t>(bin_count - 1)
            );
            m_hist.bins[idx]++;
        }
    }

    m_stats.histogram = m_hist;
}



dataExport::Snapshot GenerationAccumulator::finalize(
    size_t gen,
    const dataExport::BestIndividualUpdate& best_update
) {
    constexpr int HIST_BINS = 32;

    // -------------------------
    // 1) Finalisation des stats
    // -------------------------
    initial_position.finalize(HIST_BINS);

    for (auto& acc : impulses)
        acc.finalize(HIST_BINS);

    for (auto& acc : delta_times)
        acc.finalize(HIST_BINS);

    scores.finalize(HIST_BINS);

    // -------------------------
    // 2) Construction du snapshot
    // -------------------------
    dataExport::Snapshot snapshot;
    snapshot.generation = gen;
    snapshot.best_update = best_update;

    snapshot.stats.initial_position = initial_position.stats();

    snapshot.stats.impulsions.resize(impulses.size());
    for (size_t i = 0; i < impulses.size(); ++i)
        snapshot.stats.impulsions[i] = impulses[i].stats();

    snapshot.stats.delta_times.resize(delta_times.size());
    for (size_t i = 0; i < delta_times.size(); ++i)
        snapshot.stats.delta_times[i] = delta_times[i].stats();

    const auto s = scores.stats();
    snapshot.stats.best_score = s.max;
    snapshot.stats.worst_score = s.min;
    snapshot.stats.mean_score = s.mean;
    snapshot.stats.median_score = s.median;
    snapshot.stats.variance_score = s.variance;
    snapshot.stats.standard_deviation_score = s.standard_deviation;

    snapshot.validity = calculate_validity_stats(
        individual_validity,
        trajectory_validity
    );

    // -------------------------
    // 3) Reset pour la génération suivante
    // -------------------------
    initial_position = VectorAccumulator();
    impulses.clear();
    delta_times.clear();
    scores = ScalarAccumulator();
    individual_validity.clear();
    trajectory_validity.clear();

    return snapshot;
}
