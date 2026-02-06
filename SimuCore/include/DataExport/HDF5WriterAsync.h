#pragma once

#include <pch.h>
#include <DataExport/Snapshot.h>
#include <H5Cpp.h>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <memory>
#include <vector>
#include <string>

class HDF5WriterAsync {
public:
    explicit HDF5WriterAsync(const std::string& filename)
        : m_filename(filename), m_stop(false)
    {
        m_worker = std::thread([this] { this->process(); });
    }

    ~HDF5WriterAsync() {
        {
            std::unique_lock lock(m_mutex);
            m_stop = true;
        }
        m_cv.notify_one();
        if (m_worker.joinable())
            m_worker.join();
    }

    void enqueue(std::shared_ptr<dataExport::Snapshot> snapshot) {
        {
            std::unique_lock lock(m_mutex);
            m_queue.push(snapshot);
        }
        m_cv.notify_one();
    }

private:
    std::string m_filename;
    std::queue<std::shared_ptr<dataExport::Snapshot>> m_queue;
    std::mutex m_mutex;
    std::condition_variable m_cv;
    std::atomic<bool> m_stop;
    std::thread m_worker;

    void process() {
        H5::H5File file(m_filename, H5F_ACC_TRUNC);

        while (true) {
            std::shared_ptr<dataExport::Snapshot> snapshot;
            {
                std::unique_lock lock(m_mutex);
                m_cv.wait(lock, [this] { return !m_queue.empty() || m_stop; });

                if (m_stop && m_queue.empty())
                    break;

                snapshot = m_queue.front();
                m_queue.pop();
            }

            if (snapshot) {
                write_snapshot(file, *snapshot);
            }
        }

        file.close();
    }

    void write_snapshot(H5::H5File& file, const dataExport::Snapshot& snapshot) {
        std::string gen_group_name = "generation_" + std::to_string(snapshot.generation);
        H5::Group gen_group = file.createGroup(gen_group_name);

        write_scores(gen_group, snapshot);
        write_initial_position(gen_group, snapshot);
        write_impulsions(gen_group, snapshot);
        write_delta_times(gen_group, snapshot);
        write_validity(gen_group, snapshot);
        write_best_individual(gen_group, snapshot);
    }

    // --- Helper functions ---
    void write_scores(H5::Group& group, const dataExport::Snapshot& snapshot) {
        H5::Group scores_group = group.createGroup("scores");
        double scores[6] = {
            snapshot.stats.best_score,
            snapshot.stats.worst_score,
            snapshot.stats.mean_score,
            snapshot.stats.median_score,
            snapshot.stats.variance_score,
            snapshot.stats.standard_deviation_score
        };
        H5::DataSpace ds(1, std::array<hsize_t, 1>{6}.data());
        scores_group.createDataSet("values", H5::PredType::NATIVE_DOUBLE, ds).write(scores, H5::PredType::NATIVE_DOUBLE);
    }

    void write_initial_position(H5::Group& group, const dataExport::Snapshot& snapshot) {
        H5::Group pos_group = group.createGroup("initial_position");
        double pos[3] = {
            snapshot.stats.initial_position.x.mean,
            snapshot.stats.initial_position.y.mean,
            snapshot.stats.initial_position.z.mean
        };
        H5::DataSpace ds(1, std::array<hsize_t, 1>{3}.data());
        pos_group.createDataSet("mean", H5::PredType::NATIVE_DOUBLE, ds).write(pos, H5::PredType::NATIVE_DOUBLE);

        const auto& hx = snapshot.stats.initial_position.x.histogram;
        const auto& hy = snapshot.stats.initial_position.y.histogram;
        const auto& hz = snapshot.stats.initial_position.z.histogram;
        if (!hx.bins.empty()) {
            H5::DataSpace ds_hist(1, std::array<hsize_t, 1>{hx.bin_count}.data());
            pos_group.createDataSet("x_histogram", H5::PredType::NATIVE_UINT32, ds_hist).write(hx.bins.data(), H5::PredType::NATIVE_UINT32);
        }
        if (!hy.bins.empty()) {
            H5::DataSpace ds_hist(1, std::array<hsize_t, 1>{hy.bin_count}.data());
            pos_group.createDataSet("y_histogram", H5::PredType::NATIVE_UINT32, ds_hist).write(hy.bins.data(), H5::PredType::NATIVE_UINT32);
        }
        if (!hz.bins.empty()) {
            H5::DataSpace ds_hist(1, std::array<hsize_t, 1>{hz.bin_count}.data());
            pos_group.createDataSet("z_histogram", H5::PredType::NATIVE_UINT32, ds_hist).write(hz.bins.data(), H5::PredType::NATIVE_UINT32);
        }
    }

    void write_impulsions(H5::Group& group, const dataExport::Snapshot& snapshot) {
        H5::Group imp_group = group.createGroup("impulsions");
        for (size_t i = 0; i < snapshot.stats.impulsions.size(); ++i) {
            std::string imp_name = "impulsion_" + std::to_string(i);
            H5::Group single_imp = imp_group.createGroup(imp_name);

            double mean[3] = {
                snapshot.stats.impulsions[i].x.mean,
                snapshot.stats.impulsions[i].y.mean,
                snapshot.stats.impulsions[i].z.mean
            };
            H5::DataSpace ds_mean(1, std::array<hsize_t, 1>{3}.data());
            single_imp.createDataSet("mean", H5::PredType::NATIVE_DOUBLE, ds_mean).write(mean, H5::PredType::NATIVE_DOUBLE);

            const auto& hx = snapshot.stats.impulsions[i].x.histogram;
            const auto& hy = snapshot.stats.impulsions[i].y.histogram;
            const auto& hz = snapshot.stats.impulsions[i].z.histogram;
            if (!hx.bins.empty()) {
                H5::DataSpace ds_hist(1, std::array<hsize_t, 1>{hx.bin_count}.data());
                single_imp.createDataSet("x_histogram", H5::PredType::NATIVE_UINT32, ds_hist).write(hx.bins.data(), H5::PredType::NATIVE_UINT32);
            }
            if (!hy.bins.empty()) {
                H5::DataSpace ds_hist(1, std::array<hsize_t, 1>{hy.bin_count}.data());
                single_imp.createDataSet("y_histogram", H5::PredType::NATIVE_UINT32, ds_hist).write(hy.bins.data(), H5::PredType::NATIVE_UINT32);
            }
            if (!hz.bins.empty()) {
                H5::DataSpace ds_hist(1, std::array<hsize_t, 1>{hz.bin_count}.data());
                single_imp.createDataSet("z_histogram", H5::PredType::NATIVE_UINT32, ds_hist).write(hz.bins.data(), H5::PredType::NATIVE_UINT32);
            }
        }
    }

    void write_delta_times(H5::Group& group, const dataExport::Snapshot& snapshot) {
        H5::Group dt_group = group.createGroup("delta_times");
        std::vector<double> dt_means(snapshot.stats.delta_times.size());
        for (size_t i = 0; i < snapshot.stats.delta_times.size(); ++i)
            dt_means[i] = snapshot.stats.delta_times[i].mean;

        if (!dt_means.empty()) {
            H5::DataSpace ds(1, std::array<hsize_t, 1>{dt_means.size()}.data());
            dt_group.createDataSet("mean", H5::PredType::NATIVE_DOUBLE, ds).write(dt_means.data(), H5::PredType::NATIVE_DOUBLE);
        }
    }

    void write_validity(H5::Group& group, const dataExport::Snapshot& snapshot) {
        const auto& v = snapshot.validity;
        if (v.population_size == 0) return;

        H5::Group val_group = group.createGroup("validity");

        // -------------------------
        // Scalaires globaux
        // -------------------------
        auto write_scalar = [&](const char* name, std::size_t value) {
            hsize_t dim = 1;
            H5::DataSpace ds(1, &dim);
            val_group.createDataSet(name, H5::PredType::NATIVE_HSIZE, ds)
                .write(&value, H5::PredType::NATIVE_HSIZE);
            };

        write_scalar("population_size", v.population_size);
        write_scalar("valid_individuals", v.valid_individuals);
        write_scalar("invalid_individuals", v.invalid_individuals);
        write_scalar("valid_trajectories", v.valid_trajectories);
        write_scalar("invalid_trajectories", v.invalid_trajectories);

        // -------------------------
        // Histogrammes d'échecs
        // -------------------------
        auto write_enum_histogram = [&](const char* name, const auto& map) {
            if (map.empty()) return;

            std::vector<uint8_t> keys;
            std::vector<std::size_t> counts;

            for (const auto& [k, c] : map) {
                keys.push_back(static_cast<uint8_t>(k));
                counts.push_back(c);
            }

            hsize_t dim = keys.size();
            H5::DataSpace ds(1, &dim);

            H5::Group g = val_group.createGroup(name);
            g.createDataSet("codes", H5::PredType::NATIVE_UINT8, ds)
                .write(keys.data(), H5::PredType::NATIVE_UINT8);
            g.createDataSet("counts", H5::PredType::NATIVE_HSIZE, ds)
                .write(counts.data(), H5::PredType::NATIVE_HSIZE);
            };

        write_enum_histogram("individual_failures", v.individual_failure_counts);
        write_enum_histogram("trajectory_failures", v.trajectory_failure_counts);
    }

    void write_best_individual(H5::Group& group, const dataExport::Snapshot& snapshot) {
        if (!snapshot.best_update.is_new_best) return;

        H5::Group best_group = group.createGroup("best_individual");
        const auto& traj = snapshot.best_update.trajectory;

        auto create_vector_dataset = [&](const std::string& name, const std::vector<double>& data) {
            if (data.empty()) return;
            H5::DataSpace ds(1, std::array<hsize_t, 1>{data.size()}.data());
            best_group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds).write(data.data(), H5::PredType::NATIVE_DOUBLE);
            };

        create_vector_dataset("t", traj.t);
        create_vector_dataset("px", traj.px);
        create_vector_dataset("py", traj.py);
        create_vector_dataset("pz", traj.pz);
        create_vector_dataset("vx", traj.vx);
        create_vector_dataset("vy", traj.vy);
        create_vector_dataset("vz", traj.vz);
    }
};
