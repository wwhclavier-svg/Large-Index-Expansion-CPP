#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <map>
#include <ctime>
#include <unistd.h>

#include "IBPMatrixLoader_Binary.hpp"
#include "LayerRecursion.hpp"
#include "RingDataLoader.hpp"
#include "RelationSolver.hpp"
#include "Combinatorics.hpp"
#include "LinearSolver.hpp"
#include "firefly/FFInt.hpp"
#include "SeriesCoefficientIO.hpp"
#include "json.hpp"

using namespace std;
using namespace firefly;
using json = nlohmann::json;

static string getHostname() {
    char buf[256] = {0};
    if (gethostname(buf, sizeof(buf)) == 0) return string(buf);
    return "unknown";
}

static string getTimestamp() {
    time_t now = time(nullptr);
    struct tm* t = localtime(&now);
    char buf[64];
    strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", t);
    return string(buf);
}

struct LevelResult {
    int lev, deg;
    int num_variables, active_variables;
    int solution_dimension;
    double time_s;
    int system_rows, system_cols;
    int sampling_points_used;
    bool converged;
    int stable_order;
    int num_relations;
    int independent_pairs;
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <family> [order] [lev_min] [lev_max] [deg_max]" << endl;
        cerr << "  [--output <dir>]   output directory (default: results/lie)" << endl;
        cerr << "  [--topsector]      load only top sector" << endl;
        cerr << "  [--mode <n>]       ansatz mode (0=Pyramid default)" << endl;
        cerr << "  [--sector <s>]     ExtendedPyramid sector" << endl;
        return 1;
    }

    string family = argv[1];
    int order = (argc > 2) ? stoi(argv[2]) : 4;
    int lev_min = (argc > 3) ? stoi(argv[3]) : 1;
    int lev_max = (argc > 4) ? stoi(argv[4]) : 2;
    int deg_max = (argc > 5) ? stoi(argv[5]) : 2;
    int deg_min = 0;
    bool top_sector_only = false;
    RelationSolver::AnsatzMode ansatz_mode = RelationSolver::AnsatzMode::Pyramid;
    vector<int> ext_sector;
    string outDir = "results/lie";
    for (int i = 1; i < argc; ++i) {
        if (string(argv[i]) == "--output" && i+1 < argc) outDir = argv[++i];
        if (string(argv[i]) == "--topsector") top_sector_only = true;
        else if (string(argv[i]) == "--mode" && i+1 < argc) {
            int m = stoi(argv[++i]);
            if (m == 0) ansatz_mode = RelationSolver::AnsatzMode::Pyramid;
            else if (m == 1) ansatz_mode = RelationSolver::AnsatzMode::DotPyramid;
            else if (m == 2) ansatz_mode = RelationSolver::AnsatzMode::Star;
            else if (m == 3) ansatz_mode = RelationSolver::AnsatzMode::ExtendedPyramid;
        }
        else if (string(argv[i]) == "--sector" && i+1 < argc) {
            string s = argv[++i];
            ext_sector.clear();
            for (char c : s) ext_sector.push_back(c == '1' ? 1 : 0);
        }
    }
    const int incre = 2;

    json result_json;
    result_json["family"] = family;
    result_json["timestamp"] = getTimestamp();
    result_json["machine"] = getHostname();
    result_json["config"]["order"] = order;
    result_json["config"]["incre"] = incre;
    result_json["config"]["lev_min"] = lev_min;
    result_json["config"]["lev_max"] = lev_max;
    result_json["config"]["deg_min"] = deg_min;
    result_json["config"]["deg_max"] = deg_max;
    result_json["config"]["ansatz_mode"] = static_cast<int>(ansatz_mode);
    result_json["config"]["top_sector_only"] = top_sector_only;

    initBinomial();

    string binIBPFile = "data/IBPMat_" + family + ".bin";
    string binRingFile = "data/RingData_" + family + ".bin";
    string coeffCacheFile = "data/ExpansionCache_" + family + "_k" + to_string(order) + ".bin";

    auto t_total_start = chrono::high_resolution_clock::now();

    try {
        // ---- Load IBP matrices ----
        auto t0 = chrono::high_resolution_clock::now();
        set<size_t> keep_indices;
        const set<size_t>* keep_ptr = nullptr;
        vector<vector<seriesCoefficient<FFInt>>> allResults;
        vector<AlgebraData::RingCell<FFInt>> ringData;
        int ne = 0, nb = 0;

        if (top_sector_only) {
            auto allRingData = AlgebraData::LoadBinary<FFInt>(binRingFile, FFInt::p);
            int max_sum = -1;
            for (size_t i = 0; i < allRingData.size(); ++i) {
                int sum = 0;
                for (int v : allRingData[i].limitSector) sum += v;
                if (sum > max_sum) { max_sum = sum; keep_indices.clear(); keep_indices.insert(i); }
                else if (sum == max_sum) keep_indices.insert(i);
            }
            keep_ptr = &keep_indices;
            auto ibpmatlist = loadAllIBPMatricesBinary<FFInt>(binIBPFile, false, keep_ptr);
            if (ibpmatlist.empty()) { cerr << "Error: No IBP matrices." << endl; return 1; }
            ne = ibpmatlist[0].ne; nb = ibpmatlist[0].nb;
            ringData.clear();
            for (size_t idx : keep_indices) ringData.push_back(move(allRingData[idx]));
            auto t_expand = chrono::high_resolution_clock::now();
            allResults = batchProcessRecursion<FFInt>(ibpmatlist, order, incre);
            auto t1 = chrono::high_resolution_clock::now();
            result_json["phases"]["load_IBP"]["time_s"] = chrono::duration<double>(t_expand - t0).count();
            result_json["phases"]["expand"]["time_s"] = chrono::duration<double>(t1 - t_expand).count();
            result_json["phases"]["expand"]["cache_hit"] = false;
        } else {
            auto ibpmatlist = loadAllIBPMatricesBinary<FFInt>(binIBPFile);
            if (ibpmatlist.empty()) { cerr << "Error: No IBP matrices." << endl; return 1; }
            ne = ibpmatlist[0].ne; nb = ibpmatlist[0].nb;
            auto t_load = chrono::high_resolution_clock::now();
            result_json["phases"]["load_IBP"]["time_s"] = chrono::duration<double>(t_load - t0).count();
            result_json["phases"]["load_IBP"]["num_regimes"] = ibpmatlist.size();

            // Try cache, else compute
            bool cache_hit = false;
            try {
                allResults = SeriesIO::loadAllResults<FFInt>(coeffCacheFile);
                cache_hit = true;
            } catch (...) {
                auto t_expand = chrono::high_resolution_clock::now();
                allResults = batchProcessRecursion<FFInt>(ibpmatlist, order, incre);
                auto t2 = chrono::high_resolution_clock::now();
                result_json["phases"]["expand"]["time_s"] = chrono::duration<double>(t2 - t_expand).count();
                SeriesIO::saveAllResults(allResults, coeffCacheFile);
            }
            result_json["phases"]["expand"]["cache_hit"] = cache_hit;
            result_json["phases"]["expand"]["k_max"] = order;

            ringData = AlgebraData::LoadBinary<FFInt>(binRingFile, FFInt::p);
        }

        result_json["config"]["ne"] = ne;
        result_json["config"]["nb"] = nb;
        result_json["config"]["prime"] = static_cast<int>(FFInt::p);

        // ---- Extract ring data ----
        vector<vector<int>> sector_list;
        vector<vector<vector<FFInt>>> A_list, Ainv_list;
        for (const auto& ring : ringData) {
            sector_list.push_back(ring.limitSector);
            A_list.push_back(ring.A_list);
            Ainv_list.push_back(ring.Ainv_list);
        }

        // ---- Relation solving ----
        auto t_solve_start = chrono::high_resolution_clock::now();
        RelationSolver::AdaptiveSamplingConfig base_config;
        base_config.min_nu = 3;
        base_config.max_nu = 50;
        base_config.nullity_stable_threshold = 5;
        base_config.check_interval = 5;
        base_config.verification_points = 5;
        base_config.plateau_size = 1;
        base_config.random_min = 0;
        base_config.random_max = static_cast<double>(FFInt::p - 1);

        auto all_results = RelationSolver::reconstructAllRelations<FFInt>(
            allResults, sector_list, A_list, Ainv_list,
            ne, lev_min, lev_max, deg_max, ansatz_mode, base_config, ext_sector
        );

        // ---- Collect results ----
        vector<LevelResult> levels;
        json counts;
        json stability;
        for (const auto& ld : all_results) {
            LevelResult lr;
            lr.lev = ld.lev;
            lr.deg = ld.deg;
            lr.num_variables = ld.total_vars;
            lr.active_variables = ld.active_vars;
            lr.solution_dimension = ld.coeff.getNumSolutions();
            lr.time_s = 0;
            lr.system_rows = 0;
            lr.system_cols = ld.active_vars;
            lr.sampling_points_used = 0;
            lr.converged = ld.stable_order > 0;
            lr.stable_order = ld.stable_order;
            lr.num_relations = ld.num_relations;
            lr.independent_pairs = ld.independent_pairs.size();
            levels.push_back(lr);

            string lev_key = "lev" + to_string(lr.lev);
            string deg_key = "deg" + to_string(lr.deg);
            counts[lev_key][deg_key] = lr.solution_dimension;
        }

        auto t_solve_end = chrono::high_resolution_clock::now();
        double solve_time = chrono::duration<double>(t_solve_end - t_solve_start).count();
        auto t_total_end = chrono::high_resolution_clock::now();
        double total_time = chrono::duration<double>(t_total_end - t_total_start).count();

        json levels_json = json::array();
        for (const auto& lr : levels) {
            levels_json.push_back({
                {"lev", lr.lev},
                {"deg", lr.deg},
                {"num_variables", lr.num_variables},
                {"active_variables", lr.active_variables},
                {"solution_dimension", lr.solution_dimension},
                {"system_rows", lr.system_rows},
                {"system_cols", lr.system_cols},
                {"sampling_points_used", lr.sampling_points_used},
                {"converged", lr.converged},
                {"stable_order", lr.stable_order},
                {"time_s", lr.time_s}
            });
        }

        result_json["phases"]["relation_solve"]["total_time_s"] = solve_time;
        result_json["phases"]["relation_solve"]["levels"] = levels_json;
        result_json["total_time_s"] = total_time;

        // ---- Relation counts table ----
        result_json["relation_counts"]["family"] = family;
        result_json["relation_counts"]["ne"] = ne;
        result_json["relation_counts"]["counts"] = counts;

        // ---- Completeness check (stub) ----
        result_json["completeness"]["status"] = "not_checked";

    } catch (const exception& e) {
        result_json["error"] = e.what();
        result_json["total_time_s"] = chrono::duration<double>(
            chrono::high_resolution_clock::now() - t_total_start).count();
    }

    // ---- Write output JSON ----
    if (outDir.back() == '/') outDir.pop_back();
    (void)system(("mkdir -p " + outDir).c_str());
    string outFile = outDir + "/" + family + ".json";
    ofstream out(outFile);
    out << result_json.dump(2) << endl;
    out.close();
    cout << "Benchmark results written to: " << outFile << endl;
    cout << "Total time: " << result_json["total_time_s"] << " s" << endl;

    return 0;
}
