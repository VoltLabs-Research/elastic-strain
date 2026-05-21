#include <volt/cli/common.h>
#include <volt/elastic_strain_service.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/structures/crystal_topology_registry.h>
#include <oneapi/tbb/global_control.h>
#include <tbb/info.h>

#include <algorithm>
#include <fstream>
#include <set>

using namespace Volt;
using namespace Volt::CLI;

LatticeStructureType parseCrystalStructure(const std::string& str) {
    if (str == "FCC") return LATTICE_FCC;
    if (str == "BCC") return LATTICE_BCC;
    if (str == "HCP") return LATTICE_HCP;
    if (str == "SC")  return LATTICE_SC;
    if (str == "CUBIC_DIAMOND") return LATTICE_CUBIC_DIAMOND;
    if (str == "HEX_DIAMOND")   return LATTICE_HEX_DIAMOND;
    spdlog::warn("Unknown crystal structure '{}', defaulting to BCC.", str);
    return LATTICE_BCC;
}

void showUsage(const std::string& name) {
    printUsageHeader(name, "Volt - Elastic Strain Analysis");
    std::cerr
        << "  --clusters_table <path>       Path to *_clusters.table exported upstream.\n"
        << "  --clusters_transitions <path> Path to *_cluster_transitions.table exported upstream.\n"
        << "  --crystal_structure <type>     Crystal structure. (BCC|FCC|HCP|...) [default: BCC]\n"
        << "  --lattice_dir <path>          Directory containing lattice topology YAMLs.\n"
        << "  --lattice_constant <float>     Lattice constant a₀. [required]\n"
        << "  --ca_ratio <float>             c/a ratio for HCP/hex crystals. [default: 1.0]\n"
        << "  --push_forward                 Push to spatial frame (Euler strain). [default: false]\n"
        << "  --calc_deformation_gradient     Compute deformation gradient F. [default: true]\n"
        << "  --calc_strain_tensors           Compute strain tensors. [default: true]\n"
        << "  --threads <int>               Max worker threads (TBB/OMP). [default: auto]\n";
    printHelpOption();
}

int main(int argc, char* argv[]){
    if(argc < 2){
        showUsage(argv[0]);
        return 1;
    }

    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);

    if(hasOption(opts, "--help") || filename.empty()){
        showUsage(argv[0]);
        return filename.empty() ? 1 : 0;
    }

    if(!hasOption(opts, "--lattice_constant")){
        spdlog::error("--lattice_constant is required for elastic strain analysis.");
        showUsage(argv[0]);
        return 1;
    }

    if(!hasOption(opts, "--threads")) {
        const int maxAvailableThreads = static_cast<int>(oneapi::tbb::info::default_concurrency());
        int physicalCores = 0;
        std::ifstream cpuinfo("/proc/cpuinfo");
        if(cpuinfo.is_open()) {
            std::set<std::pair<int, int>> physicalCoreIds;
            int fallbackCpuCores = 0;
            int physicalId = -1;
            int coreId = -1;
            std::string line;
            while(std::getline(cpuinfo, line)) {
                if(line.empty()) {
                    if(physicalId >= 0 && coreId >= 0) {
                        physicalCoreIds.emplace(physicalId, coreId);
                    }
                    physicalId = -1;
                    coreId = -1;
                    continue;
                }
                if(line.rfind("physical id", 0) == 0) {
                    physicalId = std::stoi(line.substr(line.find(':') + 1));
                } else if(line.rfind("core id", 0) == 0) {
                    coreId = std::stoi(line.substr(line.find(':') + 1));
                } else if(line.rfind("cpu cores", 0) == 0) {
                    fallbackCpuCores = std::max(fallbackCpuCores, std::stoi(line.substr(line.find(':') + 1)));
                }
            }
            if(physicalId >= 0 && coreId >= 0) {
                physicalCoreIds.emplace(physicalId, coreId);
            }
            physicalCores = !physicalCoreIds.empty()
                ? static_cast<int>(physicalCoreIds.size())
                : fallbackCpuCores;
        }
        int defaultThreads = maxAvailableThreads;
        if(physicalCores > 0) {
            defaultThreads = std::min(maxAvailableThreads, physicalCores);
        }
        opts["--threads"] = std::to_string(std::max(1, defaultThreads));
    }

    const int requestedThreads = getInt(opts, "--threads");
    oneapi::tbb::global_control parallelControl(
        oneapi::tbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(std::max(1, requestedThreads))
    );
    initLogging("volt-elastic-strain");
    spdlog::info("Using {} threads (OneTBB)", requestedThreads);

    const std::string latticeDirectory = getString(opts, "--lattice_dir", "");
    if(!latticeDirectory.empty()){
        setCrystalTopologySearchRoot(latticeDirectory);
        spdlog::info("Using lattice directory: {}", latticeDirectory);
    }

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    ElasticStrainService analyzer;
    analyzer.setClustersTablePath(getString(opts, "--clusters_table"));
    analyzer.setClusterTransitionsPath(getString(opts, "--clusters_transitions"));
    analyzer.setInputCrystalStructure(parseCrystalStructure(getString(opts, "--crystal_structure", "BCC")));
    analyzer.setParameters(
        getDouble(opts, "--lattice_constant", 1.63),
        getDouble(opts, "--ca_ratio", 1.0),
        getBool(opts, "--push_forward", false),
        getBool(opts, "--calc_deformation_gradient", true),
        getBool(opts, "--calc_strain_tensors", true)
    );

    spdlog::info("Starting elastic strain analysis...");
    json result = analyzer.compute(frame, outputBase);

    if(result.value("is_failed", false)){
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }
    
    spdlog::info("Elastic strain analysis completed.");
    return 0;
}
