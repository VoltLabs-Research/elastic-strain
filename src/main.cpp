#include <volt/cli/common.h>
#include <volt/elastic_strain_service.h>
#include <volt/plugin/option_reader.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/structures/crystal_topology_registry.h>
#include <oneapi/tbb/global_control.h>
#include <tbb/info.h>

#include <algorithm>
#include <fstream>
#include <set>

using namespace Volt;
using namespace Volt::CLI;
using namespace Volt::Plugin;

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

static PluginDescriptor buildDescriptor() {
    return {
        "volt-elastic-strain",
        "Elastic Strain Analysis",
        {
            {"--clusters_table", "path", "Clusters table exported by an upstream structure-identification step.", "", {}, ""},
            {"--clusters_transitions", "path", "Cluster transitions table exported upstream.", "", {}, ""},
            {"--neighbor_lattice", "path", "Per-atom neighbor topology parquet exported upstream.", "", {}, ""},
            {"--crystal_structure", "enum", "Crystal structure to match against.", "BCC",
             {"BCC", "FCC", "HCP", "SC", "CUBIC_DIAMOND", "HEX_DIAMOND"}, ""},
            {"--lattice_dir", "path", "Directory containing lattice topology YAMLs.", "", {},
             "share/volt/lattices"},
            {"--lattice_constant", "float", "Lattice constant a0. Required.", "", {}, ""},
            {"--ca_ratio", "float", "c/a ratio for HCP/hex crystals.", "1.0", {}, ""},
            {"--push_forward", "bool", "Push to spatial frame (Euler strain).", "false", {}, ""},
            {"--calc_deformation_gradient", "bool", "Compute deformation gradient F.", "true", {}, ""},
            {"--calc_strain_tensors", "bool", "Compute strain tensors.", "true", {}, ""},
        }
    };
}

int main(int argc, char* argv[]){
    const PluginDescriptor descriptor = buildDescriptor();

    if(argc < 2){
        showPluginUsage(argv[0], descriptor);
        return 1;
    }

    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);

    if(auto exitCode = handleIntrospection(argv[0], descriptor, opts, filename)){
        return *exitCode;
    }

    if(!hasOption(opts, "--lattice_constant")){
        spdlog::error("--lattice_constant is required for elastic strain analysis.");
        showPluginUsage(argv[0], descriptor);
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

    const OptionReader options(descriptor, opts);

    const std::string latticeDirectory = options.text("--lattice_dir");
    if(!latticeDirectory.empty()){
        setCrystalTopologySearchRoot(latticeDirectory);
        spdlog::info("Using lattice directory: {}", latticeDirectory);
    }

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    ElasticStrainService analyzer;
    analyzer.setClustersTablePath(options.text("--clusters_table"));
    analyzer.setClusterTransitionsPath(options.text("--clusters_transitions"));
    analyzer.setNeighborLatticePath(options.text("--neighbor_lattice"));
    analyzer.setInputCrystalStructure(parseCrystalStructure(options.text("--crystal_structure")));
    analyzer.setParameters(
        options.number("--lattice_constant"),
        options.number("--ca_ratio"),
        options.boolean("--push_forward"),
        options.boolean("--calc_deformation_gradient"),
        options.boolean("--calc_strain_tensors")
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
