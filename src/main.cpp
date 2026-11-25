#include <volt/plugin/plugin_entry.h>
#include <volt/elastic_strain_service.h>
#include <volt/structures/crystal_topology_registry.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = ElasticStrainService;

static const std::vector<OptionBinding<S>> bindings = {
    opt({"--clusters_table", "path", "Clusters table exported by an upstream structure-identification step.", "", {}, ""},
        std::function<void(S&, const OptsMap&)>([](S& s, const OptsMap& opts) {
            s.setClustersTablePath(CLI::getString(opts, "--clusters_table", ""));
        })),
    opt({"--clusters_transitions", "path", "Cluster transitions table exported upstream.", "", {}, ""},
        std::function<void(S&, const OptsMap&)>([](S& s, const OptsMap& opts) {
            s.setClusterTransitionsPath(CLI::getString(opts, "--clusters_transitions", ""));
        })),
    opt({"--neighbor_lattice", "path", "Per-atom neighbor topology parquet exported upstream.", "", {}, ""},
        std::function<void(S&, const OptsMap&)>([](S& s, const OptsMap& opts) {
            s.setNeighborLatticePath(CLI::getString(opts, "--neighbor_lattice", ""));
        })),
    opt({"--crystal_structure", "enum", "Crystal structure to match against.", "BCC",
         {"BCC", "FCC", "HCP", "SC", "CUBIC_DIAMOND", "HEX_DIAMOND"}, ""},
        std::function<void(S&, const OptsMap&)>([](S& s, const OptsMap& opts) {
            const std::string name = CLI::getString(opts, "--crystal_structure", "BCC");
            LatticeStructureType structure = LATTICE_BCC;
            if (name == "FCC") structure = LATTICE_FCC;
            else if (name == "HCP") structure = LATTICE_HCP;
            else if (name == "SC") structure = LATTICE_SC;
            else if (name == "CUBIC_DIAMOND") structure = LATTICE_CUBIC_DIAMOND;
            else if (name == "HEX_DIAMOND") structure = LATTICE_HEX_DIAMOND;
            else if (name != "BCC") spdlog::warn("Unknown crystal structure '{}', defaulting to BCC.", name);
            s.setInputCrystalStructure(structure);
        })),
    opt({"--lattice_dir", "path", "Directory containing lattice topology YAMLs.", "", {}, "share/volt/lattices"},
        std::function<void(S&, const OptsMap&)>([](S&, const OptsMap& opts) {
            const std::string directory = CLI::getString(opts, "--lattice_dir", "");
            if (directory.empty()) return;
            setCrystalTopologySearchRoot(directory);
            spdlog::info("Using lattice directory: {}", directory);
        })),
    opt({"--lattice_constant", "float", "Lattice constant a0. Required.", "", {}, ""},
        std::function<void(S&, const OptsMap&)>([](S& s, const OptsMap& opts) {
            s.setLatticeConstant(CLI::getDouble(opts, "--lattice_constant", 0.0));
        })),
    opt("--ca_ratio", "c/a ratio for HCP/hex crystals.", 1.0, &S::setCaRatio),
    opt("--push_forward", "Push to spatial frame (Euler strain).", false, &S::setPushForward),
    opt("--calc_deformation_gradient", "Compute deformation gradient F.", true, &S::setCalculateDeformationGradient),
    opt("--calc_strain_tensors", "Compute strain tensors.", true, &S::setCalculateStrainTensors),
};

VOLT_PLUGIN_MAIN((PluginDescriptor{
        .name = "volt-elastic-strain",
        .description = "Elastic Strain Analysis",
        .options = optionsMeta(bindings),
        .defaultThreads = DefaultThreads::PhysicalCores,
    }),
    serviceRun(bindings, {"--lattice_constant"}))
