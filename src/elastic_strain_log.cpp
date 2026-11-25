#include <volt/elastic_strain_log.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/core/type_counts.h>
#include <volt/elastic_strain_engine.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <limits>
#include <string>

namespace Volt{

void ElasticStrainLog::dump(const StructureAnalysis& loaded){
    // The cluster graph always holds the null cluster 0.
    const std::size_t clusterCount = loaded.clusterGraph().clusters().empty()
        ? 0 : loaded.clusterGraph().clusters().size() - 1;
    spdlog::info("Structure identification loaded in {}: {} clusters", phaseElapsed(), clusterCount);
}

void ElasticStrainLog::dump(const ElasticStrainEngine& engine){
    spdlog::info("Elastic strain computed in {}", phaseElapsed());
    _results.restart();
    const ParticleProperty* volumetric = engine.volumetricStrains();
    if(!volumetric || volumetric->size() == 0){
        return;
    }
    const std::size_t n = volumetric->size();
    double totalVolumetric = 0.0;
    double minVolumetric = std::numeric_limits<double>::infinity();
    double maxVolumetric = -std::numeric_limits<double>::infinity();
    for(std::size_t i = 0; i < n; i++){
        const double value = volumetric->getDouble(i);
        totalVolumetric += value;
        minVolumetric = std::min(minVolumetric, value);
        maxVolumetric = std::max(maxVolumetric, value);
    }
    spdlog::info(
        "Volumetric strain: min {:.6f}, mean {:.6f}, max {:.6f}",
        minVolumetric, totalVolumetric / n, maxVolumetric
    );
}

void ElasticStrainLog::dump(const ElasticStrainService::AtomStructures& structures){
    const StructureIdentificationExport::AtomsByStructure& byStructure = structures.byStructure;
    TypeCounts line{"Structure types", {}, byStructure.atoms, TypeCounts::Share{
        static_cast<std::size_t>(byStructure.clusteredAtoms), "atoms in clusters"}};
    for(int st : byStructure.order){
        line.counts.emplace_back(st == static_cast<int>(StructureType::OTHER)
            ? "Other" : StructureIdentificationExport::structureTypeNameForExport(st),
            byStructure.atomsOfType[static_cast<std::size_t>(st)].size());
    }
    dump(line);
    if(byStructure.atoms > 0 && byStructure.atomsOfType[static_cast<std::size_t>(structures.referenceStructure)].empty()){
        spdlog::warn(
            "No atom is in a cluster of the {} reference structure, so every elastic strain is zero; check crystal_structure",
            latticeStructureTypeName(structures.referenceStructure)
        );
    }
}

void ElasticStrainLog::dump(const ElasticStrainService::Outputs&){
    spdlog::info("Per-atom results written in {}", _results.elapsed());
}

}
