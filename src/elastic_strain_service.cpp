#include <volt/elastic_strain_service.h>
#include <volt/elastic_strain_engine.h>
#include <volt/elastic_strain_log.h>
#include <volt/core/reconstructed_structure.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/progress.h>
#include <volt/utilities/json_utils.h>
#include <volt/utilities/written_table.h>
#include <utility>

namespace Volt{

using namespace Volt::Particles;

const OutputFile<ElasticStrainService, ElasticStrainService::Outputs> ElasticStrainService::outputFiles[] = {
    {"Writing the elastic strain table", 3.0, nullptr, [](const Outputs& out){
        out.log.dump(JsonUtils::writeSummaryTable(out.result, out.base, "elastic_strain", "Elastic strain"));
    }},
    {"Writing the atoms table", 6.0, nullptr, [](const Outputs& out){
        const AtomStructures structures = atomStructures(out);
        out.log.dump(structures);
        const AtomsTable atoms{outputPath(out.base, "atoms")};
        JsonUtils::writeJsonToParquet(atomsTable(out, structures), atoms.path, false);
        out.log.dump(atoms);
    }},
};

std::vector<int> ElasticStrainService::atomStructureTypes(const StructureAnalysis& analysis, std::size_t atomCount){
    std::vector<int> structureTypes(atomCount, static_cast<int>(StructureType::OTHER));
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        Cluster* cluster = analysis.atomCluster(static_cast<int>(atomIndex));
        if(cluster && 0 <= cluster->structure && cluster->structure < static_cast<int>(StructureType::NUM_STRUCTURE_TYPES)){
            structureTypes[atomIndex] = cluster->structure;
        }
    }
    return structureTypes;
}

ElasticStrainService::ElasticStrainService()
    : _inputCrystalStructure(LATTICE_BCC),
      _latticeConstant(1.63f),
      _caRatio(1.0),
      _pushForward(false),
      _calculateDeformationGradient(true),
      _calculateStrainTensors(true){}

void ElasticStrainService::setInputCrystalStructure(LatticeStructureType structure){
    _inputCrystalStructure = structure;
}

void ElasticStrainService::setClustersTablePath(std::string path){
    _clustersTablePath = std::move(path);
}

void ElasticStrainService::setClusterTransitionsPath(std::string path){
    _clusterTransitionsPath = std::move(path);
}

void ElasticStrainService::setNeighborLatticePath(std::string path){
    _neighborLatticePath = std::move(path);
}

void ElasticStrainService::setLatticeConstant(double latticeConstant){
    _latticeConstant = latticeConstant;
}

void ElasticStrainService::setCaRatio(double caRatio){
    _caRatio = caRatio;
}

void ElasticStrainService::setPushForward(bool pushForward){
    _pushForward = pushForward;
}

void ElasticStrainService::setCalculateDeformationGradient(bool calculateDeformationGradient){
    _calculateDeformationGradient = calculateDeformationGradient;
}

void ElasticStrainService::setCalculateStrainTensors(bool calculateStrainTensors){
    _calculateStrainTensors = calculateStrainTensors;
}

json ElasticStrainService::compute(const LammpsParser::Frame &frame, const std::string &outputFilename){
    if(_clustersTablePath.empty() || _clusterTransitionsPath.empty()){
        return AnalysisResult::failure(
            "ElasticStrain requires --clusters_table and --clusters_transitions"
        );
    }
    if(_neighborLatticePath.empty()){
        return AnalysisResult::failure(
            "ElasticStrain requires --neighbor_lattice (per-atom neighbor topology parquet)"
        );
    }

    ElasticStrainLog log;
    FrameAdapter::PreparedAnalysisInput prepared;
    std::string frameError;
    if(!FrameAdapter::prepareAnalysisInput(frame, prepared, &frameError))
        return AnalysisResult::failure(frameError);

    auto positions = std::move(prepared.positions);

    ReconstructedStructureContext context(positions.get(), frame.simulationCell);
    context.inputCrystalType = _inputCrystalStructure;

    StructureAnalysis analysis(context);
    // Weighed by a 51k-atom FCC run (157 ms, 2 ms, 0.9 s, 0.47 s, 0.91 s): the per-atom JSON
    // and its exports dominate.
    ProgressSteps steps({
        {"Loading the upstream structure identification", 1.0},
        {"Computing elastic strain", 0.05},
        {"Building the per-atom results", 6.0},
        {"Writing outputs", totalWeight(outputSteps(outputFiles, *this))},
    });
    std::string reconstructionError;
    if(!ReconstructedStructureLoader::load(
        frame,
        _neighborLatticePath,
        {_clustersTablePath, _clusterTransitionsPath},
        analysis,
        context,
        &reconstructionError
    )){
        return AnalysisResult::failure(reconstructionError);
    }
    rebuildImportedClusterParentHierarchy(analysis);
    log.dump(analysis);

    steps.next();
    ElasticStrainEngine engine(
        analysis,
        context,
        static_cast<LatticeStructureType>(_inputCrystalStructure),
        _calculateDeformationGradient,
        _calculateStrainTensors,
        _latticeConstant,
        _caRatio,
        _pushForward
    );
    engine.perform();
    log.dump(engine);

    steps.next();
    const std::vector<int> structureTypes = atomStructureTypes(analysis, static_cast<std::size_t>(frame.natoms));
    json result = perAtomResults(frame, analysis, engine, structureTypes);

    steps.next();
    const Outputs outputs{outputFilename, frame, analysis, engine, structureTypes, result, _inputCrystalStructure, log};
    if(!outputFilename.empty()){
        writeOutputFiles(outputFiles, *this, outputs);
    }
    log.dump(outputs);
    return result;
}

json ElasticStrainService::perAtomResults(
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    const ElasticStrainEngine& engine,
    const std::vector<int>& structureTypes
){
    auto volumetric = engine.volumetricStrains();
    auto strainTensor = engine.strainTensors();
    auto defGrad = engine.deformationGradients();

    const size_t n = static_cast<size_t>(frame.natoms);
    double totalVolumetric = 0.0;
    for(size_t i = 0; i < n; i++){
        if(volumetric){
            totalVolumetric += volumetric->getDouble(i);
        }
    }

    json result;
    result["main_listing"] = {
        { "average_volumetric_strain", n > 0 ? totalVolumetric / n : 0.0 }
    };

    json perAtom = json::array();
    for(size_t i = 0; i < n; i++){
        json a;
        a["id"] = frame.ids[i];

        Cluster* cluster = analysis.atomCluster(static_cast<int>(i));
        a["cluster_id"] = cluster ? cluster->id : 0;
        const int structureType = structureTypes[i];
        a["structure_id"] = structureType;
        a["structure_name"] = StructureIdentificationExport::structureTypeNameForExport(structureType);

        if(volumetric) a["volumetric_strain"] = volumetric->getDouble(i);

        if(strainTensor){
            json tensor = json::array();
            for(int c = 0; c < 6; c++) tensor.push_back(strainTensor->getDoubleComponent(i, c));
            a["strain_tensor"] = std::move(tensor);
        }

        if(defGrad){
            json grad = json::array();
            for(int c = 0; c < 9; c++) grad.push_back(defGrad->getDoubleComponent(i, c));
            a["deformation_gradient"] = std::move(grad);
        }

        perAtom.push_back(std::move(a));
    }
    result["per-atom-properties"] = std::move(perAtom);
    return result;
}

ElasticStrainService::AtomStructures ElasticStrainService::atomStructures(const Outputs& outputs){
    return {
        StructureIdentificationExport::groupAtomsByStructure(outputs.structureTypes, [&](std::size_t atomIndex){
            Cluster* cluster = outputs.analysis.atomCluster(static_cast<int>(atomIndex));
            return cluster && cluster->id != 0;
        }),
        outputs.referenceStructure
    };
}

json ElasticStrainService::atomsTable(const Outputs& outputs, const AtomStructures& structures){
    auto volumetric = outputs.engine.volumetricStrains();
    auto strainTensor = outputs.engine.strainTensors();
    auto defGrad = outputs.engine.deformationGradients();
    json table = StructureIdentificationExport::atomsByStructureTable(outputs.frame, structures.byStructure,
        [&](json& atom, std::size_t atomIndex){
            Cluster* cluster = outputs.analysis.atomCluster(static_cast<int>(atomIndex));
            atom["cluster_id"] = cluster ? cluster->id : 0;
            if(volumetric){
                atom["volumetric_strain"] = volumetric->getDouble(atomIndex);
            }
            if(strainTensor){
                json tensor = json::array();
                for(int c = 0; c < 6; c++)
                    tensor.push_back(strainTensor->getDoubleComponent(atomIndex, c));
                atom["strain_tensor"] = std::move(tensor);
            }
            if(defGrad){
                json grad = json::array();
                for(int c = 0; c < 9; c++)
                    grad.push_back(defGrad->getDoubleComponent(atomIndex, c));
                atom["deformation_gradient"] = std::move(grad);
            }
        });
    table["per-atom-properties"] = outputs.result["per-atom-properties"];
    return table;
}

}
