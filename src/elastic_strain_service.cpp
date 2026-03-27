#include <volt/elastic_strain_service.h>
#include <volt/elastic_strain_engine.h>
#include <volt/core/reconstructed_structure.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/analysis/cluster_graph_export.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <utility>

namespace Volt{

using namespace Volt::Particles;

namespace {

std::string structureTypeNameForExport(int structureType){
    switch(static_cast<StructureType>(structureType)){
        case StructureType::SC:
            return "SC";
        case StructureType::FCC:
            return "FCC";
        case StructureType::HCP:
            return "HCP";
        case StructureType::BCC:
            return "BCC";
        case StructureType::CUBIC_DIAMOND:
            return "CUBIC_DIAMOND";
        case StructureType::HEX_DIAMOND:
            return "HEX_DIAMOND";
        case StructureType::ICO:
            return "ICO";
        case StructureType::GRAPHENE:
            return "GRAPHENE";
        case StructureType::CUBIC_DIAMOND_FIRST_NEIGH:
            return "CUBIC_DIAMOND_FIRST_NEIGH";
        case StructureType::CUBIC_DIAMOND_SECOND_NEIGH:
            return "CUBIC_DIAMOND_SECOND_NEIGH";
        case StructureType::HEX_DIAMOND_FIRST_NEIGH:
            return "HEX_DIAMOND_FIRST_NEIGH";
        case StructureType::HEX_DIAMOND_SECOND_NEIGH:
            return "HEX_DIAMOND_SECOND_NEIGH";
        case StructureType::OTHER:
        case StructureType::NUM_STRUCTURE_TYPES:
        default:
            return "OTHER";
    }
}

std::vector<int> buildAtomStructureTypes(const StructureAnalysis& analysis, std::size_t atomCount){
    std::vector<int> structureTypes(atomCount, static_cast<int>(StructureType::OTHER));
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        Cluster* cluster = analysis.atomCluster(static_cast<int>(atomIndex));
        if(cluster){
            structureTypes[atomIndex] = cluster->structure;
        }
    }
    return structureTypes;
}

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

void ElasticStrainService::setParameters(
    double latticeConstant,
    double caRatio,
    bool pushForward,
    bool calculateDeformationGradient,
    bool calculateStrainTensors
){
    _latticeConstant = latticeConstant;
    _caRatio = caRatio;
    _pushForward = pushForward;
    _calculateDeformationGradient = calculateDeformationGradient;
    _calculateStrainTensors = calculateStrainTensors;
}

json ElasticStrainService::compute(const LammpsParser::Frame &frame, const std::string &outputFilename){
    auto startTime = std::chrono::high_resolution_clock::now();

    if(_clustersTablePath.empty() || _clusterTransitionsPath.empty()){
        return AnalysisResult::failure(
            "ElasticStrain requires --clusters-table and --clusters-transitions"
        );
    }

    FrameAdapter::PreparedAnalysisInput prepared;
    std::string frameError;
    if(!FrameAdapter::prepareAnalysisInput(frame, prepared, &frameError))
        return AnalysisResult::failure(frameError);

    auto positions = std::move(prepared.positions);

    ReconstructedStructureContext context(positions.get(), frame.simulationCell);
    context.inputCrystalType = _inputCrystalStructure;

    StructureAnalysis analysis(context);
    std::string reconstructionError;
    if(!ReconstructedStructureLoader::load(
        frame,
        {_clustersTablePath, _clusterTransitionsPath},
        analysis,
        context,
        &reconstructionError
    )){
        return AnalysisResult::failure(reconstructionError);
    }
    rebuildImportedClusterParentHierarchy(analysis);

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

    auto volumetric = engine.volumetricStrains();
    auto strainTensor = engine.strainTensors();
    auto defGrad = engine.deformationGradients();

    const size_t n = static_cast<size_t>(frame.natoms);
    const std::vector<int> atomStructureTypes = buildAtomStructureTypes(analysis, n);
    double totalVolumetric = 0.0;
    for(size_t i = 0; i < n; i++){
        if(volumetric) totalVolumetric += volumetric->getDouble(i);
    }

    json result;
    result["main_listing"] = {
        { "average_volumetric_strain", n > 0 ? totalVolumetric / n : 0.0 }
    };

    json perAtom = json::array();
    for(size_t i = 0; i < n; i++){
        json a;
        a["id"] = frame.ids[i];

        if(volumetric) a["volumetric_strain"] = volumetric->getDouble(i);

        if(strainTensor){
            json tensor = json::array();
            for(int c = 0; c < 6; c++) tensor.push_back(strainTensor->getDoubleComponent(i, c));
            a["strain_tensor"] = tensor;
        }

        if(defGrad){
            json grad = json::array();
            for(int c = 0; c < 9; c++) grad.push_back(defGrad->getDoubleComponent(i, c));
            a["deformation_gradient"] = grad;
        }

        perAtom.push_back(a);
    }
    result["per-atom-properties"] = perAtom;

    if(!outputFilename.empty()){
        const std::string outputPath = outputFilename + "_elastic_strain.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(result, outputPath, false)){
            spdlog::info("Elastic strain msgpack written to {}", outputPath);
        }else{
            spdlog::warn("Could not write elastic strain msgpack: {}", outputPath);
        }

        // --- atoms.msgpack export (Structure Identification exposure) ---
        {
            constexpr int K = static_cast<int>(StructureType::NUM_STRUCTURE_TYPES);
            std::vector<std::string> names(K);
            for(int st = 0; st < K; st++)
                names[st] = structureTypeNameForExport(st);

            std::vector<std::vector<size_t>> structureAtomIndices(K);
            for(size_t i = 0; i < n; ++i){
                const int raw = atomStructureTypes[i];
                const int st = (0 <= raw && raw < K) ? raw : static_cast<int>(StructureType::OTHER);
                structureAtomIndices[static_cast<size_t>(st)].push_back(i);
            }

            std::vector<int> structureOrder;
            structureOrder.reserve(K);
            for(int st = 0; st < K; st++){
                if(!structureAtomIndices[static_cast<size_t>(st)].empty())
                    structureOrder.push_back(st);
            }
            std::sort(structureOrder.begin(), structureOrder.end(),
                [&](int a, int b){ return names[a] < names[b]; });

            json atomsByStructure;
            for(int st : structureOrder){
                json atomsArray = json::array();
                for(size_t atomIndex : structureAtomIndices[static_cast<size_t>(st)]){
                    const Point3& pos = frame.positions[atomIndex];
                    atomsArray.push_back({
                        {"id", frame.ids[atomIndex]},
                        {"pos", {pos.x(), pos.y(), pos.z()}}
                    });
                }
                atomsByStructure[names[st]] = atomsArray;
            }

            json exportWrapper;
            exportWrapper["export"] = json::object();
            exportWrapper["export"]["AtomisticExporter"] = atomsByStructure;
            const std::string atomsPath = outputFilename + "_atoms.msgpack";
            if(JsonUtils::writeJsonMsgpackToFile(exportWrapper, atomsPath, false)){
                spdlog::info("Exported atoms data to: {}", atomsPath);
            }else{
                spdlog::warn("Could not write atoms msgpack: {}", atomsPath);
            }
        }
    }

    spdlog::info("Elastic strain analysis completed.");
    return result;
}
}
