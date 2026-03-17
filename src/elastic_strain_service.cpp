#include <volt/elastic_strain_service.h>
#include <volt/elastic_strain_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/utilities/concurrence/parallel_system.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>
#include <algorithm>

namespace Volt{

using namespace Volt::Particles;

ElasticStrainService::ElasticStrainService()
    : _inputCrystalStructure(LATTICE_BCC),
      _identificationMode(StructureAnalysis::Mode::PTM),
      _rmsd(0.10f),
      _latticeConstant(1.63f),
      _caRatio(1.0),
      _pushForward(false),
      _calculateDeformationGradient(true),
      _calculateStrainTensors(true){}

void ElasticStrainService::setInputCrystalStructure(LatticeStructureType structure){
    _inputCrystalStructure = structure;
}

void ElasticStrainService::setIdentificationMode(StructureAnalysis::Mode mode){
    _identificationMode = mode;
}

void ElasticStrainService::setRMSD(float rmsd){
    _rmsd = rmsd;
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

    if(frame.natoms <= 0)
        return AnalysisResult::failure("Invalid number of atoms");

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions)
        return AnalysisResult::failure("Failed to create position property");

    std::vector<Matrix3> preferredOrientations;
    preferredOrientations.push_back(Matrix3::Identity());

    auto structureTypes = std::make_shared<ParticleProperty>(frame.natoms, DataType::Int, 1, 0, true);
    ElasticStrainEngine engine(
        positions.get(),
        structureTypes.get(),
        frame.simulationCell,
        static_cast<LatticeStructureType>(_inputCrystalStructure),
        std::move(preferredOrientations),
        _calculateDeformationGradient,
        _calculateStrainTensors,
        _latticeConstant,
        _caRatio,
        _pushForward,
        _identificationMode,
        _rmsd
    );

    engine.perform();

    auto volumetric = engine.volumetricStrains();
    auto strainTensor = engine.strainTensors();
    auto defGrad = engine.deformationGradients();

    const size_t n = static_cast<size_t>(frame.natoms);
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
            const StructureAnalysis& sa = engine.structureAnalysis();
            constexpr int K = static_cast<int>(StructureType::NUM_STRUCTURE_TYPES);
            std::vector<std::string> names(K);
            for(int st = 0; st < K; st++)
                names[st] = sa.getStructureTypeName(st);

            std::vector<std::vector<size_t>> structureAtomIndices(K);
            for(size_t i = 0; i < n; ++i){
                const int raw = structureTypes->getInt(i);
                const int st = (0 <= raw && raw < K) ? raw : 0;
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

