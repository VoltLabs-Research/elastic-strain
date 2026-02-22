#include <volt/elastic_strain_service.h>
#include <volt/elastic_strain_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/utilities/concurrence/parallel_system.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>

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

    if(frame.natoms <= 0){
        return AnalysisResult::failure("Invalid number of atoms");
    }

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions){
        return AnalysisResult::failure("Failed to create position property");
    }

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

    json result = AnalysisResult::success();
    AnalysisResult::addTiming(result, startTime);

    if(!outputFilename.empty()){
        const std::string outputPath = outputFilename + "_elastic_strain.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(result, outputPath, false)){
            spdlog::info("Elastic strain msgpack written to {}", outputPath);
        }else{
            spdlog::warn("Could not write elastic strain msgpack: {}", outputPath);
        }
    }

    return result;
}
}
