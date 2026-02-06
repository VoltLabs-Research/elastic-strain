#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <nlohmann/json.hpp>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/analysis/structure_analysis.h>
#include <string>

namespace Volt{
using json = nlohmann::json;
    
class ElasticStrainService{
public:
    ElasticStrainService();

    void setInputCrystalStructure(LatticeStructureType structure);
    void setIdentificationMode(StructureAnalysis::Mode mode);
    void setRMSD(float rmsd);

    void setParameters(
        double latticeConstant,
        double caRatio,
        bool pushForward,
        bool calculateDeformationGradient,
        bool calculateStrainTensors
    );

    json compute(
        const LammpsParser::Frame& frame,
        const std::string& outputFilename = ""
    );

private:
    LatticeStructureType _inputCrystalStructure;
    StructureAnalysis::Mode _identificationMode;
    float _rmsd;

    double _latticeConstant;
    double _caRatio;
    bool _pushForward;
    bool _calculateDeformationGradient;
    bool _calculateStrainTensors;
};

}
