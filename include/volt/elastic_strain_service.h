#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <nlohmann/json.hpp>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_structure_types.h>
#include <string>

namespace Volt{
using json = nlohmann::json;
    
class ElasticStrainService{
public:
    ElasticStrainService();

    void setInputCrystalStructure(LatticeStructureType structure);
    void setClustersTablePath(std::string path);
    void setClusterTransitionsPath(std::string path);

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
    std::string _clustersTablePath;
    std::string _clusterTransitionsPath;

    double _latticeConstant;
    double _caRatio;
    bool _pushForward;
    bool _calculateDeformationGradient;
    bool _calculateStrainTensors;
};

}
