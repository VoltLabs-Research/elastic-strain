#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <nlohmann/json.hpp>
#include <volt/core/particle_property.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/plugin/output_files.h>
#include <volt/structures/crystal_structure_types.h>
#include <cstddef>
#include <string>
#include <vector>

namespace Volt{
using json = nlohmann::json;

class ElasticStrainEngine;
class ElasticStrainLog;
class StructureAnalysis;

class ElasticStrainService{
public:
    // What the files are written from: where they go, the strain the engine computed on the upstream structures, the
    // structure type of each atom, the per-atom results, and the log that names each file once it is written.
    struct Outputs{
        const std::string& base;
        const LammpsParser::Frame& frame;
        const StructureAnalysis& analysis;
        const ElasticStrainEngine& engine;
        const std::vector<int>& structureTypes;
        const json& result;
        LatticeStructureType referenceStructure;
        ElasticStrainLog& log;
    };

    // The atoms by structure, and the structure the strain is measured against.
    struct AtomStructures{
        StructureIdentificationExport::AtomsByStructure byStructure;
        LatticeStructureType referenceStructure = LATTICE_OTHER;
    };

    ElasticStrainService();

    void setInputCrystalStructure(LatticeStructureType structure);
    void setClustersTablePath(std::string path);
    void setClusterTransitionsPath(std::string path);
    void setNeighborLatticePath(std::string path);
    void setLatticeConstant(double latticeConstant);
    void setCaRatio(double caRatio);
    void setPushForward(bool pushForward);
    void setCalculateDeformationGradient(bool calculateDeformationGradient);
    void setCalculateStrainTensors(bool calculateStrainTensors);

    json compute(
        const LammpsParser::Frame& frame,
        const std::string& outputFilename = ""
    );

private:
    // The files in the order they are written. None is optional, so no flag of the service asks for them.
    static const OutputFile<ElasticStrainService, Outputs> outputFiles[];

    // The structure type of each atom's cluster; OTHER for an atom of no cluster or of a type with no export name.
    static std::vector<int> atomStructureTypes(const StructureAnalysis& analysis, std::size_t atomCount);
    static json perAtomResults(
        const LammpsParser::Frame& frame,
        const StructureAnalysis& analysis,
        const ElasticStrainEngine& engine,
        const std::vector<int>& structureTypes
    );
    static AtomStructures atomStructures(const Outputs& outputs);
    static json atomsTable(const Outputs& outputs, const AtomStructures& structures);

    LatticeStructureType _inputCrystalStructure;
    std::string _clustersTablePath;
    std::string _clusterTransitionsPath;
    std::string _neighborLatticePath;

    double _latticeConstant;
    double _caRatio;
    bool _pushForward;
    bool _calculateDeformationGradient;
    bool _calculateStrainTensors;
};

}
