#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/core/particle_property.h>
#include <volt/core/simulation_cell.h>
#include <volt/analysis/analysis_context.h>
#include <volt/structures/cluster_graph.h>

#include <memory>
#include <vector>

namespace Volt{

class ElasticStrainEngine{
public:
    ElasticStrainEngine(
        ParticleProperty* positions,
        ParticleProperty* structures,
        const SimulationCell& simcell,
        LatticeStructureType inputCrystalStructure,
        std::vector<Matrix3>&& preferredCrystalOrientations,
        bool calculateDeformationGradients,
        bool calculateStrainTensors,
        double latticeConstant,
        double caRatio,
        bool pushStrainTensorsForward,
        StructureAnalysis::Mode identificationMode = StructureAnalysis::Mode::PTM,
        double rmsd = 0.12
    );

    void perform();

    // Returns the array of atom cluster IDs
    ParticleProperty* atomClusters() const{
        return _context.atomClusters ? _context.atomClusters.get() : nullptr;
    }

    // Returns the created cluster graph
    ClusterGraph* clusterGraph(){
        return &_structureAnalysis.clusterGraph();
    }

    const StructureAnalysis& structureAnalysis() const{
        return _structureAnalysis;
    }

    // Returns the property storage that contains the computed per-particle volumetric strain values
    ParticleProperty* volumetricStrains() const{
        return _volumetricStrains.get();
    }

    // Returns the property storage that contains the computed per-particle strain tensors
    ParticleProperty* strainTensors() const{
        return _strainTensors.get();
    }

    // Returns the property storage that contains the computed per-particle deformation gradient tensors
    ParticleProperty* deformationGradients() const{
        return _deformationGradients.get();
    }

private:
    double _latticeConstant;
    double _axialScaling;
    LatticeStructureType _inputCrystalStructure;
    bool _pushStrainTensorsForward;

    AnalysisContext   _context;  
    StructureAnalysis _structureAnalysis;

    std::unique_ptr<ParticleProperty> _volumetricStrains;
    std::unique_ptr<ParticleProperty> _strainTensors;
    std::unique_ptr<ParticleProperty> _deformationGradients;
};

}
