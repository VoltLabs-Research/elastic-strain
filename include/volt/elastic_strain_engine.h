#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/core/particle_property.h>
#include <volt/core/simulation_cell.h>

#include <memory>

namespace Volt{

class ElasticStrainEngine{
public:
    ElasticStrainEngine(
        StructureAnalysis& structureAnalysis,
        StructureContext& context,
        LatticeStructureType inputCrystalStructure,
        bool calculateDeformationGradients,
        bool calculateStrainTensors,
        double latticeConstant,
        double caRatio,
        bool pushStrainTensorsForward
    );

    void perform();

    ParticleProperty* atomClusters() const{
        return _context.atomClusters ? _context.atomClusters.get() : nullptr;
    }

    const StructureAnalysis& structureAnalysis() const{
        return _structureAnalysis;
    }

    ParticleProperty* volumetricStrains() const{
        return _volumetricStrains.get();
    }

    ParticleProperty* strainTensors() const{
        return _strainTensors.get();
    }

    ParticleProperty* deformationGradients() const{
        return _deformationGradients.get();
    }

private:
    double _latticeConstant;
    double _axialScaling;
    LatticeStructureType _inputCrystalStructure;
    bool _pushStrainTensorsForward;

    StructureContext& _context;
    StructureAnalysis& _structureAnalysis;

    std::unique_ptr<ParticleProperty> _volumetricStrains;
    std::unique_ptr<ParticleProperty> _strainTensors;
    std::unique_ptr<ParticleProperty> _deformationGradients;
};

}
