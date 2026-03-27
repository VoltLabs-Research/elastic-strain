#include <volt/elastic_strain_engine.h>

#include <cmath>
#include <cassert>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace Volt{

ElasticStrainEngine::ElasticStrainEngine(
    StructureAnalysis& structureAnalysis,
    StructureContext& context,
    LatticeStructureType inputCrystalStructure,
    bool calculateDeformationGradients,
    bool calculateStrainTensors,
    double latticeConstant,
    double caRatio,
    bool pushStrainTensorsForward
)
    : _latticeConstant(latticeConstant)
    , _axialScaling(1.0)
    , _inputCrystalStructure(inputCrystalStructure)
    , _pushStrainTensorsForward(pushStrainTensorsForward)
    , _context(context)
    , _structureAnalysis(structureAnalysis)
    , _volumetricStrains(std::make_unique<ParticleProperty>(
          context.atomCount(), DataType::Double, 1, 0, false))
    , _strainTensors(calculateStrainTensors
          ? std::make_unique<ParticleProperty>(
                context.atomCount(), DataType::Double, 6, 0, false)
          : nullptr)
    , _deformationGradients(calculateDeformationGradients
          ? std::make_unique<ParticleProperty>(
                context.atomCount(), DataType::Double, 9, 0, false)
          : nullptr)
{
    if(inputCrystalStructure == LatticeStructureType::LATTICE_FCC ||
       inputCrystalStructure == LatticeStructureType::LATTICE_BCC ||
       inputCrystalStructure == LatticeStructureType::LATTICE_CUBIC_DIAMOND){
        _axialScaling = 1.0;
    }else{
        _latticeConstant *= std::sqrt(2.0);
        _axialScaling = caRatio / std::sqrt(8.0 / 3.0);
    }
}

void ElasticStrainEngine::perform(){
    const std::size_t N = _context.atomCount();

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, N),
        [this](const tbb::blocked_range<std::size_t>& r){
        for(std::size_t particleIndex = r.begin(); particleIndex < r.end(); ++particleIndex){

        Cluster* localCluster = _structureAnalysis.atomCluster(static_cast<int>(particleIndex));
        if(!localCluster || localCluster->id == 0){
            _volumetricStrains->setDouble(particleIndex, 0.0);
            if(_strainTensors){
                for(size_t c = 0; c < 6; ++c){
                    _strainTensors->setDoubleComponent(particleIndex, c, 0.0);
                }
            }
            if(_deformationGradients){
                for(size_t c = 0; c < 9; ++c){
                    _deformationGradients->setDoubleComponent(particleIndex, c, 0.0);
                }
            }
            continue;
        }

        Matrix3 idealUnitCellTM(
            _latticeConstant, 0.0,            0.0,
            0.0,              _latticeConstant, 0.0,
            0.0,              0.0,            _latticeConstant * _axialScaling
        );

        Cluster* parentCluster = localCluster;
        ClusterTransition* parentTransition = localCluster->parentTransition;
        while(parentTransition != nullptr){
            idealUnitCellTM = idealUnitCellTM * parentTransition->tm;
            parentCluster = parentTransition->cluster2;
            parentTransition = parentCluster ? parentCluster->parentTransition : nullptr;
        }

        if(parentCluster == localCluster && localCluster->structure != _inputCrystalStructure){
            parentCluster = nullptr;
        }

        if(!parentCluster){
            _volumetricStrains->setDouble(particleIndex, 0.0);
            if(_strainTensors){
                for(size_t c = 0; c < 6; ++c){
                    _strainTensors->setDoubleComponent(particleIndex, c, 0.0);
                }
            }
            if(_deformationGradients){
                for(size_t c = 0; c < 9; ++c){
                    _deformationGradients->setDoubleComponent(particleIndex, c, 0.0);
                }
            }
            continue;
        }

        if(parentCluster->structure != _inputCrystalStructure){
            _volumetricStrains->setDouble(particleIndex, 0.0);
            if(_strainTensors){
                for(size_t c = 0; c < 6; ++c){
                    _strainTensors->setDoubleComponent(particleIndex, c, 0.0);
                }
            }
            if(_deformationGradients){
                for(size_t c = 0; c < 9; ++c){
                    _deformationGradients->setDoubleComponent(particleIndex, c, 0.0);
                }
            }
            continue;
        }

        Matrix_3<double> orientationV = Matrix_3<double>::Zero();
        Matrix_3<double> orientationW = Matrix_3<double>::Zero();

        int numneighs = _structureAnalysis.numberOfNeighbors(static_cast<int>(particleIndex));
        for(int n = 0; n < numneighs; ++n){
            int neighborAtomIndex = _structureAnalysis.getNeighbor(static_cast<int>(particleIndex), n);

            Vector3 latticeVector =
                idealUnitCellTM * _structureAnalysis.neighborLatticeVector(static_cast<int>(particleIndex), n);

            const Vector3& spatialVector =
                _context.simCell.wrapVector(
                    _context.positions->getPoint3(neighborAtomIndex) -
                    _context.positions->getPoint3(particleIndex)
                );

            for(size_t r = 0; r < 3; ++r){
                for(size_t c = 0; c < 3; ++c){
                    orientationV(r,c) += static_cast<double>(latticeVector[c] * latticeVector[r]);
                    orientationW(r,c) += static_cast<double>(latticeVector[c] * spatialVector[r]);
                }
            }
        }

        Matrix_3<double> elasticF = orientationW * orientationV.inverse();

        if(_deformationGradients){
            for(size_t col = 0; col < 3; ++col){
                for(size_t row = 0; row < 3; ++row){
                    _deformationGradients->setDoubleComponent(
                        particleIndex, col*3 + row, static_cast<double>(elasticF(row,col)));
                }
            }
        }

        SymmetricTensor2T<double> elasticStrain;
        if(!_pushStrainTensorsForward){
            elasticStrain = (Product_AtA(elasticF) - SymmetricTensor2T<double>::Identity()) * 0.5;
        }else{
            Matrix_3<double> inverseF;
            if(!elasticF.inverse(inverseF)){
                _volumetricStrains->setDouble(particleIndex, 0.0);
                if(_strainTensors){
                    for(size_t c = 0; c < 6; ++c){
                        _strainTensors->setDoubleComponent(particleIndex, c, 0.0);
                    }
                }
                continue;
            }
            elasticStrain = (SymmetricTensor2T<double>::Identity() - Product_AtA(inverseF)) * 0.5;
        }

        if(_strainTensors){
            _strainTensors->setSymmetricTensor2(particleIndex, (SymmetricTensor2)elasticStrain);
        }

        double volumetricStrain =
            (elasticStrain(0,0) + elasticStrain(1,1) + elasticStrain(2,2)) / 3.0;
        assert(std::isfinite(volumetricStrain));
        _volumetricStrains->setDouble(particleIndex, volumetricStrain);
    }
    });
}

}
