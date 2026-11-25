#pragma once

#include <volt/core/run_log.h>
#include <volt/elastic_strain_service.h>

namespace Volt{

class ElasticStrainEngine;
class StructureAnalysis;

// The execution log of one elastic strain analysis. Each phase hands it what it built, and the log writes what that
// holds and how long the phase took. Scripts parse these lines, so their wording stays as it is.
class ElasticStrainLog : public RunLog{
public:
    using RunLog::dump;

    // The structure identification loaded from the upstream contract.
    void dump(const StructureAnalysis& loaded);
    // The strain, with the range of the volumetric strain.
    void dump(const ElasticStrainEngine& engine);
    void dump(const ElasticStrainService::AtomStructures& structures);
    // The outputs once written; the line times the per-atom results with their files.
    void dump(const ElasticStrainService::Outputs& outputs);

private:
    // Since the strain was computed, when the per-atom results began.
    Stopwatch _results;
};

}
