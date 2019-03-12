#include "GaudiKernel/DeclareFactoryEntries.h"
#include "SYJ_Lambda_c/Sigma0PionPi0.h"
#include "SYJ_Lambda_c/Sigma0PionEta.h"
#include "SYJ_Lambda_c/Sigma0PipPipPim.h"
#include "SYJ_Lambda_c/Sigma0Pip.h"

DECLARE_ALGORITHM_FACTORY( Sigma0PionPi0 )
DECLARE_ALGORITHM_FACTORY( Sigma0PionEta )
DECLARE_ALGORITHM_FACTORY( Sigma0PipPipPim )
DECLARE_ALGORITHM_FACTORY( Sigma0Pip )

DECLARE_FACTORY_ENTRIES( SYJ_Lambda_c ) {
DECLARE_ALGORITHM( Sigma0PionPi0 );
DECLARE_ALGORITHM( Sigma0PionEta );
DECLARE_ALGORITHM( Sigma0PipPipPim );
DECLARE_ALGORITHM( Sigma0Pip);
}

