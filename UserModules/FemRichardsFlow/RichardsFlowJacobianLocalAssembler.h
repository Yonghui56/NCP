#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Compound.h"

#include "Ogs6FemData.h"

class RichardsFlowJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
    RichardsFlowJacobianLocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects)
    : _feObjects(*feObjects)
    {};

    virtual ~RichardsFlowJacobianLocalAssembler() {};

    void assembly(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &/*e*/, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/, MathLib::LocalMatrix &/*localJ*/)
    {
    }

private:
    FemLib::LagrangeFeObjectContainer _feObjects;
};
