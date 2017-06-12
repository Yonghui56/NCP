<?xml version="1.0"?>
<ogs6>
<coupling>
    <P algorithm="Serial" convergence="FemFunctionConvergenceCheck" max_itr="1" epsilon="1e-4">
        <out>NONWET_PRESSURE</out>
        <out>CAPILLARY_PRESSURE</out>
        <problems>
            <M name="CMP_CapPresFormular" type="CMP_CapPresFormular">
                <out>NONWET_PRESSURE</out>
				<out>CAPILLARY_PRESSURE</out>
            </M>
        </problems>
    </P>
</coupling>
</ogs6> 

