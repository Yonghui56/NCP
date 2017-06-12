<?xml version="1.0"?>
<ogs6>
<coupling>
    <P algorithm="Serial" convergence="FemFunctionConvergenceCheck" max_itr="1" epsilon="1e-4">
        <out>MEAN_PRESSURE</out>
        <out>MOLAR_FRACTION</out>
        <problems>
            <M name="CMP_2P2C" type="CMP_2P2C">
                <out>MEAN_PRESSURE</out>
				<out>MOLAR_FRACTION</out>
            </M>
        </problems>
    </P>
</coupling>
</ogs6> 

