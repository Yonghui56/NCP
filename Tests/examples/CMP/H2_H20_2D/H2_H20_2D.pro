<?xml version="1.0"?>
<ogs6>
<coupling>
    <P algorithm="Serial" convergence="FemFunctionConvergenceCheck" max_itr="1" epsilon="1e-4">
        <out>LIQUID_PRESSURE</out>
        <out>MOLAR_FRACTION</out>
		<out>SATURATION</out>
        <problems>
            <M name="CMP_GlobalComplementaryForm" type="CMP_GlobalComplementaryForm">
                <out>LIQUID_PRESSURE</out>
				<out>MOLAR_FRACTION</out>
				<out>SATURATION</out>
            </M>
        </problems>
    </P>
</coupling>
</ogs6> 

