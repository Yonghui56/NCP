<?xml version="1.0"?>
<ogs6>
	<coupling>
		<P algorithm="Serial" convergence="FemFunctionConvergenceCheck" max_itr="1" epsilon="1e-4">
			<out>DISPLACEMENT</out>
			<out>gp_strain</out>
			<out>gp_stress</out>
			<out>STRAIN</out>
			<out>STRESS</out>
			<problems>
				<M name="DEFORMATION">
					<out>DISPLACEMENT</out>
				</M>
				<M name="ELEMENT_STRESS_STRAIN">
					<in>DISPLACEMENT</in>
					<out>gp_strain</out>
					<out>gp_stress</out>
				</M>
				<M name="NODAL_STRESS_STRAIN">
					<in>gp_strain</in>
					<in>gp_stress</in>
					<out>STRAIN</out>
					<out>STRESS</out>
				</M>
			</problems>
		</P>
	</coupling>
</ogs6>

