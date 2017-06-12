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
				<M type="INCREMENTAL_DEFORMATION" name="DEFORMATION">
					<out>DISPLACEMENT</out>
					<out>gp_strain</out>
					<out>gp_stress</out>
				</M>
				<M type="NODAL_STRESS_STRAIN">
					<in>gp_strain</in>
					<in>gp_stress</in>
					<out>STRAIN</out>
					<out>STRESS</out>
				</M>
			</problems>
		</P>
	</coupling>
</ogs6>

