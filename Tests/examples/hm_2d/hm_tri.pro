<?xml version="1.0"?>
<ogs6>
<coupling>
	<P algorithm="Serial" convergence="FemFunctionConvergenceCheck" max_itr="1" epsilon="1e-4">
		<out>DISPLACEMENT</out>
		<out>PRESSURE1</out>
		<out>VELOCITY</out>
		<out>gp_strain</out>
		<out>gp_stress</out>
		<out>STRAIN</out>
		<out>STRESS</out>
		<problems>
			<M name="DEFORMATION_FLOW">
				<out>DISPLACEMENT</out>
				<out>PRESSURE1</out>
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
			<M name="PRESSURE_TO_ELEMENT_VELOCITY">
				<in>PRESSURE1</in>
				<out>VELOCITY</out>
			</M>
			<!--
			<M name="PRESSURE_TO_HEAD">
				<in>PRESSURE1</in>
				<out>HEAD</out>
			</M>
			<M name="HEAD_TO_ELEMENT_VELOCITY">
				<in>HEAD</in>
				<out>VELOCITY</out>
			</M>
			-->
		</problems>
	</P>
</coupling>
</ogs6>

