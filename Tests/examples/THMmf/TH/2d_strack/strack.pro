<?xml version="1.0"?>
<ogs6>
	<coupling>
		<P algorithm="Serial" convergence="FemFunctionConvergenceCheck" max_itr="1" epsilon="1e-4">
			<out>PRESSURE1</out>
			<out>VELOCITY</out>
			<out>TEMPERATURE1</out>
			<problems>
				<M name="LIQUID_FLOW">
					<out>PRESSURE1</out>
				</M>
				<M name="PRESSURE_TO_ELEMENT_VELOCITY">
					<in>PRESSURE1</in>
					<out>VELOCITY</out>
				</M>
				<M name="HEAT_TRANSPORT">
					<in>VELOCITY</in>
					<out>TEMPERATURE1</out>
				</M>
			</problems>
		</P>
	</coupling>
</ogs6>

