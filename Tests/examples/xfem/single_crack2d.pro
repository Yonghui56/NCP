<?xml version="1.0"?>
<ogs6>
<!--
	<meshList>
		<mesh id="0" file="single_crack2d.msh">
	</meshList>
	<timeStepList>
		<timeStep id="0" uni="SECOND" start="0" end="1" control="STEPS">
			<step dt="1" repeat="1"/>
		</timeStep>
	</timeStepList>
-->
	<processList>
		<process type="XFEM_EXAMPLE_CRACK1" name="XFEM_EXAMPLE_CRACK1" MeshID="0" TimeGroupID="0" />
	</processList>
	<coupling>
		<M name="XFEM_EXAMPLE_CRACK1" />
	</coupling>
	<outputList>
		<output dataType="PVD" meshID="0" geoType="DOMAIN" geoName="" timeType="STEPS" timeSteps="1">
			<nodeValue name="DISPLACEMENT" />
			<nodeValue name="EXACT_U" />
		</output>
	</outputList>
</ogs6>

