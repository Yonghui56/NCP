GeoSys-ST: Source Terms ------------------------------------------------
#STOP

 $PCS_TYPE
  CMP_PressureForm
 $PRIMARY_VARIABLE
  LIQUID_PRESSURE
 $GEO_TYPE
  POLYLINE PLY_LEFT
 $DIS_TYPE   
  CONSTANT_NEUMANN 0.00005 
 $DIS_TYPE  
  CONSTANT_NEUMANN 0.0002

 $DIS_TYPE  
  RIVER 2
  4  3 1.00e-6 10.0 1.3 1.0
  5  3 1.00e-6 10.0 1.3 1.0