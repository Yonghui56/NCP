
OGS_ADD_PROCESS_SYS_SOLVER(LIQUID_FLOW, FunctionLiquidPressure, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
OGS_ADD_PROCESS_SYS_SOLVER(GROUNDWATER_FLOW, FunctionHead, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
OGS_ADD_PROCESS_SYS(HEAD_TO_ELEMENT_VELOCITY, FunctionHeadToElementVelocity, DiscreteLib::DiscreteSystem);
OGS_ADD_PROCESS_SYS(PRESSURE_TO_ELEMENT_VELOCITY, FunctionPressureToElementVelocity, DiscreteLib::DiscreteSystem);
OGS_ADD_PROCESS_SYS(PRESSURE_TO_HEAD, FunctionPressureToHead, DiscreteLib::DiscreteSystem);
OGS_ADD_PROCESS_SYS_SOLVER(MASS_TRANSPORT, FunctionConcentration, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
OGS_ADD_PROCESS_SYS_SOLVER(HEAT_TRANSPORT, FunctionTemperature, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
/*OGS_ADD_PROCESS_SYS_SOLVER(KIN_REACT_GIA, FunctionConcentrations, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
//OGS_ADD_PROCESS_SYS_SOLVER(REACT_GIA, FunctionReductConc, DiscreteLib::DiscreteSystem, MathLib::EigenDenseLinearEquation);
OGS_ADD_PROCESS_SYS_SOLVER(REACT_TRANS_OPS, FunctionOPSConc, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
//OGS_ADD_PROCESS_SYS_SOLVER(REACT_GIA, FunctionReductConc, DiscreteLib::DiscreteSystem, MathLib::DenseLinearEquation);   // RZ: using direct solver for global problem.
OGS_ADD_PROCESS_SYS_SOLVER(REACT_GIA, FunctionReductConc, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);*/
//OGS_ADD_PROCESS_SYS_SOLVER(DEFORMATION, FunctionDisplacement, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
//OGS_ADD_PROCESS_SYS(ELEMENT_STRESS_STRAIN, FunctionElementStressStrain, DiscreteLib::DiscreteSystem);
//OGS_ADD_PROCESS_SYS(NODAL_STRESS_STRAIN, FunctionNodalStressStrain, DiscreteLib::DiscreteSystem);
//OGS_ADD_PROCESS_SYS_SOLVER(DEFORMATION_FLOW, FunctionDisplacementPressure, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
/*
OGS_ADD_PROCESS(XFEM_EXAMPLE_CRACK1, xfem::FunctionXFEM_EXAMPLE_CRACK1);
OGS_ADD_PROCESS_SYS_SOLVER(LIQUID_FLOW, THMmf::Pmf, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
OGS_ADD_PROCESS_SYS_SOLVER(HEAT_TRANSPORT, THMmf::Tmf, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
OGS_ADD_PROCESS_SYS_SOLVER(DEFORMATION, THMmf::Umf, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
OGS_ADD_PROCESS_SYS_SOLVER(INCREMENTAL_DEFORMATION, FunctionIncrementalDisplacement, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
*/
OGS_ADD_PROCESS_SYS_SOLVER(RICHARDS_FLOW, FunctionRichards, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
/* this following process is component based multphase flow */
OGS_ADD_PROCESS_SYS_SOLVER(CMP_2P2C, FunctionCMP_2P2C, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);//
OGS_ADD_PROCESS_SYS_SOLVER(CMP_CapPresFormular, FunctionCMP_CapPresFormular, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);//
OGS_ADD_PROCESS_SYS_SOLVER(CMP_PressureForm, FunctionCMP_PressureForm, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);//MathLib::EigenDenseLinearEquation
OGS_ADD_PROCESS_SYS_SOLVER(CMP_TotalDensityForm, FunctionCMP_TotalDensityForm, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);//MathLib::EigenDenseLinearEquation
OGS_ADD_PROCESS_SYS_SOLVER(CMP_GlobalComplementaryForm, FunctionCMP_GlobalComplementaryForm, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);//MathLib::EigenDenseLinearEquation