/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file THMCSimulator.cpp
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#include "THMCSimulator.h"

#include <iostream>

// external library
#include "logog.hpp"
#include "tclap/CmdLine.h"
#ifdef USE_LIS
#include "lis.h"
#endif
#ifdef USE_PETSC
#include "petscksp.h"
#endif

// internal library
#include "BaseLib/CodingTools.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Options.h"
#include "BaseLib/OptionsXMLReader.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#include "NumLib/TransientCoupling/TransientCouplingStructureBuilder.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "FemIO/ogs5/Ogs5FemIO.h"

// this module
#include "SimulationInfo.h"
#include "FormatterCustom.h"
#include "Ogs6FemData.h"
#include "Ogs5ToOgs6.h"
#include "TimeSteppingControllerWithOutput.h"

namespace ogs6
{

////////////////////////////////////////////////////////////////////////////////
// Variables
////////////////////////////////////////////////////////////////////////////////
static logog::Cout *logogCout;
static logog::LogFile *logog_file;
static FormatterCustom *custom_format;
static bool isOgsInitCalled = false;
static bool isOgsExitCalled = false;
////////////////////////////////////////////////////////////////////////////////

void ogsInit(int argc, char* argv[])
{
    if (isOgsInitCalled) return;
    isOgsInitCalled = true;

    LOGOG_INITIALIZE();
    custom_format = new FormatterCustom();
    logogCout = new logog::Cout();
    logogCout->SetFormatter(*custom_format);
    logog_file = NULL;

#ifdef USE_LIS
    lis_initialize((LIS_INT*)&argc, &argv);
#endif
#ifdef USE_PETSC
    char petsc_help[] = "Using PETSc package\n";
    PetscInitialize(&argc, &argv,(char *)0,petsc_help);
#endif
}

void ogsExit()
{
    if (isOgsExitCalled) return;
    isOgsExitCalled = true;

#ifdef USE_LIS
    lis_finalize();
#endif
#ifdef USE_PETSC
    PetscFinalize();
#endif

    INFO("exit ogs6.");
    BaseLib::releaseObject(custom_format, logogCout, logog_file);
    LOGOG_SHUTDOWN();
}

THMCSimulator::THMCSimulator(int argc, char* argv[])
: _sim_info(NULL), _cpl_system(NULL)
{
    try {
        // Command line parser
        TCLAP::CmdLine cmd("ogs6", ' ', "0.1");
        TCLAP::ValueArg<std::string> input_arg("i", "input", "input file", false, "", "string");
        cmd.add( input_arg );
        TCLAP::ValueArg<std::string> output_dir_arg("o", "output", "output directory", false, "", "string");
        cmd.add( output_dir_arg );
        TCLAP::ValueArg<std::string> logfile_arg("l", "log", "log file", false, "", "string");
        cmd.add( logfile_arg );
        TCLAP::ValueArg<unsigned> verbosity_arg("v", "verbose", "level of verbosity [0 very low information, 1 much information]", false, 0, "number");
        cmd.add( verbosity_arg );
        TCLAP::SwitchArg pcs_arg("m", "modules", "list available modules", false);
        cmd.add( pcs_arg );
        cmd.parse( argc, argv ); // process can exit in this function

        // initialize
        ogsInit(argc, argv);

        // get parsed data
        // log file
        if (! logfile_arg.getValue().empty()) {
            if (!logog_file) delete logog_file;
            std::string log_file = logfile_arg.getValue();
            BaseLib::truncateFile(log_file); // do this not to append log into an existing file
            logog_file = new logog::LogFile(log_file.c_str());
            logog_file->SetFormatter( *custom_format );
        }

        SimulationInfo::outputHeader();
        // list modules
        const unsigned flag_list_modules (pcs_arg.getValue());
        if (flag_list_modules!=0) {
            ProcessBuilder::getInstance()->output();
        }

        const bool is_input_file_given = !input_arg.getValue().empty();
        if (is_input_file_given) {
            INFO("->Parsing input arguments");
            INFO("* project path     : %s", input_arg.getValue().c_str());
            if (! logfile_arg.getValue().empty()) {
                INFO("* log file path    : %s", logfile_arg.getValue().c_str());
            }

            // data output directory
            std::string output_dir_path = "";
            if (! output_dir_arg.getValue().empty()) {
                output_dir_path = output_dir_arg.getValue();
            }
            INFO("* output directory : %s", output_dir_path.c_str());

            if (! input_arg.getValue().empty()) {
                const std::string proj_path = input_arg.getValue();
                if (checkInputFiles(proj_path)) {
                    _sim_info = new SimulationInfo(proj_path, output_dir_path);
                } else {
                    ERR("***Error: Cannot find a project - %s", proj_path.c_str());
                }
            }
        }

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }

}

THMCSimulator::~THMCSimulator()
{
    BaseLib::releaseObject(_sim_info, _cpl_system);
    ogsExit();
}

bool THMCSimulator::checkInputFiles(const std::string& proj_path)
{
    // meanwhile OGS5 files are default
//    if(!BaseLib::IsFileExisting(proj_path+".pcs"))
//    {
//        ERR("Cannot find a PCS file - %s.pcs", proj_path.c_str());
//        return false;
//    }
    if(!BaseLib::IsFileExisting(proj_path+ ".pro"))
    {
        ERR("Cannot find a property file - %s.pro", proj_path.c_str());
        return false;
    }

    return true;
}

int THMCSimulator::execute()
{
    if (!_sim_info) return 0;

    BaseLib::Options opAll;
    Ogs6FemData* ogs6fem = Ogs6FemData::getInstance();
    const std::string proj_path = _sim_info->getProjectPath();
    ogs6fem->project_name = _sim_info->getProjectName();
    ogs6fem->output_dir = _sim_info->getOutputDirPath();

    //-------------------------------------------------------------------------
    // Read files
    //-------------------------------------------------------------------------
    INFO("->Reading a property file...");
    // coupling
    BaseLib::addXMLtoOptions(proj_path+".pro", opAll);
    BaseLib::Options* opOgs6 = opAll.getSubGroup("ogs6");
    if (opOgs6 == NULL) {
        ERR("***Error: tag <ogs6> was not found in the property file.");
        return 0;
    }
    // ogs5fem
    if (BaseLib::IsFileExisting(proj_path+".pcs")) {
        INFO("->Reading OGS5 input files...");
        ogs5::Ogs5FemData ogs5femdata;
        INFO("-------------------------------------------------");
        ogs5::Ogs5FemIO::read(proj_path, ogs5femdata);
        INFO("-------------------------------------------------");
        if (!Ogs5ToOgs6::convert(ogs5femdata, *ogs6fem, *opOgs6)) {
            ERR("***Error: Failure during conversion of ogs5 to ogs6.");
            return 0;
        }
        // output converted setting
        std::string str_conversion_logfile = ogs6fem->output_dir + "/converted_setting.log";
        std::ofstream of(str_conversion_logfile.c_str());
        opOgs6->printout(of);
        of.close();
    }


    // ddc

    //-------------------------------------------------------------------------
    // Setup simulation
    //-------------------------------------------------------------------------
    //INFO("Setting up simulation...");

    // construct mesh
    INFO("->Constructing meshes... %d mesh loaded", ogs6fem->list_mesh.size());
    for (size_t i=0; i<ogs6fem->list_mesh.size(); i++) {
        MeshLib::IMesh* msh = ogs6fem->list_mesh[i];
        msh->constructGeometricProperty();
        INFO("->mesh id %d: dim=%d, nodes=%d, elements=%d", i, msh->getDimension(), msh->getNumberOfNodes(), msh->getNumberOfElements());
    }

    // construct coupling system
    INFO("->Generating coupling system...");
    typedef class NumLib::TemplateCouplingStrucutreBuilder
        <
        NumLib::ITransientCoupledSystem,
        ProcessLib::Process,
        NumLib::AsyncPartitionedSystem,
        NumLib::TransientPartitionedAlgorithmFactory
        > CoupledProcessStrucutreBuilder;

    CoupledProcessStrucutreBuilder cpl_builder;
    if (_cpl_system!=NULL) delete _cpl_system;
    _cpl_system = cpl_builder.build(opOgs6, *GeoProcessBuilder::getInstance());
    std::vector<std::string> &list_mono_system_name = cpl_builder.getListOfMonolithicSystemName();
    if (list_mono_system_name.size()==0) {
        ERR("***Error: no active process is selected.");
        return 0;
    }
    if (!_cpl_system->check()) {
        ERR("***Error while checking coupled system");
        return 0;
    }

    // list up monolithic processes
    INFO("->Initializing all processes...");
    std::vector<ProcessLib::Process*> &list_mono_system = cpl_builder.getListOfMonolithicSystem();
    BaseLib::Options* opPcsParaList = opOgs6->getSubGroup("processParameters");
    for (size_t i=0; i<list_mono_system.size(); i++) {
        std::string &pcs_name = list_mono_system_name[i];
        ProcessLib::Process* pcs = list_mono_system[i];
        pcs->setProcessName(pcs_name);
        INFO("PCS %d: name=%s, type=%s (IN=%d, OUT=%d)", i, pcs_name.c_str(), pcs->getProcessType().c_str(), pcs->getNumberOfInputParameters(), pcs->getNumberOfOutputParameters());
        for (size_t j=0; j<pcs->getNumberOfInputParameters(); j++)
            INFO("* IN  %d: %s", j, pcs->getInputParameterName(j).c_str());
        for (size_t j=0; j<pcs->getNumberOfOutputParameters(); j++)
            INFO("* OUT %d: %s", j, pcs->getOutputParameterName(j).c_str());
        ogs6fem->list_pcs.insert(pcs_name, pcs);
        BaseLib::Options* opPCSList = opOgs6->getSubGroup("processList");
        BaseLib::Options* opPCS = NULL;
        if (opPCSList!=NULL) {
            std::vector<BaseLib::Options*> vec_opPCS = opPCSList->getSubGroupList("process");
            for (size_t i=0; i<vec_opPCS.size(); i++) {
                BaseLib::Options* opVal = vec_opPCS[i];
                if (opVal->getOption("name").compare(pcs_name)==0)
                    opPCS = opVal;
            }
        }
        if (opPCS==NULL) {
            INFO("* Process option not found.");
            opPCS = opPCSList->addSubGroup("process");
            opPCS->addOption("name", pcs_name);
            opPCS->addOption("type", pcs->getProcessType());
        }
        if (opPcsParaList) {
            BaseLib::Options* opPcsPara = opPcsParaList->getSubGroup(pcs_name);
            if (opPcsPara && opPcsPara->hasOption("TimeGroupID")) {
                size_t time_id = opPcsPara->getOptionAsNum<size_t>("TimeGroupID");
                opPCS->setOptionAsNum<size_t>("TimeGroupID", time_id);
            }
        }

        bool isPcsReady = pcs->initialize(opPCS!=NULL ? *opPCS : *opOgs6);
        if (!isPcsReady) {
            ERR("***Error while setting up processes");
            return 0;
        }
    }

    INFO("->Setting time stepping...");
    TimeSteppingControllerWithOutput timestepping(&ogs6fem->outController);
    timestepping.setTransientSystem(*_cpl_system);

    //TODO the following calculation should be done in TimeSteppingController
    double t_start = std::numeric_limits<double>::max();
    double t_end = -1 * std::numeric_limits<double>::max();
    if (ogs6fem->list_tim.size() > 0) {
        for (size_t i=0; i<ogs6fem->list_tim.size(); i++) {
            t_start = std::min(t_start, ogs6fem->list_tim[i]->getBeginning());
            t_end = std::max(t_end, ogs6fem->list_tim[i]->getEnd());
        }
    } else {
        INFO("Time step configuration not found.");
        t_start = 0.0;
        t_end = 1.0;
    }

    INFO("->Outputting the initial values...");
    ogs6fem->outController.outputData(NumLib::TimeStep(t_start));

    //-------------------------------------------------------------------------
    // Run simulation
    //-------------------------------------------------------------------------
    INFO("->Simulation is ready! start=%f, end=%f", t_start, t_end);
    BaseLib::RunTime runTime;
    runTime.start();

    timestepping.setBeginning(t_start); //TODO really need this? start, end is already given in timestep function
    size_t n_timesteps = timestepping.solve(t_end);

    runTime.stop();

    INFO("");
    INFO("->Simulation is finished.\n");
    INFO("#############################################################");
    INFO("*** Summary of this simulation");
    INFO("total time step : %d", n_timesteps);
    INFO("elapsed time   : %g sec", runTime.elapsed());
    INFO("#############################################################");
    INFO("");

    return 0;
}

} //end ogs6
