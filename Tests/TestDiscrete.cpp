/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestDiscrete.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>

#include <vector>

#include "BaseLib/CodingTools.h"
#include "BaseLib/BidirectionalMap.h"

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquation/DenseLinearEquation.h"
#include "MathLib/LinAlg/LinearEquation/SparseLinearEquation.h"
#ifdef USE_LIS
#include "MathLib/LinAlg/LinearEquation/LisLinearEquation.h"
#endif

#include "GeoLib/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"


#include "DiscreteLib/Core/IDiscreteLinearEquation.h"
#include "DiscreteLib/Core/IElemenetWiseLinearEquationLocalAssembler.h"
#include "DiscreteLib/Core/ElementWiseLinearEquationUpdater.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "DiscreteLib/Serial/DiscreteVector.h"
#include "DiscreteLib/Serial/SequentialElementWiseLinearEquationAssembler.h"
#include "DiscreteLib/SerialNodeDdc/SerialNodeDdcDistributedDiscreteSystem.h"
#include "DiscreteLib/SerialNodeDdc/SerialNodeDdcSharedDiscreteSystem.h"
#include "DiscreteLib/DDC/SequentialGlobaLocalMappingTable.h"
#include "DiscreteLib/DDC/RandomGlobalLocalMappingTable.h"
#include "DiscreteLib/DDC/SubDomain.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/Utils/SparsityBuilderFromNodeConnectivity.h"
#include "DiscreteLib/Utils/SparsityBuilderDDC.h"

#include "TestUtil.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace DiscreteLib;

struct DiscreteExample1
{
    std::vector<size_t> list_dirichlet_bc_id;
    std::vector<double> list_dirichlet_bc_value;
    static const size_t dim_eqs = 9;
    std::vector<double> exH;

    DiscreteExample1()
    {
        size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        list_dirichlet_bc_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        list_dirichlet_bc_value.resize(6);
        fill(list_dirichlet_bc_value.begin(), list_dirichlet_bc_value.begin()+3, .0);
        fill(list_dirichlet_bc_value.begin()+3, list_dirichlet_bc_value.end(), 1.0);
        exH.resize(9);
        for (size_t i=0; i<9; i++) {
            if (i%3==0) exH[i] = 1.0;
            if (i%3==1) exH[i] = 0.5;
            if (i%3==2) exH[i] = 0.;
        }
    }

    void setLocalDirichletBC(const BaseLib::BidirectionalMap<size_t, size_t> &map_global2localNodeId, std::vector<size_t> &local_dirichlet_bc_id, std::vector<double> &local_dirichlet_bc_value)
    {
        for (size_t i=0; i<list_dirichlet_bc_id.size(); i++) {
            if (map_global2localNodeId.countInA(list_dirichlet_bc_id[i])>0) {
                size_t local_id = map_global2localNodeId.mapAtoB(list_dirichlet_bc_id[i]);
                local_dirichlet_bc_id.push_back(local_id);
                local_dirichlet_bc_value.push_back(list_dirichlet_bc_value[i]);
            }
        }
    }


    class TestElementAssembler : public DiscreteLib::IElemenetWiseLinearEquationLocalAssembler
    {
        MathLib::LocalMatrix _m;
    public:
        TestElementAssembler()
        {
            _m.resize(4,4);
            _m(0,0) = 4.0; _m(0,1) = -1.0; _m(0,2) = -2.0; _m(0,3) = -1.0; 
            _m(1,1) = 4.0; _m(1,2) = -1.0; _m(1,3) = -2.0;
            _m(2,2) = 4.0; _m(2,3) = -1.0;
            _m(3,3) = 4.0;
            for (size_t i=0; i<4; i++)
                for (size_t j=0; j<i; j++) _m(i,j) = _m(j,i);
            _m *= 1.e-11/6.0;
        }
        void assembly(const MeshLib::IElement &/*e*/, MathLib::LocalEquation &eqs)
        {
            (*eqs.getA()) = _m;
        }
    };
};

DecomposedDomain* setupNDDC1()
{
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    SequentialGlobaLocalMappingTable* mapping = new SequentialGlobaLocalMappingTable(0, msh->getNumberOfNodes(), 0);
    SubDomain* dom1 = new SubDomain(msh, mapping);
    // global
    DecomposedDomain *ddc = new DecomposedDomain(DecompositionType::Node);
    ddc->addSubDomain(dom1);
    return ddc;
}

DecomposedDomain* setupNDDC2()
{
    MeshLib::IMesh *org_msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DecomposedDomain *ddc = new DecomposedDomain(DecompositionType::Node);
    {
        int dom1_eles[] = {0, 1, 2, 3};
        int dom1_ghost_nodes[] = {5, 6, 7, 8};
        std::set<size_t> list_ghost_nodes(dom1_ghost_nodes, dom1_ghost_nodes + 4);
        std::vector<size_t> dom1_e(dom1_eles, dom1_eles+4);
        BaseLib::BidirectionalMap<size_t, size_t> msh_node_id_mapping;
        MeshLib::IMesh* local_msh;
        MeshGenerator::generateSubMesh(*org_msh, dom1_e, local_msh, msh_node_id_mapping);
        SequentialGlobaLocalMappingTable* mapping = new SequentialGlobaLocalMappingTable(0, local_msh->getNumberOfNodes(), 0);
        SubDomain* dom = new SubDomain(local_msh, mapping, &list_ghost_nodes);
        ddc->addSubDomain(dom);
    }
    {
        MeshLib::IMesh* local_msh;
        int dom1_eles[] = {1, 2, 3};
        int dom1_ghost_nodes[] = {1, 2, 3, 4}; 
        size_t n_ghost = 4;
        std::vector<size_t> dom1_e(dom1_eles, dom1_eles+3);
        BaseLib::BidirectionalMap<size_t, size_t>* msh_node_id_mapping = new BaseLib::BidirectionalMap<size_t, size_t>();
        MeshGenerator::generateSubMesh(*org_msh, dom1_e, local_msh, *msh_node_id_mapping);
        std::set<size_t> list_ghost;
        for (size_t i=0; i<n_ghost; i++) {
            list_ghost.insert(msh_node_id_mapping->mapAtoB(dom1_ghost_nodes[i]));
        }
        RandomGlobalLocalMappingTable* mapping = new RandomGlobalLocalMappingTable(msh_node_id_mapping);
        SubDomain* dom = new SubDomain(local_msh, mapping, &list_ghost);
        ddc->addSubDomain(dom);
    }
    return ddc;
}

TEST(Discrete, NDDCSSVec1)
{
    DecomposedDomain* ddc = setupNDDC1();
    DecomposedMesh* msh = new DecomposedMesh(0, ddc);

    // discrete system
    SerialNodeDdcSharedDiscreteSystem dis(msh);
    dis.initialize(ddc);

    // vector
    IDiscreteVector<double>* v = dis.createVector<double>(ddc->getTotalNumberOfDecomposedObjects());
    for (size_t i=0; i<v->size(); ++i)
        (*v)[i] = i;
    ASSERT_EQ(9u, v->size());
    double expected[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_DOUBLE_ARRAY_EQ(expected, *v, 9);
}

TEST(Discrete, NDDCSSVec2)
{
    DecomposedDomain* ddc = setupNDDC2();
    DecomposedMesh* msh = new DecomposedMesh(0, ddc);
    // discrete system
    SerialNodeDdcSharedDiscreteSystem dis(msh);
    dis.initialize(ddc);

    // vector
    IDiscreteVector<double>* v = dis.createVector<double>(ddc->getTotalNumberOfDecomposedObjects());
    for (size_t i=0; i<v->size(); ++i)
        (*v)[i] = i;
    ASSERT_EQ(9u, v->size());
    double expected[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_DOUBLE_ARRAY_EQ(expected, *v, 9);
}

#ifdef USE_LIS
TEST(Discrete, NDDCSSEqs2)
{
    DiscreteExample1 ex1;
    DiscreteExample1::TestElementAssembler ele_assembler;
    LisLinearEquation lis;
    lis.getOption().ls_method = LIS_option::CG;
    lis.getOption().ls_precond = LIS_option::NONE;
    DecomposedDomain* ddc = setupNDDC2();
    DecomposedMesh* msh = new DecomposedMesh(0, ddc);

    // discrete system
    SerialNodeDdcSharedDiscreteSystem dis(msh);
    dis.initialize(ddc);
    // dof
    DofEquationIdTable dofManager;
    dofManager.addVariableDoFs(0, 0, ddc->getTotalNumberOfDecomposedObjects());
    dofManager.setNumberingType(DofNumberingType::BY_VARIABLE);
    dofManager.construct();
    // eqs
    IDiscreteLinearEquation* linear_eq = dis.createLinearEquation<LisLinearEquation, SparsityBuilderFromNodeConnectivity>(&lis, &dofManager);
    linear_eq->initialize();
    linear_eq->setPrescribedDoF(0, ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
    typedef DiscreteLib::ElementWiseLinearEquationUpdater<DiscreteExample1::TestElementAssembler,LisLinearEquation> MyUpdater;
    typedef SerialNodeDdcSharedDiscreteSystem::MyLinearEquationAssembler<MyUpdater,LisLinearEquation>::type MyGlobalAssembler;
    MyUpdater updater(msh, &ele_assembler);
    MyGlobalAssembler assembler(&updater);
    linear_eq->construct(assembler);
    linear_eq->solve();

    std::cout << "expected = ";
    for (size_t i=0; i<9; i++) std::cout << ex1.exH[i] << ", ";
    std::cout << "actual = ";
    for (size_t i=0; i<9; i++) std::cout << linear_eq->getLocalX()[i] << ", ";

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getLocalX(), 9, 1.e-5);
}

TEST(Discrete, NDDCSDEqs2)
{
    DiscreteExample1 ex1;
    DiscreteExample1::TestElementAssembler ele_assembler;
    LisLinearEquation lis;
    lis.getOption().ls_method = LIS_option::CG;
    lis.getOption().ls_precond = LIS_option::NONE;
    DecomposedDomain* ddc = setupNDDC2();
    DecomposedMesh* msh = new DecomposedMesh(0, ddc);
    
    // discrete system
    SerialNodeDdcDistributedDiscreteSystem dis(msh);
    dis.initialize(ddc);
    // dof
    DofEquationIdTable dofManager;
    dofManager.addVariableDoFs(0, 0, ddc->getTotalNumberOfDecomposedObjects());
    dofManager.setNumberingType(DofNumberingType::BY_VARIABLE);
    dofManager.construct();
    // eqs
    IDiscreteLinearEquation* linear_eq = dis.createLinearEquation<LisLinearEquation, SparsityBuilderFromNodeConnectivity>(&lis, &dofManager);
    linear_eq->initialize();
    linear_eq->setPrescribedDoF(0, ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
    typedef DiscreteLib::ElementWiseLinearEquationUpdater<DiscreteExample1::TestElementAssembler,LisLinearEquation> MyUpdater;
    typedef SerialNodeDdcSharedDiscreteSystem::MyLinearEquationAssembler<MyUpdater,LisLinearEquation>::type MyGlobalAssembler;
    MyUpdater updater(msh, &ele_assembler);
    MyGlobalAssembler assembler(&updater);
    linear_eq->construct(assembler);
    //linear_eq->getLinearEquation()->printout();
    linear_eq->solve();

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getLocalX(), 9, 1.e-5);
}
#endif

TEST(Discrete, VecSingle1)
{
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(msh);
    DiscreteVector<double> *v = dis.createVector<double>(msh->getNumberOfNodes());

    double expected[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    std::copy(expected, expected+9, v->begin());
    
    ASSERT_EQ(9u, v->size());
    ASSERT_DOUBLE_ARRAY_EQ(expected, *v, 9);
}


TEST(Discrete, VecDecomposed2)
{
#if 0
    MeshLib::IMesh* org_msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    const size_t n_dom = 2;
    std::vector<NodeDecomposedDiscreteSystem*> vec_dis(n_dom);
    DecomposedDiscreteVector<double> global_v(10);
    global_v.decompose(n_dom);

    omp_set_num_threads(n_dom);

#pragma omp parallel shared(global_v, org_msh, vec_dis, std::cout) default(none)
    {
        int iam = omp_get_thread_num();
        std::cout << "threads " << iam << std::endl;

        DecomposedDiscreteVectorLeaf<double>* local_v = global_v.createLocal(iam, iam*5, (iam+1)*5);
        for (size_t i=local_v->getRangeBegin(); i<local_v->getRangeEnd(); i++)
            (*local_v)[i] = iam*10 + i;
        Base::BidirectionalMap<size_t, size_t> msh_node_id_mapping;
        std::set<size_t> ghost_nodes;
        MeshLib::IMesh* local_msh;
        if (iam==0) {
            int dom1_eles[] = {0, 1, 2, 3};
            int dom1_ghost_nodes[] = {5, 6, 7, 8};
            std::vector<size_t> dom1_e(dom1_eles, dom1_eles+4);
            MeshGenerator::generateSubMesh(*org_msh, dom1_e, local_msh, msh_node_id_mapping);
            for (size_t i=0; i<4; i++) {
                ghost_nodes.insert(msh_node_id_mapping.mapAtoB(dom1_ghost_nodes[i]));
            }
        } else {
            int dom2_eles[] = {1, 2, 3};
            int dom2_ghost_nodes[] = {0, 1, 2, 3, 4};
            std::vector<size_t> dom2_e(dom2_eles, dom2_eles+3);
            MeshGenerator::generateSubMesh(*org_msh, dom2_e, local_msh, msh_node_id_mapping);
            for (size_t i=0; i<5; i++) {
                ghost_nodes.insert(msh_node_id_mapping.mapAtoB(dom2_ghost_nodes[i]));
            }
        }

        NodeDecomposedDiscreteSystem* dis = new NodeDecomposedDiscreteSystem(*local_msh, msh_node_id_mapping, ghost_nodes);
        vec_dis[iam] = dis;
        DecomposedDiscreteVector<double> *v = dis->createVector<double>(local_msh->getNumberOfNodes());
        for (size_t i=v->getRangeBegin(); i<v->getRangeEnd(); i++)
            (*v)[i] = iam*10 + i;
    }

    for (size_t i=0; i<global_v.size(); i++) {
        std::cout << global_v[i] << " ";
    }
    std::cout << std::endl;
#endif
}

#ifdef USE_LIS
TEST(Discrete, Lis1)
{
    DiscreteExample1 ex1;
    DiscreteExample1::TestElementAssembler ele_assembler;
    LisLinearEquation lis;
    lis.getOption().ls_method = LIS_option::CG;
    lis.getOption().ls_precond = LIS_option::NONE;
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);

    // define discrete system
    DiscreteSystem dis(msh);
    {
        // DoF?
        DofEquationIdTable dofManager;
        dofManager.addVariableDoFs(msh->getID(), 0, msh->getNumberOfNodes());
        dofManager.construct(DofNumberingType::BY_VARIABLE);
        // create a linear problem
        IDiscreteLinearEquation *linear_eq = dis.createLinearEquation<LisLinearEquation, SparsityBuilderFromNodeConnectivity>(&lis, &dofManager);
        // solve the equation
        linear_eq->initialize();
        linear_eq->setPrescribedDoF(0, ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
        typedef DiscreteLib::ElementWiseLinearEquationUpdater<DiscreteExample1::TestElementAssembler,LisLinearEquation> MyUpdater;
        typedef SerialNodeDdcSharedDiscreteSystem::MyLinearEquationAssembler<MyUpdater,LisLinearEquation>::type MyGlobalAssembler;
        MyUpdater updater(msh, &ele_assembler);
        MyGlobalAssembler assembler(&updater);
        linear_eq->construct(assembler);
        //linear_eq->getLinearEquation()->printout();
        linear_eq->solve();

        ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getLocalX(), 9, 1.e-5);
    }
}
#endif

//# OpenMP ###################################################################################################
#if 0
#ifdef _OPENMP
TEST(Discrete, OMP_vec1)
{
    //MeshLib::IMesh* org_msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    const size_t n_dom = 2;
    OMPGlobalDiscreteVector<double> global_v(10, n_dom);

    omp_set_num_threads(n_dom);

    #pragma omp parallel shared(global_v, /*org_msh,*/ std::cout) default(none)
    {
        int iam = omp_get_thread_num();
        //std::cout << "threads " << iam << std::endl;

        OMPLocalDiscreteVector<double>* local_v = global_v.createLocal(iam, iam*5, (iam+1)*5);
        for (size_t i=local_v->getRangeBegin(); i<local_v->getRangeEnd(); i++)
            local_v->global(i) = iam*10 + i;
    }

    double expected[] = {0, 1, 2, 3, 4, 15, 16, 17, 18, 19};
    ASSERT_DOUBLE_ARRAY_EQ(expected, global_v, 10);
}

#endif
#endif


#if 0
TEST(Discrete, OMP_eqs1)
{
    DiscreteExample1 ex1;
    struct NodeDDC
    {
        BaseLib::BidirectionalMap<size_t, size_t> map_global2localNodeId;
        std::set<size_t> ghost_nodes;
    };
    MeshLib::IMesh *org_msh = MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0);
    const size_t n_dom = 2;
    OmpNodeDdcDiscreteSystem global_dis(*org_msh, n_dom);

    omp_set_num_threads(n_dom);

    // define decomposed discrete system
    #pragma omp parallel shared(ex1, global_dis, org_msh, std::cout) default(none)
    {
        int iam = omp_get_thread_num();
        //std::cout << "threads " << iam << std::endl;

        MeshLib::IMesh* local_msh = 0;
        NodeDDC dom;
        if (iam==0) {
            // define a local problem with ddc
            int dom1_eles[] = {0, 1, 2, 3};
            int dom1_ghost_nodes[] = {5, 6, 7, 8};
            std::vector<size_t> dom1_e(dom1_eles, dom1_eles+4);

            MeshGenerator::generateSubMesh(*org_msh, dom1_e, local_msh, dom.map_global2localNodeId);
            for (size_t i=0; i<4; i++) {
                dom.ghost_nodes.insert(dom.map_global2localNodeId.mapAtoB(dom1_ghost_nodes[i]));
            }
        } else {
            int dom2_eles[] = {1, 2, 3};
            int dom2_ghost_nodes[] = {0, 1, 2, 3, 4};
            std::vector<size_t> dom2_e(dom2_eles, dom2_eles+3);
            MeshGenerator::generateSubMesh(*org_msh, dom2_e, local_msh, dom.map_global2localNodeId);
            for (size_t i=0; i<5; i++) {
                dom.ghost_nodes.insert(dom.map_global2localNodeId.mapAtoB(dom2_ghost_nodes[i]));
            }
        }
        std::vector<size_t> list_dirichlet_bc_id;
        std::vector<double> list_dirichlet_bc_value;
        ex1.setLocalDirichletBC(dom.map_global2localNodeId, list_dirichlet_bc_id, list_dirichlet_bc_value);

        // discrete system
        OmpNodeDdcLocalDiscreteSystem *local_dis = global_dis.createLocal(*local_msh, dom.map_global2localNodeId, dom.ghost_nodes);

        //// create dof map
        //DofMapManager dofManager;
        //size_t dofId = dofManager.addDoF(local_msh->getNumberOfNodes(), dom.ghost_nodes);
        //dofManager.construct(DofMapManager::BY_POINT);

        // create a linear problem
        LisLinearEquation lis;
        lis.getOption().ls_method = LIS_option::CG;
        lis.getOption().ls_precond = LIS_option::NONE;
        DofEquationIdTable dofManager;
        //size_t dofId = dofManager.addVariableDoF(org_msh->getNumberOfNodes());
        //dofManager.construct(DofNumberingType::BY_POINT);
        //IDiscreteLinearEquation *linear_eq = local_dis->createLinearEquation<CRSLisSolver, SparsityBuilderFromNodeConnectivityWithInactiveDoFs>(lis, dofManager);
        //linear_eq->setPrescribedDoF(dofId, list_dirichlet_bc_id, list_dirichlet_bc_value);
        //// construct and solve
        //DiscreteExample1::TestElementAssembler ele_assembler;
        //linear_eq->construct(ElementBasedAssembler(ele_assembler));    
    }


     IDiscreteLinearEquation *global_eq;
     global_eq->solve();
}

#endif
