/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DofEquationIdTable.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <map>
#include <cassert>
#include <algorithm>

#include "BaseLib/CodingTools.h"

#include "IEquationIdStorage.h"
#include "SequentialEquationIdStorage.h"
#include "RandomEquationIdStorage.h"


namespace DiscreteLib
{

/**
 * \brief Type of numbering DoF's equation index
 */
struct DofNumberingType
{
    enum type
    {
        BY_VARIABLE,
        BY_POINT
    };
};


/**
 * \brief DoF and equation index table
 *
 * This class contains a table which relates DoFs and index in a equation.
 * A DoF is implicitly defined as a unique combination of the followings 
 * - variable type
 * - mesh id
 * - discrete point id in the mesh (e.g. node id, element id)
 * 
 */
class DofEquationIdTable
{
public:

    ///
    DofEquationIdTable()
    : _total_dofs(0), _total_dofs_without_ghosts(0), _is_seq_address(false),
      _numbering_type(DofNumberingType::BY_POINT), _local_numbering_type(DofNumberingType::BY_VARIABLE)
    {};

    ///
    virtual ~DofEquationIdTable()
    {
        for (size_t i=0; i<_map_var2dof.size(); i++) {
            BaseLib::releaseObjectsInStdMap(_map_var2dof[i]);
        }
        _map_var2dof.clear();
    }

    //---------------------------------------------------------------------------------------------
    // setup DoFs
    //---------------------------------------------------------------------------------------------
    /// add DoFs
    /// @param var_id       variable type
    /// @param mesh_id      mesh id
    /// @param list_dof_pt_id   a list of discrete point id defined as DoFs
    size_t addVariableDoFs(size_t mesh_id, const std::vector<size_t> &list_dof_pt_id)
    {
        RandomEquationIdStorage* address = new RandomEquationIdStorage(list_dof_pt_id);
        size_t var_id = _map_var2dof.size();
        _map_var2dof.resize(_map_var2dof.size()+1);
        _map_var2dof[var_id][mesh_id] = address;
        return var_id;
    }

    /// add DoFs
    /// @param var_id           variable type
    /// @param mesh_id          mesh id
    /// @param dof_pt_id_begin  the beginning id of discrete points
    /// @param dof_pt_count     the number of discrete points
    size_t addVariableDoFs(size_t mesh_id, size_t dof_pt_id_begin, size_t dof_pt_count)
    {
        RandomEquationIdStorage* address = new RandomEquationIdStorage(dof_pt_id_begin, dof_pt_count);
        size_t var_id = _map_var2dof.size();
        _map_var2dof.resize(_map_var2dof.size()+1);
        _map_var2dof[var_id][mesh_id] = address;
        _map_msh2var[mesh_id].push_back(var_id);
        return var_id;
    }

    /// add DoFs
    void addVariableDoFs(size_t n_var, size_t mesh_id, size_t dof_pt_id_begin, size_t dof_pt_count)
    {
        _is_seq_address = true;
        for (size_t i=0; i<n_var; i++) {
            SequentialEquationIdStorage* address = new SequentialEquationIdStorage(dof_pt_id_begin, dof_pt_count);
            size_t var_id = _map_var2dof.size();
            _map_var2dof.resize(_map_var2dof.size()+1);
            _map_var2dof[var_id][mesh_id] = address;
            _map_msh2var[mesh_id].push_back(var_id);
        }
    }

    /// deactivate DoFs
    void deactivateDoFs(size_t var_id, size_t mesh_id, size_t pt_id)
    {
        getPointEquationIdTable(var_id, mesh_id)->activate(pt_id, false);
    }

    /// set ghost points for overlapped regions appeared with domain-decomposition methods
    void setGhostPoints(size_t mesh_id, std::vector<size_t> &list_pt_id)
    {
        _ghost_pt[mesh_id].assign(list_pt_id.begin(), list_pt_id.end());
    }

    //---------------------------------------------------------------------------------------------
    // construct the table
    //---------------------------------------------------------------------------------------------
    /// numbering DoFs
    virtual void construct(long offset=0);

    //---------------------------------------------------------------------------------------------
    // get summary of the table
    //---------------------------------------------------------------------------------------------
    //# Mesh
    /// get the number of meshes
    size_t getNumberOfMeshes() const {return _map_msh2var.size();};

    //# Variables
    /// get the number of registered variables
    size_t getNumberOfVariables() const { return _map_var2dof.size(); }
    /// get the number of variables associated with the given mesh
    size_t getNumberOfVariables(size_t mesh_id) {return _map_msh2var[mesh_id].size();};

    //# Ghost 
    /// is this point a ghost?
    bool isGhostPoint(size_t mesh_id, size_t pt_id) const
    {
        if (_ghost_pt.size()==0) return false;
        assert(_ghost_pt.count(mesh_id)>0);
        const std::vector<size_t>& ghosts = _ghost_pt.find(mesh_id)->second;
        return std::find(ghosts.begin(), ghosts.end(), pt_id)!=ghosts.end();
    }
    /// get the number of ghost points
    size_t getNumberOfGhostPoints(size_t mesh_id) {return _ghost_pt[mesh_id].size();}
    /// get the total number of ghost points
    size_t getTotalNumberOfGhostPoints() const
    {
        size_t count = 0;
        for (std::map<size_t, std::vector<size_t> >::const_iterator itr=_ghost_pt.begin(); itr!=_ghost_pt.end(); ++itr) {
            count += itr->second.size();
        }
        return count;
    }

    //# DoFs
    /// is this DoF active?
    bool isActiveDoF(size_t var_id, size_t mesh_id, size_t pt_id) const
    {
        return getPointEquationIdTable(var_id, mesh_id)->isActive(pt_id);
    }
    /// get the total number of active DoFs
    size_t getTotalNumberOfActiveDoFs() const { return _total_dofs; }
    /// get the total number of active DoFs excluding ghost points
    size_t getTotalNumberOfActiveDoFsWithoutGhost() const {return _total_dofs_without_ghosts;};


    //---------------------------------------------------------------------------------------------
    // equation address mapping
    //---------------------------------------------------------------------------------------------
    /**
     * return equation id
     *
     * @param var_id
     * @param mesh_id
     * @param pt_id
     * @return equation id or npos if no match
     */
    size_t mapEqsID(size_t var_id, size_t mesh_id, size_t pt_id) const
    {
        const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
        return add->address(pt_id);
    }

    /**
     * return a list of equation id
     *
     * @param var_id
     * @param mesh_id
     * @param pt_id
     * @return equation id or npos if no match
     */
    void mapEqsID(size_t var_id, size_t mesh_id, const std::vector<size_t> &pt_id, std::vector<size_t> &eqs_id) const
    {
        eqs_id.resize(pt_id.size());
        const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
        for (size_t i=0; i<pt_id.size(); i++)
            eqs_id[i] = add->address(pt_id[i]);
    }

    /**
     * return a list of equation id for all variables
     *
     * @param mesh_id
     * @param pt_id
     * @param eqs_id
     */
    void mapEqsID(size_t mesh_id, const std::vector<size_t> &pt_id, std::vector<size_t> &eqs_id) const
    {
        if (_map_msh2var.count(mesh_id)==0) return;

        const std::vector<size_t> &list_var = _map_msh2var.find(mesh_id)->second;
        const size_t n_pt = pt_id.size();
        const size_t n_dof_per_pt = list_var.size();
        const size_t n_total_dof = n_pt * n_dof_per_pt;
        eqs_id.resize(n_total_dof);

        if (_local_numbering_type==DofNumberingType::BY_VARIABLE) {
            for (size_t i=0; i<list_var.size(); i++) {
                size_t var_id = list_var[i];
                const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
                for (size_t j=0; j<pt_id.size(); j++) {
                    eqs_id[j+i*n_pt] = add->address(pt_id[j]);
                }
            }
        } else {
            for (size_t i=0; i<list_var.size(); i++) {
                size_t var_id = list_var[i];
                const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
                for (size_t j=0; j<pt_id.size(); j++) {
                    eqs_id[i + j*n_dof_per_pt] = add->address(pt_id[j]);
                }
            }
        }
    }

    /**
     * return a list of equation id for all variables
     *
     * Length of a created list is the number of points * the number of variables.
     * a created list may contain entries with -1 if no valid address is found
     * for corresponding key (var id, point id). In addition to that, the second
     * may have -1 for ghost points.
     *
     * @param mesh_id
     * @param pt_id
     * @param eqs_id
     * @param eqs_id_without_ghost
     */
    void mapEqsID(
            size_t mesh_id, const std::vector<size_t> &list_pt_id,
            std::vector<size_t> &list_eqs_id,
            std::vector<size_t> &list_eqs_id_without_ghost
            ) const
    {
        if (_map_msh2var.count(mesh_id)==0) return;

        const std::vector<size_t> &list_var = _map_msh2var.find(mesh_id)->second;
        const size_t n_pt = list_pt_id.size();
        const size_t n_dof_per_pt = list_var.size();
        const size_t n_total_dof = n_pt * n_dof_per_pt;
        list_eqs_id.resize(n_total_dof);
        list_eqs_id_without_ghost.resize(n_total_dof);

        for (size_t i_var=0; i_var<list_var.size(); i_var++) {
            const size_t var_id = list_var[i_var];
            const IEquationIdStorage* table = getPointEquationIdTable(var_id, mesh_id);
            for (size_t i_pt=0; i_pt<n_pt; i_pt++) {
                const size_t pt_id = list_pt_id[i_pt];
                size_t pos;
                if (_local_numbering_type==DofNumberingType::BY_VARIABLE) {
                    pos = i_pt + i_var * n_pt;
                } else {
                    pos = i_var + i_pt * n_dof_per_pt;
                }
                size_t this_eqs_id = table->address(pt_id);
                list_eqs_id[pos] = isGhostPoint(mesh_id, pt_id) ? BaseLib::index_npos : this_eqs_id;
                list_eqs_id_without_ghost[pos] = this_eqs_id;
            }
        }
    }
    
    /**
     * return a list of equation id for all variables
     *
     * Length of a created list is sum of the number of active points for each variable.
     * The second list may have -1 for ghost points.
     *
     * @param mesh_id
     * @param pt_id
     * @param eqs_id
     * @param eqs_id_without_ghost
     */
    void mapEqsIDreduced(
            size_t mesh_id, const std::vector<size_t> &list_pt_id,
            std::vector<size_t> &list_eqs_id,
            std::vector<size_t> &list_eqs_id_without_ghost
            ) const
    {
        if (_map_msh2var.count(mesh_id)==0) return;

        const std::vector<size_t> &list_var = _map_msh2var.find(mesh_id)->second;
        const size_t n_pt = list_pt_id.size();
        const size_t n_dof_per_pt = list_var.size();
        const size_t n_total_dof = n_pt * n_dof_per_pt;
        list_eqs_id.reserve(n_total_dof);
        list_eqs_id_without_ghost.reserve(n_total_dof);

        if (_local_numbering_type==DofNumberingType::BY_VARIABLE) {
            for (size_t i_var=0; i_var<list_var.size(); i_var++) {
                const size_t var_id = list_var[i_var];
                const IEquationIdStorage* table = getPointEquationIdTable(var_id, mesh_id);
                std::vector<size_t> active_pt_list;
                std::vector<size_t> active_eqs_id_list;
                for (size_t i_pt=0; i_pt<n_pt; i_pt++) {
                    const size_t pt_id = list_pt_id[i_pt];
                    if (!table->hasKey(pt_id)) continue;
                    size_t this_eqs_id = table->address(pt_id);
                    list_eqs_id.push_back(isGhostPoint(mesh_id, pt_id) ? BaseLib::index_npos : this_eqs_id);
                    list_eqs_id_without_ghost.push_back(this_eqs_id);
                }
            }
        } else { // by point
            std::vector<std::vector<size_t> > var_active_pt_list(list_var.size());
            std::vector<std::vector<size_t> > var_active_eqs_id_list(list_var.size());
            for (size_t i_pt=0; i_pt<n_pt; i_pt++) {
                const size_t pt_id = list_pt_id[i_pt];
                for (size_t i_var=0; i_var<list_var.size(); i_var++) {
                    const size_t var_id = list_var[i_var];
                    const IEquationIdStorage* table = getPointEquationIdTable(var_id, mesh_id);
                    if (!table->hasKey(pt_id)) continue;
                    size_t this_eqs_id = table->address(pt_id);
                    list_eqs_id.push_back(isGhostPoint(mesh_id, pt_id) ? BaseLib::index_npos : this_eqs_id);
                    list_eqs_id_without_ghost.push_back(this_eqs_id);
                }
            }
        }
    }

    /**
     *
     * @param mesh_id
     * @param pt_id
     * @param local_table
     */
    void createLocalMappingTable(
        size_t mesh_id, const std::vector<size_t> &list_pt_id,
        DofEquationIdTable &local_table
        ) const
    {
        const std::vector<size_t> &list_var = _map_msh2var.find(mesh_id)->second;
        const size_t n_pt = list_pt_id.size();
        
        if (_local_numbering_type==DofNumberingType::BY_VARIABLE) {
            size_t active_entry_counter = 0;
            for (size_t i_var=0; i_var<list_var.size(); i_var++) {
                const size_t var_id = list_var[i_var];
                const IEquationIdStorage* table = getPointEquationIdTable(var_id, mesh_id);
                std::vector<size_t> active_pt_list;
                std::vector<size_t> active_eqs_id_list;
                for (size_t i_pt=0; i_pt<n_pt; i_pt++) {
                    const size_t pt_id = list_pt_id[i_pt];
                    if (!table->hasKey(pt_id)) continue;
                    active_pt_list.push_back(pt_id);
                    active_eqs_id_list.push_back(active_entry_counter++);
                }
                local_table.addVariableDoFs(mesh_id, active_pt_list);
                IEquationIdStorage *local_storage = local_table.getPointEquationIdTable(i_var, mesh_id);
                for (size_t j=0; j<active_pt_list.size(); j++) {
                    local_storage->set(active_pt_list[j], active_eqs_id_list[j]);
                }
            }
        } else { // by point
            size_t active_entry_counter = 0;
            std::vector<std::vector<size_t> > var_active_pt_list(list_var.size());
            std::vector<std::vector<size_t> > var_active_eqs_id_list(list_var.size());
            for (size_t i_pt=0; i_pt<n_pt; i_pt++) {
                const size_t pt_id = list_pt_id[i_pt];
                for (size_t i_var=0; i_var<list_var.size(); i_var++) {
                    const size_t var_id = list_var[i_var];
                    const IEquationIdStorage* table = getPointEquationIdTable(var_id, mesh_id);
                    if (!table->hasKey(pt_id)) continue;
                    var_active_pt_list[i_var].push_back(pt_id);
                    var_active_eqs_id_list[i_var].push_back(active_entry_counter++);
                }
            }

            for (size_t i_var=0; i_var<var_active_pt_list.size(); i_var++) {
                local_table.addVariableDoFs(mesh_id, var_active_pt_list[i_var]);
                IEquationIdStorage *local_storage = local_table.getPointEquationIdTable(i_var, mesh_id);
                for (size_t j=0; j<var_active_pt_list[i_var].size(); j++) {
                    local_storage->set(var_active_pt_list[i_var][j], var_active_eqs_id_list[i_var][j]);
                }
            }
        }

        local_table.setNumberingType(_local_numbering_type);
        //local_table.construct();
    }
    /**
     *
     * @param eqs_id
     * @param var_id
     * @param mesh_id
     * @param pt_id
     */
    void mapDoF(size_t eqs_id, size_t &var_id, size_t &mesh_id, size_t &pt_id) const
    {
        var_id = BaseLib::index_npos;
        mesh_id = BaseLib::index_npos;
        pt_id = BaseLib::index_npos;
        for (size_t i=0; i<_map_var2dof.size(); i++) {
            const std::map<size_t, IEquationIdStorage*> &obj = _map_var2dof[i];
            for (std::map<size_t, IEquationIdStorage*>::const_iterator itr=obj.begin(); itr!=obj.end(); ++itr) {
                const IEquationIdStorage* pt2dof = itr->second;
                if (!pt2dof->hasValue(eqs_id)) continue;
                pt_id = pt2dof->key(eqs_id);
                var_id = i;
                mesh_id = itr->first;
                return;
            }
        }
    }

    /**
     *
     * @param var_id
     * @param mesh_id
     * @return
     */
    const IEquationIdStorage* getPointEquationIdTable(size_t var_id, size_t mesh_id) const 
    {
        std::map<size_t, IEquationIdStorage*>::const_iterator itr = _map_var2dof[var_id].find(mesh_id);
        if (itr!=_map_var2dof[var_id].end()) {
            return itr->second;
        } else {
            return 0;
        }
    }

    /**
     *
     * @param var_id
     * @param mesh_id
     * @return
     */
    IEquationIdStorage* getPointEquationIdTable(size_t var_id, size_t mesh_id) { return _map_var2dof[var_id][mesh_id]; }

    /**
     *
     * @param numbering
     */
    void setNumberingType(DofNumberingType::type numbering) {_numbering_type = numbering;};

    /**
     *
     * @return
     */
    DofNumberingType::type getNumberingType() const {return _numbering_type;};

    /**
     *
     * @param numbering
     */
    void setLocalNumberingType(DofNumberingType::type numbering) {_local_numbering_type = numbering;};

    /**
     *
     * @return
     */
    DofNumberingType::type getLocalNumberingType() const {return _local_numbering_type;};

protected:
    void setTotalDoFs(size_t n) {_total_dofs = n;};
    //std::map<size_t, std::vector<DofMap*> >& getMapMsh2Dof() {return _map_msh2dof;}
    /// get the DoF mapping for the given variable

private:
    DISALLOW_COPY_AND_ASSIGN(DofEquationIdTable);


private:
    std::vector<std::map<size_t, IEquationIdStorage*> > _map_var2dof; // var -> (mesh, address)
    //std::vector<Base::IMappedAddress*> _map_var2dof;
    std::map<size_t, std::vector<size_t> > _map_msh2var; // mesh -> (var)
    std::map<size_t, std::vector<size_t> > _ghost_pt; // mesh -> ghost points

    size_t _total_dofs;
    size_t _total_dofs_without_ghosts;
    bool _is_seq_address;

    DofNumberingType::type _numbering_type;
    DofNumberingType::type _local_numbering_type;
};

} //end
