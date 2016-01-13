/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * module-collision
 * AUTHOR: G. Douglas Baldwin
        Copyright (C) 2015 all rights reserved.
 */

#ifndef MODULE_COLLISION_H
#define MODULE_COLLISION_H

#include "intersect.h"

class Contact {
public:
    Contact(std::pair<fcl::Vec3f, fcl::Vec3f> pt_pair, const StructDispNode* pNode1, const StructDispNode* pNode2, doublereal penetration_ratio);
    ~Contact(void);
    Vec3 Arm1;
    Vec3 f1;
    Vec3 f2;
    Vec3 Ft;
    doublereal Fn_Norm;
    Vec3 tangent;
};

class CollisionObjectData {
public:
    CollisionObjectData(const StructNode* pNode, fcl::CollisionObject* pObject, std::string material);
    ~CollisionObjectData(void);
    const StructNode* pNode;
    fcl::CollisionObject* pObject;
    std::string material;
};

class Collision :
public ConstitutiveLaw1DOwner {
private:
	const StructDispNode* pNode1;
	const StructDispNode* pNode2;
    fcl::CollisionObject* pObject1;
    fcl::CollisionObject* pObject2;
    const BasicScalarFunction* pSF;
    const doublereal penetration_ratio;
    integer iR;
    integer iC;
    int iNumRowsNode;
    int iNumColsNode;
    std::vector<doublereal> dEpsilonPrime;
    std::vector<Contact> contacts;
    void AssMat(FullSubMatrixHandler& WM, doublereal dCoef, Contact& contact);
    void AssVec(SubVectorHandler& WorkVec, doublereal dCoef, Contact& contact);
    FCL::Func func;
    std::vector<Vec3> tangents;
public:
    Collision(FCL::Func func,
        const ConstitutiveLaw1D* pCL, const BasicScalarFunction* pSFTmp, const doublereal penetration_ratio,
        const CollisionObjectData* pD1, const CollisionObjectData* pD2, integer* iRow, integer* iCol);
    void Intersect(void);
    void ClearContacts(void);
    void ClearAndSetTangents(void);
    std::ostream& OutputAppend(std::ostream& out) const;

    VariableSubMatrixHandler&
    AssJac(VariableSubMatrixHandler& WorkMat,
        doublereal dCoef,
        const VectorHandler& XCurr, 
        const VectorHandler& XPrimeCurr);

    SubVectorHandler& 
    AssRes(SubVectorHandler& WorkVec,
        doublereal dCoef,
        const VectorHandler& XCurr, 
        const VectorHandler& XPrimeCurr);
};

class CollisionWorld
: virtual public Elem, public UserDefinedElem {
private:
    integer iNumRows;
    integer iNumCols;
    fcl::BroadPhaseCollisionManager* collision_manager;
    std::map<const FCL::ObjectPair, Collision*> objectpair_collision_map;
    std::set<const Node*> nodes;
    std::ostringstream ss;
    FCL::FuncMatrix func_matrix;
public:
    CollisionWorld(unsigned uLabel, const DofOwner *pDO,
        DataManager* pDM, MBDynParser& HP);
    ~CollisionWorld(void);
    void Output(OutputHandler& OH) const;
    void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
    unsigned int iGetNumPrivData(void) const;
    unsigned int iGetPrivDataIdx(const char *s) const;
    doublereal dGetPrivData(unsigned int i) const;
    int iGetNumConnectedNodes(void) const;
    void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
    std::ostream& Restart(std::ostream& out) const;
    unsigned int iGetInitialNumDof(void) const;

    void
    AfterPredict(VectorHandler& X, VectorHandler& XP);
    
    void
    AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

    VariableSubMatrixHandler& 
    AssJac(VariableSubMatrixHandler& WorkMat,
        doublereal dCoef, 
        const VectorHandler& XCurr,
        const VectorHandler& XPrimeCurr);

    SubVectorHandler& 
    AssRes(SubVectorHandler& WorkVec,
        doublereal dCoef,
        const VectorHandler& XCurr, 
        const VectorHandler& XPrimeCurr);

    void
    SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
        SimulationEntity::Hints *ph);

    void 
    InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

    VariableSubMatrixHandler&
    InitialAssJac(VariableSubMatrixHandler& WorkMat, const VectorHandler& XCurr);

    SubVectorHandler& 
    InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
};

class CollisionObject
: virtual public Elem, public UserDefinedElem {
private:
    const StructNode* pNode;
    Vec3 f;
    Mat3x3 R;
    fcl::CollisionObject* ob;
public:
    CollisionObject(unsigned uLabel, const DofOwner *pDO,
        DataManager* pDM, MBDynParser& HP);
    ~CollisionObject(void);
    void Output(OutputHandler& OH) const;
    void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
    unsigned int iGetNumPrivData(void) const;
    unsigned int iGetPrivDataIdx(const char *s) const;
    doublereal dGetPrivData(unsigned int i) const;
    int iGetNumConnectedNodes(void) const;
    void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
    std::ostream& Restart(std::ostream& out) const;
    unsigned int iGetInitialNumDof(void) const;
    void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

    VariableSubMatrixHandler& 
    AssJac(VariableSubMatrixHandler& WorkMat,
        doublereal dCoef, 
        const VectorHandler& XCurr,
        const VectorHandler& XPrimeCurr);

    SubVectorHandler& 
    AssRes(SubVectorHandler& WorkVec,
        doublereal dCoef,
        const VectorHandler& XCurr, 
        const VectorHandler& XPrimeCurr);

    void
    SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
        SimulationEntity::Hints *ph);

    VariableSubMatrixHandler&
    InitialAssJac(VariableSubMatrixHandler& WorkMat, const VectorHandler& XCurr);

    SubVectorHandler& 
    InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
};

#endif
