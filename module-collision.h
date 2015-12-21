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

#include <boost/foreach.hpp>
#include <fcl/shape/geometric_shapes.h>
#include <fcl/shape/geometric_shapes_utility.h>
#include <fcl/narrowphase/narrowphase.h>
#include <fcl/broadphase/broadphase.h>
#include <fcl/collision.h>

class CollisionObjectData {
public:
    CollisionObjectData(const StructNode* pNode, btCollisionObject* pBulletObject, fcl::CollisionObject* pFCLObject, std::string material);
    ~CollisionObjectData(void);
    const StructNode* pNode;
    btCollisionObject* pBulletObject;
    fcl::CollisionObject* pFCLObject;
    std::string material;
};

class Collision :
public ConstitutiveLaw1DOwner  {
private:
    typedef RodWithOffset super;
	const StructDispNode* pNode1;
	const StructDispNode* pNode2;
    const BasicScalarFunction* pSF;
    const doublereal penetration;
    integer iR;
    integer iC;
    int iNumRowsNode;
    int iNumColsNode;
    std::vector<Vec3> Arm2;
    std::vector<Vec3> f1;
    std::vector<Vec3> f2;
    std::vector<doublereal> Fn_Norm;
    std::vector<Vec3> Ft;
    std::vector<Vec3> v;
    std::vector<doublereal> dEpsilonPrime;
    void AssMat(FullSubMatrixHandler& WM, doublereal dCoef, const int i);
    void AssVec(SubVectorHandler& WorkVec, doublereal dCoef, const int i);
public:
    Collision(const DofOwner* pDO, 
        const ConstitutiveLaw1D* pCL, const BasicScalarFunction* pSFTmp, const doublereal penetrationTmp,
        const StructNode* pN1, const StructNode* pN2, integer* iRow, integer* iCol);
    fcl::CollisionRequest request;
    fcl::CollisionResult result;
    bool done;
    void PushBackContact(const btVector3& point1, const btVector3& point2);
    void ClearContacts(void);
    int ContactsSize(void);
    bool SetOffsets(const int i);
    std::ostream& OutputAppend(std::ostream& out, int i) const;

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

typedef std::pair<btCollisionObject*, btCollisionObject*> ObjectPair;
typedef std::pair<fcl::CollisionObject*, fcl::CollisionObject*> FCLObjectPair;

class CollisionWorld
: virtual public Elem, public UserDefinedElem {
private:
    integer iNumRows;
    integer iNumCols;
    btDefaultCollisionConfiguration* configuration;
    btCollisionDispatcher* dispatcher;
    btBroadphaseInterface* broadphase;
    btCollisionWorld* world;
    fcl::BroadPhaseCollisionManager* fcl_broadphase;
    std::map<const ObjectPair, Collision*> objectpair_collision_map;
    std::set<const Node*> nodes;
    std::ostringstream ss;
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
    void ClearAndPushContacts(void);

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
    Mat3x3 Rh;
    btCollisionObject* ob;
    btCollisionShape* shape;
    typedef boost::shared_ptr <fcl::CollisionGeometry> CollisionGeometryPtr_t;
    fcl::CollisionObject* fcl_ob;
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
