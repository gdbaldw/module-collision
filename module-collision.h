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

typedef std::pair<Vec3, Vec3> PointPair;
enum State {UNDEFINED, SLIP, NO_SLIP};

class Collision : public RodWithOffset {
private:
    typedef RodWithOffset super;
    const BasicScalarFunction* pSF;
    const doublereal penetration;
    integer iFirstDofIndex;
    integer iR;
    integer iC;
    int iNumRowsNode;
    int iNumColsNode;
    int iNumDofs;
    doublereal dCalcEpsilon(void);
    Vec3 Arm2;
    Vec3 F_reaction;
    int index;
    std::vector<Vec3> Arm2_vector;
    std::vector<Vec3> f1_vector;
    std::vector<Vec3> f2_vector;
    std::vector<State> state_vector;
    std::vector<Vec3> Ft_vector;
public:
    Collision(const DofOwner* pDO, 
        const ConstitutiveLaw1D* pCL, const BasicScalarFunction* pSFTmp, const doublereal penetrationTmp,
        const StructNode* pN1, const StructNode* pN2);
    unsigned int iGetNumDof(void) const;
    void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
    void PushBackContact(const btVector3& point1, const btVector3& point2);
    void ClearContacts(void);
    void ResetStates(void);
    int ContactsSize(void);
    bool SetOffsets(const int i);
    std::ostream& OutputAppend(std::ostream& out) const;
    void SetIndices(integer* iIndex, integer* iRow, integer* iCol);
    VariableSubMatrixHandler&
    AssJac(VariableSubMatrixHandler& WorkMat,
        doublereal dCoef,
        const VectorHandler& XCurr, 
        const VectorHandler& XPrimeCurr);
    void AssMat(FullSubMatrixHandler& WM,
        doublereal dCoef,
        const VectorHandler& XCurr,
        const int i);

    SubVectorHandler& 
    AssRes(SubVectorHandler& WorkVec,
        doublereal dCoef,
        const VectorHandler& XCurr, 
        const VectorHandler& XPrimeCurr);
    void AssVec(SubVectorHandler& WorkVec,
        doublereal dCoef,
        const VectorHandler& XCurr,
        const int i);
};

class CollisionWorld
: virtual public Elem, public UserDefinedElem {
private:
    integer iNumDofs;
    integer iNumRows;
    integer iNumCols;
    btDefaultCollisionConfiguration* configuration;
    btCollisionDispatcher* dispatcher;
    btBroadphaseInterface* broadphase;
    btCollisionWorld* world;
    typedef std::pair<btCollisionObject*, btCollisionObject*> ObjectPair;
    std::map<ObjectPair, Collision*> objectpair_collision_map;
    std::set<const Node*> nodes;
public:
    CollisionWorld(unsigned uLabel, const DofOwner *pDO,
        DataManager* pDM, MBDynParser& HP);
    ~CollisionWorld(void);
    unsigned int iGetNumDof(void) const;
    DofOrder::Order GetEqType(unsigned int i) const;
    DofOrder::Order GetDofType(unsigned int i) const;
    void Output(OutputHandler& OH) const;
    void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
    unsigned int iGetNumPrivData(void) const;
    unsigned int iGetPrivDataIdx(const char *s) const;
    doublereal dGetPrivData(unsigned int i) const;
    int iGetNumConnectedNodes(void) const;
    void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
    std::ostream& Restart(std::ostream& out) const;
    unsigned int iGetInitialNumDof(void) const;
    void ResetStates(void);
    void ClearAndPushContacts(void);

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
    Mat3x3 Rh;
    btCollisionObject* ob;
    btCollisionShape* shape;
public:
    CollisionObject(unsigned uLabel, const DofOwner *pDO,
        DataManager* pDM, MBDynParser& HP);
    ~CollisionObject(void);

    void Output(OutputHandler& OH) const;
    void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
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
    unsigned int iGetNumPrivData(void) const;
    unsigned int iGetPrivDataIdx(const char *s) const;
    doublereal dGetPrivData(unsigned int i) const;
    int iGetNumConnectedNodes(void) const;
    void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
    void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
        SimulationEntity::Hints *ph);
    std::ostream& Restart(std::ostream& out) const;
    unsigned int iGetInitialNumDof(void) const;
    void 
    InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
       VariableSubMatrixHandler&
    InitialAssJac(VariableSubMatrixHandler& WorkMat, 
              const VectorHandler& XCurr);
       SubVectorHandler& 
    InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
};

#endif
