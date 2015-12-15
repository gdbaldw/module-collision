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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <ostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"
#include "btBulletDynamicsCommon.h"
#include <set>
#include "rodj.h"
#include <limits>
#include "module-collision.h"

Collision::Collision(const DofOwner* pDO, 
    const ConstitutiveLaw1D* pCL, const BasicScalarFunction* pSF, const doublereal penetration,
    const StructNode* pN1, const StructNode* pN2)
: Elem(-1, flag(0)),
super(-1, pDO, pCL, pN1, pN2, Zero3, Zero3, 1.0, flag(0)),
iNumRowsNode(6),
iNumColsNode(6),
pSF(pSF),
penetration(penetration)
{
    dElle = 0.0;
}

doublereal
Collision::dCalcEpsilon(void)
{
    return dElle; //depth of collision
}

unsigned int
Collision::iGetNumDof(void) const
{
    return 0;
}

void
Collision::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = 2 * iNumRowsNode;
    *piNumCols = 2 * iNumColsNode;
}

int
Collision::ContactsSize(void)
{
    return Arm2_vector.size();
}

void
Collision::ClearContacts(void)
{
    Arm2_vector.clear();
    f1_vector.clear();
    f2_vector.clear();
    Fn_Norm_vector.clear();
}

void
Collision::PushBackContact(const btVector3& point1, const btVector3& point2)
{
    const Vec3 p1(point1[0], point1[1], point1[2]);
    const Vec3 p2(point2[0], point2[1], point2[2]);
    f1_vector.push_back(dynamic_cast<const StructNode *>(pNode1)->GetRCurr().Transpose()*(p1 - pNode1->GetXCurr()));
    f2_vector.push_back(dynamic_cast<const StructNode *>(pNode2)->GetRCurr().Transpose()*(p2 - pNode2->GetXCurr()));
    Arm2_vector.push_back(dynamic_cast<const StructNode *>(pNode2)->GetRCurr().Transpose()*(p2 + (p1 - p2) * penetration - pNode2->GetXCurr()));
    Fn_Norm_vector.push_back(0.0);
}

void
Collision::InitializeTangentialForces(void)
{
    Ft_vector.clear();
    for (int i = 0; i < 4; i++) {
        Ft_vector.push_back(Zero3);
    }
}

void
Collision::SetIndices(integer* piRow, integer* piCol)
{
    iR = *piRow;
    iC = *piCol;
    *piRow += (2 * iNumRowsNode);
    *piCol += (2 * iNumColsNode);
}

bool
Collision::SetOffsets(const int i)
{
    //printf("Entering Collision::SetOffsets()\n");
    Arm2 = Arm2_vector[i];
    f1 = f1_vector[i];
    f2 = f2_vector[i];
    const StructNode* pStructNode1(dynamic_cast<const StructNode *>(pNode1));
    const StructNode* pStructNode2(dynamic_cast<const StructNode *>(pNode2));
    Mat3x3 R1(pStructNode1->GetRCurr());
    Mat3x3 R2(pStructNode2->GetRCurr());
    Vec3 p2(pNode2->GetXCurr() + R2 * f2);
    Vec3 p1(pNode1->GetXCurr() + R1 * f1);
    v = p2 - p1;
    dElle = v.Norm();
    dEpsilon = dCalcEpsilon();
    Vec3 vPrime(pNode2->GetVCurr() + (pStructNode2->GetWCurr()).Cross(R2 * f2) - pNode1->GetVCurr() - (pStructNode1->GetWCurr()).Cross(R1 * f1));
    dEpsilonPrime  = (v.Dot(vPrime)) / v.Norm();
    return v.Dot(pNode2->GetXCurr() - pNode1->GetXCurr()) < 0.0
        && std::numeric_limits<doublereal>::epsilon() < v.Norm(); // collision?
}

VariableSubMatrixHandler&
Collision::AssJac(VariableSubMatrixHandler& WorkMat,
    doublereal dCoef,
    const VectorHandler& XCurr,
    const VectorHandler& XPrimeCurr)
{
    //printf("Entering Collision::AssJac()\n");
    DEBUGCOUT("Entering Collision::AssJac()" << std::endl);
    FullSubMatrixHandler& WM = WorkMat.SetFull();
    const integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
    const integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
    const integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
    const integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
    for (int iCnt = 1; iCnt <= iNumRowsNode; iCnt++) {
        WM.PutRowIndex(iR + iCnt, iNode1FirstMomIndex + iCnt);
        WM.PutRowIndex(iR + iNumRowsNode + iCnt, iNode2FirstMomIndex + iCnt);
    }
    for (int iCnt = 1; iCnt <= iNumColsNode; iCnt++) {
        WM.PutColIndex(iC + iCnt, iNode1FirstPosIndex + iCnt);
        WM.PutColIndex(iC + iNumColsNode + iCnt, iNode2FirstPosIndex + iCnt);
    }
    for (int i = 0; i < ContactsSize(); i++) {
        if (SetOffsets(i)) {
            AssMat(WM, dCoef, XCurr, i);
        }
    }
    return WorkMat;
}

void
Collision::AssMat(FullSubMatrixHandler& WM,
    doublereal dCoef,
    const VectorHandler& XCurr,
    const int i)
{
    const StructNode* pStructNode1(dynamic_cast<const StructNode *>(pNode1));
    const StructNode* pStructNode2(dynamic_cast<const StructNode *>(pNode2));

    /* Impact */
    const Vec3 f1Tmp((pStructNode1->GetRRef()) * f1);
    const Vec3 f2Tmp((pStructNode2->GetRRef()) * f2);
    const Vec3& Omega1(pStructNode1->GetWRef());
    const Vec3& Omega2(pStructNode2->GetWRef());
    Vec3 vPrime(pNode2->GetVCurr() + Omega2.Cross(f2Tmp) - pNode1->GetVCurr() - Omega1.Cross(f1Tmp));
    doublereal dF = GetF() / ContactsSize();
    doublereal dFDE = GetFDE();
    doublereal dFDEPrime = GetFDEPrime();

    /* Vettore forza */
    Vec3 F = v*(dF/dElle);

    Mat3x3 K(v.Tens(v * (dCoef * (dFDE / dL0 - (dEpsilonPrime * dFDEPrime + dF) / dElle)/(dElle * dElle))));
    if (dFDEPrime != 0.) {
        K += v.Tens(vPrime * (dCoef * dFDEPrime / (dElle * dElle * dL0)));
    }
    doublereal d = dCoef * dF / dElle;
    for (unsigned iCnt = 1; iCnt <= 3; iCnt++) {
        K(iCnt, iCnt) += d;
    }

    Mat3x3 KPrime;
    if (dFDEPrime != 0.) {
        KPrime = v.Tens(v * (dFDEPrime / (dL0 * dElle * dElle)));
    }

    /* Termini di forza diagonali */
    Mat3x3 Tmp1(K);
    if (dFDEPrime != 0.) {
        Tmp1 += KPrime;
    }
    WM.Add(iR + 1, iC + 1, Tmp1);
    WM.Add(iR + 6 + 1, iC + 6 + 1, Tmp1);

    /* Termini di coppia, nodo 1 */
    Mat3x3 Tmp2 = f1Tmp.Cross(Tmp1);
    WM.Add(iR + 3 + 1, iC + 1, Tmp2);
    WM.Sub(iR + 3 + 1, iC + 6 + 1, Tmp2);

    /* Termini di coppia, nodo 2 */
    Tmp2 = f2Tmp.Cross(Tmp1);
    WM.Add(iR + 9 + 1, iC + 6 + 1, Tmp2);
    WM.Sub(iR + 9 + 1, iC + 1, Tmp2);

    /* termini di forza extradiagonali */
    WM.Sub(iR + 1, iC + 6 + 1, Tmp1);
    WM.Sub(iR + 6 + 1, iC + 1, Tmp1);

    /* Termini di rotazione, Delta g1 */
    Mat3x3 Tmp3 = Tmp1 * Mat3x3(MatCross, -f1Tmp);
    if (dFDEPrime != 0.) {
        Tmp3 += KPrime * Mat3x3(MatCross, f1Tmp.Cross(Omega1 * dCoef));
    }
    WM.Add(iR + 1, iC + 3 + 1, Tmp3);
    WM.Sub(iR + 6 + 1, iC + 3 + 1, Tmp3);

    /* Termini di coppia, Delta g1 */
    Tmp2 = f1Tmp.Cross(Tmp3) + Mat3x3(MatCrossCross, F, f1Tmp * dCoef);
    WM.Add(iR + 3 + 1, iC + 3 + 1, Tmp2);
    Tmp2 = f2Tmp.Cross(Tmp3);
    WM.Sub(iR + 9 + 1, iC + 3 + 1, Tmp2);

    /* Termini di rotazione, Delta g2 */
    Tmp3 = Tmp1*Mat3x3(MatCross, -f2Tmp);
    if (dFDEPrime != 0.) {
        Tmp3 += KPrime * Mat3x3(MatCross, f2Tmp.Cross(Omega2 * dCoef));
    }
    WM.Add(iR + 6 + 1, iC + 9 + 1, Tmp3);
    WM.Sub(iR + 1, iC + 9 + 1, Tmp3);

    /* Termini di coppia, Delta g2 */
    Tmp2 = f2Tmp.Cross(Tmp3) + Mat3x3(MatCrossCross, F, f2Tmp * dCoef);
    WM.Add(iR + 9 + 1, iC + 9 + 1, Tmp2);
    Tmp2 = f1Tmp.Cross(Tmp3);
    WM.Sub(iR + 3 + 1, iC + 9 + 1, Tmp2);

    /* Resistance */
    Vec3 R_Arm2(pStructNode2->GetRCurr() * Arm2);
    Vec3 R_Arm1_calc(pNode2->GetXCurr() + R_Arm2 - pNode1->GetXCurr());
    WM.Sub(iR + 4, iC + 4, Mat3x3(MatCrossCross, Ft_vector[i] * dCoef, R_Arm1_calc));
    WM.Add(iR + 10, iC + 10, Mat3x3(MatCrossCross, Ft_vector[i] * dCoef, R_Arm2));
}

SubVectorHandler& 
Collision::AssRes(SubVectorHandler& WorkVec,
    doublereal dCoef,
    const VectorHandler& XCurr, 
    const VectorHandler& XPrimeCurr)
{
    //printf("Entering Collision::AssRes()\n");
    DEBUGCOUT("Entering Collision::AssRes()" << std::endl);
    integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
    integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

    for (int iCnt = 1; iCnt <= iNumRowsNode; iCnt++) {
      WorkVec.PutRowIndex(iR + iCnt, iNode1FirstMomIndex + iCnt);
      WorkVec.PutRowIndex(iR + iNumRowsNode + iCnt, iNode2FirstMomIndex + iCnt);
    }
    bool ChangeJac(false);
    for (int i = 0; i < ContactsSize(); i++) {
        if (SetOffsets(i)) {
            try {
                AssVec(WorkVec, dCoef, XCurr, i);
            } catch (Elem::ChangedEquationStructure) {
                ChangeJac = true;
            }
        }
    }
    if (ChangeJac) {
        throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
    }
    return WorkVec;
}

void
Collision::AssVec(SubVectorHandler& WorkVec,
    doublereal dCoef,
    const VectorHandler& XCurr,
    const int i)
{
    DEBUGCOUT("RodWithOffset::AssVec()" << std::endl);
    const StructNode* pStructNode1(dynamic_cast<const StructNode *>(pNode1));
    const StructNode* pStructNode2(dynamic_cast<const StructNode *>(pNode2));
    Mat3x3 R1(pStructNode1->GetRCurr());
    Mat3x3 R2(pStructNode2->GetRCurr());
    bool ChangeJac(false);
    
    /* Impact */
    try {
        ConstitutiveLaw1DOwner::Update(dEpsilon, dEpsilonPrime);

    } catch (Elem::ChangedEquationStructure) {
        ChangeJac = true;
    }
    Fn_Norm_vector[i] = GetF() / ContactsSize();
    Vec3 Fn(v * (Fn_Norm_vector[i] / v.Norm()));
    WorkVec.Add(iR + 1, Fn);
    WorkVec.Add(iR + 3 + 1, (R1 * f1).Cross(Fn));
    WorkVec.Sub(iR + 6 + 1, Fn);
    WorkVec.Sub(iR + 9 + 1, (R2 * f2).Cross(Fn));

    /* Resistance */
    if (pSF != NULL) {
        Vec3 unit_normal(pNode2->GetXCurr() + R2 * f2 - pNode1->GetXCurr() - R1 * f1);
        unit_normal /= unit_normal.Norm();
        Vec3 R_Arm2(R2*Arm2);
        Vec3 R_Arm1_calc(pNode2->GetXCurr() + R_Arm2 - pNode1->GetXCurr());
        Vec3 Vt(pNode2->GetVCurr() - R_Arm2.Cross(pStructNode2->GetWCurr()) - pNode1->GetVCurr()
            + R_Arm1_calc.Cross(pStructNode1->GetWCurr()));
        Vt = Vt - unit_normal*Vt.Dot(unit_normal);
        const DynamicStructNode* pDynamicStructNode1(dynamic_cast<const DynamicStructNode *>(pNode1));
        const DynamicStructNode* pDynamicStructNode2(dynamic_cast<const DynamicStructNode *>(pNode2));
        const DynamicStructDispNode* pDynamicStructDispNode1(dynamic_cast<const DynamicStructDispNode *>(pNode1));
        const DynamicStructDispNode* pDynamicStructDispNode2(dynamic_cast<const DynamicStructDispNode *>(pNode2));
        Vec3 G1(Zero3);
        Vec3 B1(Zero3);
        Vec3 G2(Zero3);
        Vec3 B2(Zero3);
        if (pDynamicStructNode1 != 0) {
            G1 = pDynamicStructNode1->GetGCurr();
        }
        if (pDynamicStructDispNode1 != 0) {
            B1 = pDynamicStructDispNode1->GetBCurr();
        }
        if (pDynamicStructNode2 != 0) {
            G2 = pDynamicStructNode2->GetGCurr();
        }
        if (pDynamicStructDispNode2 != 0) {
            B2 = pDynamicStructDispNode2->GetBCurr();
        }
        Vec3 Ft1_dCoef((R_Arm1_calc.Cross(B1) + G1).Cross(R_Arm1_calc));
        Vec3 Ft2_dCoef((R_Arm2.Cross(B2) + G2).Cross(R_Arm2));
        Vec3 Ft((Ft2_dCoef - Ft1_dCoef) / dCoef);
        Ft = Ft - unit_normal*Ft.Dot(unit_normal);
        doublereal Ft_Norm = Ft.Norm();
        doublereal Ft_max_Norm = (*pSF)(Vt.Norm()) * Fn_Norm_vector[i];
        Vec3 Ft_max = Ft * (Ft_max_Norm / Ft_Norm);
        if (Ft_max_Norm < Ft_Norm) {
            Ft_vector[i] = Ft_max;
        } else {
            Ft_vector[i] = (Ft_vector[i] + Ft) * 0.5;
            if (Ft_max_Norm < Ft_vector[i].Norm()) {
                Ft_vector[i] = Ft_max;
            }
        }
        WorkVec.Add(iR + 1, Ft_vector[i]);
        WorkVec.Add(iR + 4, R_Arm1_calc.Cross(Ft_vector[i]));
        WorkVec.Sub(iR + 7, Ft_vector[i]);
        WorkVec.Sub(iR + 10, R_Arm2.Cross(Ft_vector[i]));
    }
    if (ChangeJac) {
        throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
    }
}

std::ostream&
Collision::OutputAppend(std::ostream& out, int i) const {
    out << " " << pNode1->GetLabel();
    out << " " << pNode2->GetLabel();
    for (int iCnt = 1; iCnt <= 3; iCnt++) {
        out << " " << f1_vector[i](iCnt);
    }
    for (int iCnt = 1; iCnt <= 3; iCnt++) {
        out << " " << f2_vector[i](iCnt);
    }
    for (int iCnt = 1; iCnt <= 3; iCnt++) {
        out << " " << Ft_vector[i](iCnt);
    }
    out << " " << Fn_Norm_vector[i];
    //ConstitutiveLaw1DOwner::OutputAppend(out);
}

class Universe {
public:
    std::map<unsigned, btCollisionObject*> label_object_map;
    std::map<btCollisionObject*, const StructNode*> object_node_map;
    std::map<btCollisionObject*, std::string> object_material_map;
} universe;

// CollisionWorld: begin

CollisionWorld::CollisionWorld(
    unsigned uLabel, const DofOwner *pDO,
    DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
    printf("Entering CollisionWorld::CollisionWorld()\n");
    if (HP.IsKeyWord("help")) {
        silent_cout(
            "\n"
            "Module:     Collision\n"
            "\n"
            "    This element implements a collision world\n"
            "\n"
            "    collision world,\n"
            "       (integer)<number_of_material_pairs>, (integer)<number_of_collision_objects\n"
            "           <material_pair> [,...]\n"
            "           (CollisionObject) <label> [,...]\n"
            "\n"
            "    <material_pair> ::= (str)<material1>, (str)<material2>, (ConstitutiveLaw<1D>)<const_law> [, friction function, (ScalarFunction)<SF>, (real)<penetration>]\n"
            "\n\n"
            << std::endl);

        if (!HP.IsArg()) {
            /*
             * Exit quietly if nothing else is provided
             */
            throw NoErr(MBDYN_EXCEPT_ARGS);
        }
    }
    configuration = new btDefaultCollisionConfiguration();
    dispatcher = new btCollisionDispatcher(configuration);
    broadphase = new btDbvtBroadphase();
    world = new btCollisionWorld(dispatcher, broadphase, configuration);
    int N_CL = HP.GetInt();
    int N_CO = HP.GetInt();
    ConstLawType::Type VECLType(ConstLawType::VISCOELASTIC);
    typedef std::pair<std::string, std::string> MaterialPair;
    std::map<MaterialPair, ConstitutiveLaw1D*> pCL;
    std::map<MaterialPair, const BasicScalarFunction*> pSF;
    std::map<MaterialPair, doublereal> penetration;
    for (int i; i < N_CL; i++) {
        MaterialPair material_pair(std::make_pair(HP.GetValue(TypedValue::VAR_STRING).GetString(), HP.GetValue(TypedValue::VAR_STRING).GetString()));
        pCL[material_pair] = HP.GetConstLaw1D(VECLType);
        if (HP.IsKeyWord("friction" "function")) {
            pSF[material_pair] = ParseScalarFunction(HP, pDM);
            penetration[material_pair] = HP.GetReal();
        } else {
            pSF[material_pair] = NULL;
            penetration[material_pair] = 0.0;
        }
    }
    std::set<btCollisionObject*> all_objects;
    std::set<btCollisionObject*> objects;
    for (int i; i < N_CO; i++) {
        btCollisionObject* ob = universe.label_object_map[HP.GetInt()];
        for (std::set<btCollisionObject*>::iterator it = all_objects.begin();
            it != all_objects.end(); it++) {
            MaterialPair material_pair(std::make_pair(
                universe.object_material_map[ob], universe.object_material_map[*it]));
            if (pCL.find(material_pair) == pCL.end()) {
                std::swap(material_pair.first, material_pair.second);
            }
            if (pCL.find(material_pair) != pCL.end()) {
                objectpair_collision_map[std::make_pair(ob, *it)] = new Collision(
                    pDO, pCL[material_pair], pSF[material_pair], penetration[material_pair],
                    universe.object_node_map[ob], universe.object_node_map[*it]);
                objects.insert(ob);
                objects.insert(*it);
            }
        }
        all_objects.insert(ob);
    }
    for (std::set<btCollisionObject*>::iterator it = objects.begin();
        it != objects.end(); it++) {
        nodes.insert(universe.object_node_map[*it]);
        world->addCollisionObject((*it));
    }
    iNumRows = 0;
    iNumCols = 0;
    for (std::map<ObjectPair, Collision*>::const_iterator it = objectpair_collision_map.begin();
        it != objectpair_collision_map.end(); it++) {
        it->second->SetIndices(&iNumRows, &iNumCols);
    }
    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

CollisionWorld::~CollisionWorld(void)
{
    delete world;
    delete broadphase;
    delete dispatcher;
    delete configuration;
}

unsigned int
CollisionWorld::iGetNumDof(void) const
{
    return 0;
}

DofOrder::Order
CollisionWorld::GetDofType(unsigned int i) const {
  ASSERT(i >= 0 && i < iGetNumDof());
  return DofOrder::ALGEBRAIC;
};

DofOrder::Order
CollisionWorld::GetEqType(unsigned int i) const
{
    ASSERTMSGBREAK(i >=0 and i < iGetNumDof(), 
        "INDEX ERROR in SphericalHingeJoint::GetEqType");
    return DofOrder::ALGEBRAIC;
}

void
CollisionWorld::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    //printf("Entering CollisionWorld::WorkSpaceDim()\n");
    *piNumRows = iNumRows;
    *piNumCols = iNumCols;
}

int
CollisionWorld::iGetNumConnectedNodes(void) const
{
    printf("Entering CollisionWorld::iGetNumConnectedNodes()\n");
    return nodes.size();
}

void
CollisionWorld::GetConnectedNodes(std::vector<const Node*>& connectedNodes) const
{
    printf("Entering CollisionWorld::GetConnectedNodes()\n");
    connectedNodes.resize(iGetNumConnectedNodes());
    integer i(0);
    for (std::set<const Node*>::const_iterator it = nodes.begin(); it != nodes.end(); it++) {
        connectedNodes[i++] = *it;
    }
}

void
CollisionWorld::Output(OutputHandler& OH) const
{
    //printf("Entering CollisionWorld::Output()\n");
    if (fToBeOutput()) {
        if ( OH.UseText(OutputHandler::LOADABLE) ) {
            std::ostream& os = OH.Loadable();
            os << GetLabel();
            os << ss.str();
            os << std::endl;
        }
    }
}

void
CollisionWorld::InitializeTangentialForces(void)
{
    //printf("Entering CollisionWorld::InitializeTangentialForces()\n");
    for (std::map<ObjectPair, Collision*>::const_iterator it = objectpair_collision_map.begin();
        it != objectpair_collision_map.end(); it++) {
        it->second->InitializeTangentialForces();
    }
}

void
CollisionWorld::ClearAndPushContacts(void)
{
    //printf("Entering CollisionWorld::ClearAndPushContacts()\n");
    for (std::map<ObjectPair, Collision*>::const_iterator it = objectpair_collision_map.begin();
        it != objectpair_collision_map.end(); it++) {
        it->second->ClearContacts();
    }
    world->performDiscreteCollisionDetection();
    for (int i = 0; i < dispatcher->getNumManifolds() ; i++) {
        btPersistentManifold* manifold(dispatcher->getManifoldByIndexInternal(i));
        ObjectPair object_pair = std::make_pair(
            (btCollisionObject*)(manifold->getBody0()), (btCollisionObject*)(manifold->getBody1()));
        bool swapped(false);
        if (objectpair_collision_map.count(object_pair) == 0) {
            std::swap(object_pair.first, object_pair.second);
            swapped = true;
        }
        if (objectpair_collision_map.count(object_pair)) {
            for (int k = 0; k < manifold->getNumContacts(); k++) {
                btManifoldPoint& pt(manifold->getContactPoint(k));
                btVector3 p1;
                btVector3 p2;
                if (swapped) {
                    p1 = pt.getPositionWorldOnB();
                    p2 = pt.getPositionWorldOnA();
                } else {
                    p1 = pt.getPositionWorldOnA();
                    p2 = pt.getPositionWorldOnB();
                }
                objectpair_collision_map[object_pair]->PushBackContact(p1, p2);
            }
        }
    }
}

void
CollisionWorld::SetValue(DataManager *pDM,
    VectorHandler& X, VectorHandler& XP,
    SimulationEntity::Hints *ph)
{
    //printf("Entering CollisionWorld::SetValue()\n");
    InitializeTangentialForces();
}

void
CollisionWorld::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
    //printf("Entering CollisionWorld::AfterPredict()\n");
}

void
CollisionWorld::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
    //printf("Entering CollisionWorld::AfterConvergence()\n");
    ss.str("");
    ss.clear();
    for (std::map<ObjectPair, Collision*>::const_iterator it = objectpair_collision_map.begin();
        it != objectpair_collision_map.end(); it++) {
        for (int i = 0; i < it->second->ContactsSize(); i++) {
            it->second->OutputAppend(ss, i);
        }
    }
    InitializeTangentialForces();
}

SubVectorHandler& 
CollisionWorld::AssRes(SubVectorHandler& WorkVec,
    doublereal dCoef,
    const VectorHandler& XCurr, 
    const VectorHandler& XPrimeCurr)
{
    //printf("Entering CollisionWorld::AssRes()\n");
    DEBUGCOUT("Entering CollisionWorld::AssRes()" << std::endl);
    ClearAndPushContacts();
    WorkVec.ResizeReset(iNumRows);
    integer iRow(0);
    integer iCol(0);
    bool ChangeJac(false);
    for (std::map<ObjectPair, Collision*>::const_iterator it = objectpair_collision_map.begin();
        it != objectpair_collision_map.end(); it++) {
        try {
            it->second->SetIndices(&iRow, &iCol);
            it->second->AssRes(WorkVec, dCoef, XCurr, XPrimeCurr);
        } catch (Elem::ChangedEquationStructure) {
            ChangeJac = true;
            //printf("Changed Jac\n");
        }
    }
    if (ChangeJac) {
        throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
    }
    return WorkVec;
}

VariableSubMatrixHandler& 
CollisionWorld::AssJac(VariableSubMatrixHandler& WorkMat,
    doublereal dCoef, 
    const VectorHandler& XCurr,
    const VectorHandler& XPrimeCurr)
{
    //printf("Entering CollisionWorld::AssJac()\n");
    DEBUGCOUT("Entering CollisionWorld::AssJac()" << std::endl);
    FullSubMatrixHandler& WM = WorkMat.SetFull();
    WM.ResizeReset(iNumRows, iNumCols);
    integer iRow(0);
    integer iCol(0);
    for (std::map<ObjectPair, Collision*>::const_iterator it = objectpair_collision_map.begin();
        it != objectpair_collision_map.end(); it++) {
        it->second->SetIndices(&iRow, &iCol);
        it->second->AssJac(WorkMat, dCoef, XCurr, XPrimeCurr);
    }
    return WorkMat;
}

unsigned int
CollisionWorld::iGetNumPrivData(void) const
{
    printf("Entering CollisionWorld::iGetNumPrivData()\n");
    // return number of private data
    return 0;
}

unsigned int
CollisionWorld::iGetPrivDataIdx(const char *s) const
{
    printf("Entering CollisionWorld::iGetPrivDataIdx()\n");
    // parse string and compute index of requested private data

    // shouldn't get here until private data are defined
    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

doublereal
CollisionWorld::dGetPrivData(unsigned int i) const
{
    printf("Entering CollisionWorld::dGetPrivData()\n");
    ASSERT(0 < i && i <= iGetNumPrivData());

    switch (i) {
        case 0:
        default:
            silent_cerr("collision world(" << GetLabel() << "): invalid private data index " << i << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    // shouldn't get here until private data are defined
    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

std::ostream&
CollisionWorld::Restart(std::ostream& out) const
{
    return out << "# CollisionWorld: not implemented" << std::endl;
}

unsigned int
CollisionWorld::iGetInitialNumDof(void) const
{
    printf("Entering CollisionWorld::iGetInitialNumDof()\n");
    return 0;
}

void 
CollisionWorld::InitialWorkSpaceDim(
    integer* piNumRows,
    integer* piNumCols) const
{
    printf("Entering CollisionWorld::InitialWorkSpaceDim()\n");
    *piNumRows = 0;
    *piNumCols = 0;
}

VariableSubMatrixHandler&
CollisionWorld::InitialAssJac(
    VariableSubMatrixHandler& WorkMat, 
    const VectorHandler& XCurr)
{
    // should not be called, since initial workspace is empty
    ASSERT(0);
    //printf("Entering CollisionWorld::InitialAssJac()\n");
    DEBUGCOUT("Entering CollisionWorld::InitialAssJac()" << std::endl);

    WorkMat.SetNullMatrix();

    return WorkMat;
}

SubVectorHandler& 
CollisionWorld::InitialAssRes(
    SubVectorHandler& WorkVec,
    const VectorHandler& XCurr)
{
    // should not be called, since initial workspace is empty
    ASSERT(0);
    //printf("Entering CollisionWorld::InitialAssRes()\n");
    DEBUGCOUT("Entering CollisionWorld::InitialAssRes()" << std::endl);

    WorkVec.ResizeReset(0);

    return WorkVec;
}

// CollisionWorld: end

// CollisionObject: begin

CollisionObject::CollisionObject(
    unsigned uLabel, const DofOwner *pDO,
    DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
    if (HP.IsKeyWord("help")) {
        silent_cout(
            "\n"
            "Module:     Collision\n"
            "\n"
            "    This element implements a collision object\n"
            "\n"
            "    collision object,\n"
            "        (Node) <label>,\n"
            "            (Vec3) <offset>,\n"
            "            (Mat3x3) <orientation>,\n"
            "        (str)<material>,\n"
            "        <shape> [,margin, (real)<margin>]\n"
            "\n"
            "   <shape> ::= {\n"
            "       btBoxShape, (real)<x_half_extent>, (real)<y_half_extent>, (real)<z_half_extent>\n"
            "       | btCapsuleShape, (real)<radius>, (real)<height>\n"
            "       | btConeShape, (real)<radius>, (real)<height>\n"
            "       | btSphereShape, (real)<radius>\n"
            "   }\n\n"
            << std::endl);

        if (!HP.IsArg()) {
            /*
             * Exit quietly if nothing else is provided
             */
            throw NoErr(MBDYN_EXCEPT_ARGS);
        }
    }
    pNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
    const ReferenceFrame RF(pNode);
    f = HP.GetPosRel(RF);
    Rh = HP.GetRotRel(RF);
    Vec3 x(pNode->GetXCurr() + pNode->GetRCurr()*f);
    Mat3x3 r(pNode->GetRCurr()*Rh);
    ob = new btCollisionObject();
    btTransform& transform = ob->getWorldTransform();
    transform.setOrigin(btVector3(x[0], x[1], x[2]));
    transform.setBasis(btMatrix3x3(r.dGet(1,1),r.dGet(1,2),r.dGet(1,3),r.dGet(2,1),r.dGet(2,2),r.dGet(2,3),r.dGet(3,1),r.dGet(3,2),r.dGet(3,3)));
    const TypedValue material(HP.GetValue(TypedValue::VAR_STRING));
    if (HP.IsKeyWord("btBoxShape")) {
        btScalar x(HP.GetReal());
        btScalar y(HP.GetReal());
        btScalar z(HP.GetReal());
        shape = new btBoxShape(btVector3(x, y, z));
        printf("btBoxShape");
    } else if (HP.IsKeyWord("btCapsuleShape")) {
        btScalar radius(HP.GetReal());
        btScalar height(HP.GetReal());
        shape = new btCapsuleShape(radius*.96, height*.92);
        printf("btCapsuleShape");
    } else if (HP.IsKeyWord("btConeShape")) {
        btScalar radius(HP.GetReal());
        btScalar height(HP.GetReal());
        shape = new btConeShape(radius, height);
        printf("btConeShape");
    } else if (HP.IsKeyWord("btSphereShape")) {
        btScalar radius(HP.GetReal());
        shape = new btSphereShape(radius);
        printf("btSphereShape");
    } else if (HP.IsKeyWord("btStaticPlaneShape")) {
        const btVector3 z_plane(0., 0., 1.);
        shape = new btStaticPlaneShape(z_plane, 0.);
        printf("btStaticPlaneShape");
    } else {
        silent_cerr("collision object(" << GetLabel() << "): a valid btShape type is expected at line " << HP.GetLineData() << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
    if (HP.IsKeyWord("margin")) {
        btScalar margin(HP.GetReal());
        shape->setMargin(margin);
    }
    printf(", getMargin() %f\n", shape->getMargin());
    ob->setCollisionShape(shape);
    universe.label_object_map[uLabel] = ob;
    universe.object_node_map[ob] = pNode;
    universe.object_material_map[ob] = material.GetString();
    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

CollisionObject::~CollisionObject(void)
{
    delete ob;
    delete shape;
}

void
CollisionObject::Output(OutputHandler& OH) const
{
    if (fToBeOutput()) {
        NO_OP;
    }
}

void
CollisionObject::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = 0;
    *piNumCols = 0;
}

VariableSubMatrixHandler& 
CollisionObject::AssJac(VariableSubMatrixHandler& WorkMat,
    doublereal dCoef, 
    const VectorHandler& XCurr,
    const VectorHandler& XPrimeCurr)
{
    DEBUGCOUT("Entering CollisionObject::AssJac()" << std::endl);
    WorkMat.SetNullMatrix();

    return WorkMat;
}

SubVectorHandler& 
CollisionObject::AssRes(SubVectorHandler& WorkVec,
    doublereal dCoef,
    const VectorHandler& XCurr, 
    const VectorHandler& XPrimeCurr)
{
    DEBUGCOUT("Entering CollisionObject::AssRes()" << std::endl);
    WorkVec.ResizeReset(0);

    Vec3 x(pNode->GetXCurr() + pNode->GetRCurr() * f);
    Mat3x3 r(pNode->GetRCurr() * Rh);
    btTransform& transform = ob->getWorldTransform();
    transform.setOrigin(btVector3(x[0], x[1], x[2]));
    transform.setBasis(btMatrix3x3(r.dGet(1,1),r.dGet(1,2),r.dGet(1,3),r.dGet(2,1),r.dGet(2,2),r.dGet(2,3),r.dGet(3,1),r.dGet(3,2),r.dGet(3,3)));
    
    return WorkVec;
}

unsigned int
CollisionObject::iGetNumPrivData(void) const
{
    return 0;
}

unsigned int
CollisionObject::iGetPrivDataIdx(const char *s) const
{
    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

doublereal
CollisionObject::dGetPrivData(unsigned int i) const
{
    ASSERT(i > 1 && i <= iGetNumPrivData());
    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

int
CollisionObject::iGetNumConnectedNodes(void) const
{
    return 0;
}

void
CollisionObject::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
    connectedNodes.resize(0);
}

void
CollisionObject::SetValue(DataManager *pDM,
    VectorHandler& X, VectorHandler& XP,
    SimulationEntity::Hints *ph)
{
    NO_OP;
}

std::ostream&
CollisionObject::Restart(std::ostream& out) const
{
    return out << "# CollisionObject: not implemented" << std::endl;
}

unsigned int
CollisionObject::iGetInitialNumDof(void) const
{
    return 0;
}

void 
CollisionObject::InitialWorkSpaceDim(
    integer* piNumRows,
    integer* piNumCols) const
{
    *piNumRows = 0;
    *piNumCols = 0;
}

VariableSubMatrixHandler&
CollisionObject::InitialAssJac(
    VariableSubMatrixHandler& WorkMat, 
    const VectorHandler& XCurr)
{
    // should not be called, since initial workspace is empty
    ASSERT(0);
    DEBUGCOUT("Entering CollisionObject::InitialAssJac()" << std::endl);

    WorkMat.SetNullMatrix();

    return WorkMat;
}

SubVectorHandler& 
CollisionObject::InitialAssRes(
    SubVectorHandler& WorkVec,
    const VectorHandler& XCurr)
{
    // should not be called, since initial workspace is empty
    ASSERT(0);
    DEBUGCOUT("Entering CollisionObject::InitialAssRes()" << std::endl);

    WorkVec.ResizeReset(0);

    return WorkVec;
}

// CollisionObject: end

bool module_read(void)
{
    UserDefinedElemRead *rf1 = new UDERead<CollisionWorld>;
    if (!SetUDE("collision" "world", rf1))
    {
        delete rf1;
        return false;
    }
    UserDefinedElemRead *rf2 = new UDERead<CollisionObject>;
    if (!SetUDE("collision" "object", rf2))
    {
        delete rf1;
        delete rf2;
        return false;
    }
    return true;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
    if (!module_read()) {
        silent_cerr("Module: "
            "module_init(" << module_name << ") "
            "failed" << std::endl);
        return -1;
    }
    return 0;
}

