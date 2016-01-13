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
#include <set>
#include "rodj.h"
#include <limits>
#include "module-collision.h"

Contact::Contact(std::pair<fcl::Vec3f, fcl::Vec3f> pt_pair, const StructDispNode* pNode1, const StructDispNode* pNode2, doublereal penetration_ratio)
: Ft(Zero3), Fn_Norm(0.0), tangent(Zero3)
{
    Vec3 pt1(pt_pair.first[0], pt_pair.first[1], pt_pair.first[2]);
    Vec3 pt2(pt_pair.second[0], pt_pair.second[1], pt_pair.second[2]);
    Mat3x3 R1(dynamic_cast<const StructNode *>(pNode1)->GetRCurr());
    Mat3x3 R2(dynamic_cast<const StructNode *>(pNode2)->GetRCurr());
    f1 = R1.Transpose() * (pt1 - pNode1->GetXCurr());
    f2 = R2.Transpose() * (pt2 - pNode2->GetXCurr());
    Arm1 = R1.Transpose() * (pt1 * penetration_ratio + pt2 * (1.0 - penetration_ratio) - pNode1->GetXCurr());
}

Contact::~Contact(void)
{
    NO_OP;
}

CollisionObjectData::CollisionObjectData(const StructNode* pNode,
    fcl::CollisionObject* pObject, std::string material)
: pNode(pNode),
pObject(pObject),
material(material)
{
    NO_OP;
}

CollisionObjectData::~CollisionObjectData(void)
{
    NO_OP;
}

std::map<const unsigned, CollisionObjectData*> collision_object_data;

Collision::Collision(FCL::Func func,
    const ConstitutiveLaw1D* pCL, const BasicScalarFunction* pSF, doublereal penetration_ratio,
    const CollisionObjectData* pD1, const CollisionObjectData* pD2, integer* piRow, integer* piCol)
: func(func),
ConstitutiveLaw1DOwner(pCL),
pSF(pSF),
penetration_ratio(penetration_ratio),
pNode1(pD1->pNode),
pNode2(pD2->pNode),
pObject1(pD1->pObject),
pObject2(pD2->pObject),
iNumRowsNode(6),
iNumColsNode(6)
{
    iR = *piRow;
    iC = *piCol;
    *piRow += (2 * iNumRowsNode);
    *piCol += (2 * iNumColsNode);
}

void
Collision::ClearContacts(void)
{
    contacts.clear();
}

void
Collision::Intersect(void)
{
    FCL::Vec3f_pairs pt_pairs;
    func(pObject1, pObject2, pt_pairs);
    for (std::vector<std::pair<fcl::Vec3f, fcl::Vec3f> >::iterator it = pt_pairs.begin();
        it != pt_pairs.end(); it++) {
        contacts.push_back(Contact(*it, pNode1, pNode2, penetration_ratio));
        /* For testing 
        Vec3 pt1(it->first[0], it->first[1], it->first[2]);
        Vec3 pt2(it->second[0], it->second[1], it->second[2]);
        Mat3x3 R1(dynamic_cast<const StructNode *>(pNode1)->GetRCurr());
        Mat3x3 R2(dynamic_cast<const StructNode *>(pNode2)->GetRCurr());
        Vec3 f1 = R1.Transpose() * (pt1 - pNode1->GetXCurr());
        Vec3 f2 = R2.Transpose() * (pt2 - pNode2->GetXCurr());
        printf("(%f %f %f), (%f %f %f)\n", it->first[0], it->first[1], it->first[2], it->second[0], it->second[1], it->second[2]);
        printf("(%f %f %f), (%f %f %f)\n", f1(1), f1(2), f1(3), f2(1), f2(2), f2(3));
        */
    }
    if (tangents.size() == contacts.size()) {
        for (int i = 0; i < tangents.size(); i++) {
            contacts[i].tangent = tangents[i];
        }
    }
}

void
Collision::ClearAndSetTangents(void)
{
    tangents.clear();
    const StructNode* pStructNode1(dynamic_cast<const StructNode *>(pNode1));
    const StructNode* pStructNode2(dynamic_cast<const StructNode *>(pNode2));
    const Mat3x3 R1(pStructNode1->GetRCurr());
    const Mat3x3 R2(pStructNode2->GetRCurr());
    for (std::vector<Contact>::iterator it = contacts.begin(); it != contacts.end(); it++) {
        const Vec3 Rf1(R1 * it->f1);
        const Vec3 Rf2(R2 * it->f2);
        Vec3 normal = pNode2->GetXCurr() + Rf2 - pNode1->GetXCurr() - Rf1;
        const doublereal depth = normal.Norm();
        if (std::numeric_limits<doublereal>::epsilon() < depth) {
            normal /= depth;
        } else {
            normal = Zero3;
        }
        const Vec3 R_Arm1(R1 * it->Arm1);
        const Vec3 R_Arm2(pNode1->GetXCurr() + R_Arm1 - pNode2->GetXCurr());
        Vec3 Vt(pNode2->GetVCurr() + (pStructNode2->GetWCurr()).Cross(R_Arm2) - pNode1->GetVCurr() - (pStructNode1->GetWCurr()).Cross(R_Arm1));
        Vt -= normal * Vt.Dot(normal);
        if (std::numeric_limits<doublereal>::epsilon() < Vt.Norm()) {
            tangents.push_back(Vt / Vt.Norm());
        } else {
            tangents.push_back(Zero3);
        }
    }
}

VariableSubMatrixHandler&
Collision::AssJac(VariableSubMatrixHandler& WorkMat,
    doublereal dCoef,
    const VectorHandler& XCurr,
    const VectorHandler& XPrimeCurr)
{
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
    for (std::vector<Contact>::iterator it = contacts.begin(); it != contacts.end(); it++) {
        AssMat(WM, dCoef, *it);
    }
    return WorkMat;
}

void
Collision::AssMat(FullSubMatrixHandler& WM, doublereal dCoef, Contact& contact)
{
    const StructNode* pStructNode1(dynamic_cast<const StructNode *>(pNode1));
    const StructNode* pStructNode2(dynamic_cast<const StructNode *>(pNode2));
    const Mat3x3 R1(pStructNode1->GetRCurr());
    const Mat3x3 R2(pStructNode2->GetRCurr());
    
    /* Impact */
    const Vec3 Rf1(R1 * contact.f1);
    const Vec3 Rf2(R2 * contact.f2);
    const Vec3 V(pNode2->GetVCurr() + (pStructNode2->GetWCurr()).Cross(Rf2) - pNode1->GetVCurr() - (pStructNode1->GetWCurr()).Cross(Rf1));
    Vec3 normal = pNode2->GetXCurr() + Rf2 - pNode1->GetXCurr() - Rf1;
    const doublereal depth = normal.Norm();
    const doublereal V_Norm = V.Norm();
    if (std::numeric_limits<doublereal>::epsilon() < depth) {
        normal /= depth;
    } else if (std::numeric_limits<doublereal>::epsilon() < V_Norm) {
        normal = V / V_Norm;
    } else {
        return;
    }
    const doublereal Vn_Norm = V.Dot(normal);
    ConstitutiveLaw1DOwner::Update(depth, Vn_Norm);
    doublereal Fn_Norm = GetF() / contacts.size();
    doublereal FDEPrime = GetFDEPrime();

    /* Vettore forza */
    const Vec3 Fn = normal * Fn_Norm;

    Mat3x3 K(normal.Tens() * (dCoef * (GetFDE() + (normal.Dot(V) * FDEPrime - Fn_Norm) / depth)));
    if (FDEPrime != 0.) {
        K += normal.Tens(V) * (dCoef * FDEPrime / depth);
    }
    doublereal d = dCoef * Fn_Norm / depth;
    for (unsigned iCnt = 1; iCnt <= 3; iCnt++) {
        K(iCnt, iCnt) += d;
    }

    Mat3x3 KPrime;
    if (FDEPrime != 0.) {
        KPrime = normal.Tens() * FDEPrime;
    }

    /* Termini di forza diagonali */
    Mat3x3 Tmp1(K);
    if (FDEPrime != 0.) {
        Tmp1 += KPrime;
    }
    WM.Add(iR + 1, iC + 1, Tmp1);
    WM.Add(iR + 7, iC + 7, Tmp1);

    /* Termini di coppia, nodo 1 */
    Mat3x3 Tmp2 = Rf1.Cross(Tmp1);
    WM.Add(iR + 4, iC + 1, Tmp2);
    WM.Sub(iR + 4, iC + 7, Tmp2);

    /* Termini di coppia, nodo 2 */
    Tmp2 = Rf2.Cross(Tmp1);
    WM.Add(iR + 10, iC + 7, Tmp2);
    WM.Sub(iR + 10, iC + 1, Tmp2);

    /* termini di forza extradiagonali */
    WM.Sub(iR + 1, iC + 7, Tmp1);
    WM.Sub(iR + 7, iC + 1, Tmp1);

    /* Termini di rotazione, Delta g1 */
    Mat3x3 Tmp3 = Tmp1 * Mat3x3(MatCross, -Rf1);
    if (FDEPrime != 0.) {
        Tmp3 += KPrime * Mat3x3(MatCross, Rf1.Cross((pStructNode1->GetWRef()) * dCoef));
    }
    WM.Add(iR + 1, iC + 4, Tmp3);
    WM.Sub(iR + 7, iC + 4, Tmp3);

    /* Termini di coppia, Delta g1 */
    Tmp2 = Rf1.Cross(Tmp3) + Mat3x3(MatCrossCross, Fn, Rf1 * dCoef);
    WM.Add(iR + 4, iC + 4, Tmp2);
    Tmp2 = Rf2.Cross(Tmp3);
    WM.Sub(iR + 10, iC + 4, Tmp2);

    /* Termini di rotazione, Delta g2 */
    Tmp3 = Tmp1*Mat3x3(MatCross, -Rf2);
    if (FDEPrime != 0.) {
        Tmp3 += KPrime * Mat3x3(MatCross, Rf2.Cross((pStructNode2->GetWRef()) * dCoef));
    }
    WM.Add(iR + 7, iC + 10, Tmp3);
    WM.Sub(iR + 1, iC + 10, Tmp3);

    /* Termini di coppia, Delta g2 */
    Tmp2 = Rf2.Cross(Tmp3) + Mat3x3(MatCrossCross, Fn, Rf2 * dCoef);
    WM.Add(iR + 10, iC + 10, Tmp2);
    Tmp2 = Rf1.Cross(Tmp3);
    WM.Sub(iR + 4, iC + 10, Tmp2);

    /* Resistance */
    if (pSF != NULL) {
        const Vec3 R_Arm1(R1 * contact.Arm1);
        const Vec3 R_Arm2(pNode1->GetXCurr() + R_Arm1 - pNode2->GetXCurr());
        const doublereal Ft_Norm_max = (*pSF)((V - normal * Vn_Norm).Norm()) * contact.Fn_Norm;
        Vec3 Ft = contact.tangent * Ft_Norm_max;
        WM.Sub(iR + 4, iC + 4, Mat3x3(MatCrossCross, Ft * dCoef, R_Arm1));
        WM.Add(iR + 10, iC + 10, Mat3x3(MatCrossCross, Ft * dCoef, R_Arm2));
    }

}

SubVectorHandler& 
Collision::AssRes(SubVectorHandler& WorkVec,
    doublereal dCoef,
    const VectorHandler& XCurr, 
    const VectorHandler& XPrimeCurr)
{
    DEBUGCOUT("Entering Collision::AssRes()" << std::endl);
    integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
    integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

    for (int iCnt = 1; iCnt <= iNumRowsNode; iCnt++) {
      WorkVec.PutRowIndex(iR + iCnt, iNode1FirstMomIndex + iCnt);
      WorkVec.PutRowIndex(iR + iNumRowsNode + iCnt, iNode2FirstMomIndex + iCnt);
    }
    for (std::vector<Contact>::iterator it = contacts.begin(); it != contacts.end(); it++) {
        AssVec(WorkVec, dCoef, *it);
    }
    return WorkVec;
}

void
Collision::AssVec(SubVectorHandler& WorkVec, doublereal dCoef, Contact& contact)
{
    DEBUGCOUT("RodWithOffset::AssVec()" << std::endl);
    const StructNode* pStructNode1(dynamic_cast<const StructNode *>(pNode1));
    const StructNode* pStructNode2(dynamic_cast<const StructNode *>(pNode2));
    const Mat3x3 R1(pStructNode1->GetRCurr());
    const Mat3x3 R2(pStructNode2->GetRCurr());
    
    /* Impact */
    const Vec3 Rf1(R1 * contact.f1);
    const Vec3 Rf2(R2 * contact.f2);
    const Vec3 V(pNode2->GetVCurr() + (pStructNode2->GetWCurr()).Cross(Rf2) - pNode1->GetVCurr() - (pStructNode1->GetWCurr()).Cross(Rf1));
    Vec3 normal = pNode2->GetXCurr() + Rf2 - pNode1->GetXCurr() - Rf1;
    const doublereal depth = normal.Norm();
    const doublereal V_Norm = V.Norm();
    if (std::numeric_limits<doublereal>::epsilon() < depth) {
        normal /= depth;
    } else if (std::numeric_limits<doublereal>::epsilon() < V_Norm) {
        normal = V / V_Norm;
    } else {
        return;
    }
    const doublereal Vn_Norm = V.Dot(normal);
    ConstitutiveLaw1DOwner::Update(depth, Vn_Norm);
    contact.Fn_Norm = GetF() / contacts.size();
    const Vec3 Fn(normal * contact.Fn_Norm);
    WorkVec.Add(iR + 1, Fn);
    WorkVec.Add(iR + 4, Rf1.Cross(Fn));
    WorkVec.Sub(iR + 7, Fn);
    WorkVec.Sub(iR + 10, Rf2.Cross(Fn));

    /* Resistance */
    if (pSF != NULL) {
        const Vec3 R_Arm1(R1 * contact.Arm1);
        const Vec3 R_Arm2(pNode1->GetXCurr() + R_Arm1 - pNode2->GetXCurr());
        const doublereal Ft_Norm_max = (*pSF)((V - normal * Vn_Norm).Norm()) * contact.Fn_Norm;
        contact.Ft = contact.tangent * Ft_Norm_max;
        WorkVec.Add(iR + 1, contact.Ft);
        WorkVec.Add(iR + 4, R_Arm1.Cross(contact.Ft));
        WorkVec.Sub(iR + 7, contact.Ft);
        WorkVec.Sub(iR + 10, R_Arm2.Cross(contact.Ft));
    } else {
        contact.Ft = Zero3;
    }
}

std::ostream&
Collision::OutputAppend(std::ostream& out) const {
    for (std::vector<Contact>::const_iterator it = contacts.begin(); it != contacts.end(); it++) {
        out << " " << pNode1->GetLabel();
        out << " " << pNode2->GetLabel();
        for (int iCnt = 1; iCnt <= 3; iCnt++) {
            out << " " << it->f1(iCnt);
        }
        for (int iCnt = 1; iCnt <= 3; iCnt++) {
            out << " " << it->f2(iCnt);
        }
        for (int iCnt = 1; iCnt <= 3; iCnt++) {
            out << " " << it->Ft(iCnt);
        }
        out << " " << it->Fn_Norm;
        // If desired, must recalculate depth and Vn_Norm
        //ConstitutiveLaw1DOwner::Update(depth, Vn_Norm);
        //ConstitutiveLaw1DOwner::OutputAppend(out);
    }
}

bool CollisionFunction(fcl::CollisionObject* o1, fcl::CollisionObject* o2, void* cdata_)
{
    std::map<const FCL::ObjectPair, Collision*> objectpair_collision_map(
        *(static_cast<std::map<const FCL::ObjectPair, Collision*>*>(cdata_)));
    FCL::ObjectPair object_pair(o1, o2);
    if (objectpair_collision_map.count(object_pair) == 0) {
        std::swap(object_pair.first, object_pair.second);
        if (objectpair_collision_map.count(object_pair) == 0) {
            return false;
        }
    }
    objectpair_collision_map[object_pair]->Intersect();
    return false;
}


// CollisionWorld: begin

CollisionWorld::CollisionWorld(
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
            "    This element implements a collision world\n"
            "\n"
            "    collision world,\n"
            "       [material pairs,] (integer)<number_of_material_pairs>,\n"
            "           <material_pair> [,...]\n"
            "       [collision objects,] (integer)<number_of_collision_objects>,\n"
            "           (CollisionObject) <label> [,...]\n"
            "\n"
            "    <material_pair> ::= (str)<material1>, (str)<material2>, (ConstitutiveLaw<1D>)<const_law>\n"
            "       [, friction function, (ScalarFunction)<SF>, [, penetration ratio, (real)<penetration_ratio>]]\n"
            "\n\n"
            << std::endl);

        if (!HP.IsArg()) {
            /*
             * Exit quietly if nothing else is provided
             */
            throw NoErr(MBDYN_EXCEPT_ARGS);
        }
    }
    ConstLawType::Type VECLType(ConstLawType::VISCOELASTIC);
    typedef std::pair<std::string, std::string> MaterialPair;
    std::map<MaterialPair, ConstitutiveLaw1D*> pCL;
    std::map<MaterialPair, const BasicScalarFunction*> pSF;
    std::map<MaterialPair, doublereal> penetration_ratio;
    HP.IsKeyWord("material" "pairs");
    int N = HP.GetInt();
    for (int i; i < N; i++) {
        MaterialPair material_pair(std::make_pair(HP.GetValue(TypedValue::VAR_STRING).GetString(), HP.GetValue(TypedValue::VAR_STRING).GetString()));
        pCL[material_pair] = HP.GetConstLaw1D(VECLType);
        pSF[material_pair] = NULL;
        penetration_ratio[material_pair] = 0.0;
        if (HP.IsKeyWord("friction" "function")) {
            pSF[material_pair] = ParseScalarFunction(HP, pDM);
            if (HP.IsKeyWord("penetration" "ratio")) {
                penetration_ratio[material_pair] = HP.GetReal();
                if (material_pair.first == material_pair.second) {
                    silent_cout("Identical material pair penetration ratio overridden to become 0.5" << std::endl);
                }
            }
            if (material_pair.first == material_pair.second) {
                penetration_ratio[material_pair] = 0.5;
            }
        }
    }
    func_matrix = FCL::FuncMatrix();
    std::set<CollisionObjectData*> all_objects;
    std::set<CollisionObjectData*> objects;
    iNumRows = 0;
    iNumCols = 0;
    HP.IsKeyWord("collision" "objects");
    N = HP.GetInt();
    for (int i; i < N; i++) {
        CollisionObjectData* ob_data(collision_object_data[HP.GetInt()]);
        for (std::set<CollisionObjectData*>::iterator it = all_objects.begin();
            it != all_objects.end(); it++) {
            MaterialPair material_pair(std::make_pair(ob_data->material, (*it)->material));
            if (pCL.find(material_pair) == pCL.end()) {
                std::swap(material_pair.first, material_pair.second);
            }
            if (pCL.find(material_pair) != pCL.end() && ob_data->pNode != (*it)->pNode) {
                FCL::ObjectPair object_pair(std::make_pair(ob_data->pObject, (*it)->pObject));
                FCL::Func func(func_matrix.GetFunc(object_pair));
                if (func) {
                    objectpair_collision_map[object_pair] = new Collision(func,
                        pCL[material_pair], pSF[material_pair], penetration_ratio[material_pair],
                        ob_data, *it, &iNumRows, &iNumCols);
                } else {
                    std::swap(object_pair.first, object_pair.second);
                    objectpair_collision_map[object_pair] = new Collision(func_matrix.GetFunc(object_pair),
                        pCL[material_pair], pSF[material_pair], 1.0 - penetration_ratio[material_pair],
                        *it, ob_data, &iNumRows, &iNumCols);
                }
                objects.insert(ob_data);
                objects.insert(*it);
            }
        }
        all_objects.insert(ob_data);
    }
    collision_manager = new fcl::DynamicAABBTreeCollisionManager();
    for (std::set<CollisionObjectData*>::iterator it = objects.begin();
        it != objects.end(); it++) {
        nodes.insert((*it)->pNode);
        collision_manager->registerObject((*it)->pObject);
    }
    collision_manager->setup();
    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

CollisionWorld::~CollisionWorld(void)
{
    delete collision_manager;
}

void
CollisionWorld::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = iNumRows;
    *piNumCols = iNumCols;
}

int
CollisionWorld::iGetNumConnectedNodes(void) const
{
    return nodes.size();
}

void
CollisionWorld::GetConnectedNodes(std::vector<const Node*>& connectedNodes) const
{
    connectedNodes.resize(iGetNumConnectedNodes());
    integer i(0);
    for (std::set<const Node*>::const_iterator it = nodes.begin(); it != nodes.end(); it++) {
        connectedNodes[i++] = *it;
    }
}

void
CollisionWorld::Output(OutputHandler& OH) const
{
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
CollisionWorld::SetValue(DataManager *pDM,
    VectorHandler& X, VectorHandler& XP,
    SimulationEntity::Hints *ph)
{
    NO_OP;
}

void
CollisionWorld::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
    ss.str("");
    ss.clear();
    for (std::map<FCL::ObjectPair, Collision*>::const_iterator it = objectpair_collision_map.begin();
        it != objectpair_collision_map.end(); it++) {
        it->second->OutputAppend(ss);
        it->second->ClearAndSetTangents();
    }
}

SubVectorHandler& 
CollisionWorld::AssRes(SubVectorHandler& WorkVec,
    doublereal dCoef,
    const VectorHandler& XCurr, 
    const VectorHandler& XPrimeCurr)
{
    DEBUGCOUT("Entering CollisionWorld::AssRes()" << std::endl);
    for (std::map<FCL::ObjectPair, Collision*>::const_iterator it = objectpair_collision_map.begin();
        it != objectpair_collision_map.end(); it++) {
        it->second->ClearContacts();
    }
    collision_manager->update();
    collision_manager->collide(&objectpair_collision_map, CollisionFunction);
    WorkVec.ResizeReset(iNumRows);
    for (std::map<FCL::ObjectPair, Collision*>::const_iterator it = objectpair_collision_map.begin();
        it != objectpair_collision_map.end(); it++) {
        it->second->AssRes(WorkVec, dCoef, XCurr, XPrimeCurr);
    }
    return WorkVec;
}

VariableSubMatrixHandler& 
CollisionWorld::AssJac(VariableSubMatrixHandler& WorkMat,
    doublereal dCoef, 
    const VectorHandler& XCurr,
    const VectorHandler& XPrimeCurr)
{
    DEBUGCOUT("Entering CollisionWorld::AssJac()" << std::endl);
    FullSubMatrixHandler& WM = WorkMat.SetFull();
    WM.ResizeReset(iNumRows, iNumCols);
    for (std::map<FCL::ObjectPair, Collision*>::const_iterator it = objectpair_collision_map.begin();
        it != objectpair_collision_map.end(); it++) {
        it->second->AssJac(WorkMat, dCoef, XCurr, XPrimeCurr);
    }
    return WorkMat;
}

unsigned int
CollisionWorld::iGetNumPrivData(void) const
{
    // return number of private data
    return 0;
}

unsigned int
CollisionWorld::iGetPrivDataIdx(const char *s) const
{
    // parse string and compute index of requested private data

    // shouldn't get here until private data are defined
    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

doublereal
CollisionWorld::dGetPrivData(unsigned int i) const
{
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
    return 0;
}

void 
CollisionWorld::InitialWorkSpaceDim(
    integer* piNumRows,
    integer* piNumCols) const
{
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
//            "       Box, (real)<x_half_extent>, (real)<y_half_extent>, (real)<z_half_extent>\n"
//            "       | Capsule, (real)<radius>, (real)<height>\n"
//            "       | Cone, (real)<radius>, (real)<height>\n"
            "       | Sphere, (real)<radius>\n"
            "       | Plane\n"
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
    R = HP.GetRotRel(RF);
    Vec3 x(pNode->GetXCurr() + pNode->GetRCurr()*f);
    Mat3x3 r(pNode->GetRCurr()*R);
    fcl::Vec3f translate(x[0], x[1], x[2]);
    fcl::Matrix3f rotate(r.dGet(1,1),r.dGet(1,2),r.dGet(1,3),r.dGet(2,1),r.dGet(2,2),r.dGet(2,3),r.dGet(3,1),r.dGet(3,2),r.dGet(3,3));
    const TypedValue material(HP.GetValue(TypedValue::VAR_STRING));
    /*if (HP.IsKeyWord("box")) {
        const float x(HP.GetReal());
        const float y(HP.GetReal());
        const float z(HP.GetReal());
        FCL::CollisionGeometryPtr_t fcl_shape(new fcl::Box(2 * x, 2 * y, 2 * z));
        ob = new fcl::CollisionObject(fcl_shape, rotate, translate);
    } else if (HP.IsKeyWord("capsule")) {
        const float radius(HP.GetReal());
        const float height(HP.GetReal());
        FCL::CollisionGeometryPtr_t fcl_shape(new fcl::Capsule(radius, height));
        ob = new fcl::CollisionObject(fcl_shape, rotate, translate);
    } else if (HP.IsKeyWord("cone")) {
        const float radius(HP.GetReal());
        const float height(HP.GetReal());
        FCL::CollisionGeometryPtr_t fcl_shape(new fcl::Cone(radius, height));
        ob = new fcl::CollisionObject(fcl_shape, rotate, translate);
    } else */if (HP.IsKeyWord("plane")) {
        FCL::CollisionGeometryPtr_t fcl_shape(new fcl::Plane(0., 0., 1., 0.));
        ob = new fcl::CollisionObject(fcl_shape, rotate, translate);
    } else if (HP.IsKeyWord("sphere")) {
        const float radius(HP.GetReal());
        FCL::CollisionGeometryPtr_t fcl_shape(new fcl::Sphere(radius));
        ob = new fcl::CollisionObject(fcl_shape, rotate, translate);
    } else {
        silent_cerr("collision object(" << GetLabel() << "): a valid shape is expected at line " << HP.GetLineData() << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
    collision_object_data[uLabel] = new CollisionObjectData(pNode, ob, material.GetString());
    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

CollisionObject::~CollisionObject(void)
{
    delete ob;
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
    if (ob->getNodeType() != fcl::GEOM_PLANE) {
        Vec3 x(pNode->GetXCurr() + pNode->GetRCurr() * f);
        Mat3x3 r(pNode->GetRCurr() * R);
        ob->setTransform(
            fcl::Matrix3f(r.dGet(1,1),r.dGet(1,2),r.dGet(1,3),r.dGet(2,1),r.dGet(2,2),r.dGet(2,3),r.dGet(3,1),r.dGet(3,2),r.dGet(3,3)),
            fcl::Vec3f(x[0], x[1], x[2]));
        ob->computeAABB();
    }    
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

