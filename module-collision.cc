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
#include <queue>
#include "rodj.h"
#include <limits>
#include "module-collision.h"

Impact::Impact(const DofOwner* pDO, 
    const ConstitutiveLaw1D* pCL, const BasicScalarFunction* pSFTmp, const doublereal penetrationTmp,
    const StructNode* pN1, const StructNode* pN2)
: Elem(-1, flag(0)),
super(-1, pDO, pCL, pN1, pN2, Zero3, Zero3, 1.0, flag(0)),
iNumRowsNode(6),
iNumColsNode(6),
iNumDofs(2),
pSF(pSFTmp),
penetration(penetrationTmp)
{
    dElle = 0.0;
}

doublereal
Impact::dCalcEpsilon(void)
{
    return dElle; //depth of impact
}

unsigned int
Impact::iGetNumDof(void) const
{
    return iNumDofs;
}

void
Impact::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = 2 * iNumRowsNode + iGetNumDof();
    *piNumCols = 2 * iNumColsNode + iGetNumDof();
}

int
Impact::ContactsSize(void)
{
    return Arm2_vector.size();
}

void
Impact::ClearContacts(void)
{
    Arm2_vector.clear();
    f1_vector.clear();
    f2_vector.clear();
    state_vector.clear();
    Ft_vector.clear();
}

void
Impact::PushBackContact(const btVector3& point1, const btVector3& point2)
{
    const Vec3 p1(point1[0], point1[1], point1[2]);
    const Vec3 p2(point2[0], point2[1], point2[2]);
    //printf("push pts: (%f, %f, %f), (%f, %f, %f)\n", p1(1), p1(2), p1(3), p2(1), p2(2), p2(3));
    Arm2_vector.push_back(dynamic_cast<const StructNode *>(pNode2)->GetRCurr().Transpose()*(p1 + (p2 - p1) * penetration - pNode2->GetXCurr()));
    f1_vector.push_back(dynamic_cast<const StructNode *>(pNode1)->GetRCurr().Transpose()*(p1 - pNode1->GetXCurr()));
    f2_vector.push_back(dynamic_cast<const StructNode *>(pNode2)->GetRCurr().Transpose()*(p2 - pNode2->GetXCurr()));
    state_vector.push_back(UNDEFINED);
    Ft_vector.push_back(Zero3);
    Rv_vector.push_back(Eye3);
}

void
Impact::SetFirstReactionIndex(const int i)
{
    iFirstReactionIndex = i;
}

bool
Impact::SetOffsets(const int i)
{
    index = i;
    Arm2 = Arm2_vector[i];
    f1 = f1_vector[i];
    f2 = f2_vector[i];
    state = state_vector[i];
    Ft_old = Ft_vector[i];
    Mat3x3 R1(dynamic_cast<const StructNode *>(pNode1)->GetRCurr());
    Mat3x3 R2(dynamic_cast<const StructNode *>(pNode2)->GetRCurr());
    Vec3 p2(pNode2->GetXCurr() + R2 * f2);
    Vec3 p1(pNode1->GetXCurr() + R1 * f1);
    return (p2 - p1).Dot(pNode2->GetXCurr() - pNode1->GetXCurr()) < 0.0
        && std::numeric_limits<doublereal>::epsilon() < (p2 - p1).Norm(); // impact?
}

void
Impact::SaveState(void)
{
    state_vector[index] = state;
    Ft_vector[index] = Ft_old;
}

void
Impact::PutRv(Mat3x3 Rv)
{
    Rv_vector[index] = Rv;
}

Mat3x3
Impact::GetRv(const int i)
{
    return Rv_vector[index];
}

bool
Impact::NoSlip(const int i)
{
    return state_vector[i] == NO_SLIP;
}

VariableSubMatrixHandler&
Impact::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
    //printf("Entering Impact::AssJac()\n");
	DEBUGCOUT("Entering Impact::AssJac()" << std::endl);
    integer iNumRows = 0;
    integer iNumCols = 0;
    WorkSpaceDim(&iNumRows, &iNumCols);
    WorkMat.FullSubMatrixHandler::ResizeReset(iNumRows, iNumCols);
    FullSubMatrixHandler& WM = super::AssJac(WorkMat, dCoef, XCurr, XPrimeCurr);
    WM.Resize(iNumRows, iNumCols);
    if (pSF == 0 || state != NO_SLIP) {
        //printf("Impact::AssJac() other\n");
	    for (int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		    WM.PutCoef(2 * iNumRowsNode + iCnt, 2 * iNumColsNode + iCnt, 1.0);
	    }
    } else {
        //printf("Impact::AssJac() NO_SLIP\n");
        const StructNode* pStructNode1(dynamic_cast<const StructNode *>(pNode1));
        const StructNode* pStructNode2(dynamic_cast<const StructNode *>(pNode2));
	    Vec3 R_Arm2(pStructNode2->GetRCurr() * Arm2);
	    Vec3 R_Arm1_calc(pNode2->GetXCurr() + R_Arm2 - pNode1->GetXCurr());

        Vec3 dp(v / dElle);
        //printf("dp: %f, %f, %f\n", dp(1), dp(2), dp(3));
        Vec3 cz(-dp(2), dp(1), 0.0); //cross product of dp with (0, 0, 1)
        //printf("cz: %f, %f, %f\n", cz(1), cz(2), cz(3));
        Mat3x3 R_Rv(Eye3);
        if (cz.Norm()) {
            cz /= cz.Norm();
            doublereal u(cz(1));
            doublereal v(cz(2)); //doublereal w = cz(3); = 0.0
            //printf("cz after Norm: %f, %f, %f\n", cz(1), cz(2), cz(3));
            //printf("u, v: %f, %f\n", u, v);
            doublereal rcos(dp(3));
            doublereal rsin(sqrt(dp(1) * dp(1) + dp(2) * dp(2)));
            //printf("rcos, rsin:, %f, %f\n", rcos, rsin);
            R_Rv = Mat3x3(
                rcos + u*u*(1-rcos), v*u*(1-rcos),       -v * rsin,
                u*v*(1-rcos),        rcos + v*v*(1-rcos), u * rsin,
                v * rsin,           -u * rsin,            rcos
            );    
        }
        Vec3 FTmp(R_Rv * (F_reaction*dCoef));
        Vec3 Tmp1_1(R_Rv.GetVec(1).Cross(R_Arm1_calc));
        Vec3 Tmp1_2(R_Rv.GetVec(2).Cross(R_Arm1_calc));
        Vec3 Tmp2_1(R_Arm2.Cross(R_Rv.GetVec(1)));
        Vec3 Tmp2_2(R_Arm2.Cross(R_Rv.GetVec(2)));
        for(int iCnt = 1; iCnt <= 3; iCnt++) {
            doublereal d(R_Rv.dGet(iCnt, 1));
            WM.DecCoef(0 * iNumRowsNode + iCnt, 2 * iNumColsNode + 1,     d);
            WM.IncCoef(1 * iNumRowsNode + iCnt, 2 * iNumColsNode + 1,     d); 
            WM.DecCoef(2 * iNumRowsNode + 1   , 0 * iNumColsNode + iCnt,  d); 
            WM.IncCoef(2 * iNumRowsNode + 1   , 1 * iNumColsNode + iCnt,  d); 

            d = R_Rv.dGet(iCnt, 2);
            WM.DecCoef(0 * iNumRowsNode + iCnt, 2 * iNumColsNode + 2,      d); 
            WM.IncCoef(1 * iNumRowsNode + iCnt, 2 * iNumColsNode + 2,      d); 
            WM.DecCoef(2 * iNumRowsNode + 2   , 0 * iNumColsNode + iCnt,   d); 
            WM.IncCoef(2 * iNumRowsNode + 2   , 1 * iNumColsNode + iCnt,   d); 

            d = Tmp1_1.dGet(iCnt);
            WM.IncCoef(0 * iNumRowsNode+3+iCnt, 2 * iNumColsNode + 1,      d); 
            WM.IncCoef(2 * iNumRowsNode + 1   , 0 * iNumColsNode+3+iCnt,   d); 

            d = Tmp1_2.dGet(iCnt);
            WM.IncCoef(0 * iNumRowsNode+3+iCnt, 2 * iNumColsNode + 2,      d); 
            WM.IncCoef(2 * iNumRowsNode + 2   , 0 * iNumColsNode+3+iCnt,   d); 

            d = Tmp2_1.dGet(iCnt);
            WM.IncCoef(1 * iNumRowsNode+3+iCnt, 2 * iNumColsNode + 1,      d); 
            WM.IncCoef(2 * iNumRowsNode + 1   , 1 * iNumColsNode+3+iCnt,   d); 

            d = Tmp2_2.dGet(iCnt);
            WM.IncCoef(1 * iNumRowsNode+3+iCnt, 2 * iNumColsNode + 2,      d); 
            WM.IncCoef(2 * iNumRowsNode + 2   , 1 * iNumColsNode+3+iCnt,   d);
        }
        Mat3x3 FTmpCross(MatCross, FTmp);
        WM.Add(0 * iNumRowsNode + 1, 0 * iNumColsNode + 4,    FTmpCross); 
        WM.Add(0 * iNumRowsNode + 4, 0 * iNumColsNode + 1,   -FTmpCross); 
        WM.Add(0 * iNumRowsNode + 4, 0 * iNumColsNode + 4, Mat3x3(MatCrossCross, R_Arm1_calc, FTmp)); 
        WM.Add(0 * iNumRowsNode + 4, 1 * iNumColsNode + 1,    FTmpCross); 
        WM.Add(1 * iNumRowsNode + 1, 1 * iNumColsNode + 4,   -FTmpCross); 

        Mat3x3 MTmp(MatCrossCross, FTmp, R_Arm2);
        WM.Add(0 * iNumRowsNode + 4, 1 * iNumColsNode + 4,      -MTmp); 
        WM.Add(1 * iNumRowsNode + 4, 1 * iNumColsNode + 4,      -MTmp); 
        WM.Add(1 * iNumRowsNode + 4, 0 * iNumColsNode + 4, Mat3x3(MatCrossCross, -R_Arm2, FTmp)); 
        PutRv(pStructNode1->GetRCurr().Transpose()*R_Rv);
    }
    /*
    printf("WM:\n");
    for (int r = 1; r <= 2 * iNumRowsNode + iGetNumDof(); r++) {
        printf("Row %d:", r);
        for (int c = 1; c <= 2 * iNumColsNode + iGetNumDof(); c++) {
            printf(", %f", WM.dGetCoef(r, c));
        }
        printf("\n");
    }
    */

    return WorkMat;
}

SubVectorHandler& 
Impact::AssRes(SubVectorHandler& WorkVec,
    doublereal dCoef,
    const VectorHandler& XCurr, 
    const VectorHandler& XPrimeCurr)
{
    //printf("Entering Impact::AssRes()\n");
    DEBUGCOUT("Entering Impact::AssRes()" << std::endl);
    integer iNumRows(0);
    integer iNumCols(0);
    WorkSpaceDim(&iNumRows, &iNumCols);
    WorkVec.ResizeReset(iNumRows);
    super::AssRes(WorkVec, dCoef, XCurr, XPrimeCurr);
    WorkVec.Resize(iNumRows);
    //printf("XCurr: %f, %f\n", XCurr(iFirstReactionIndex + 1), XCurr(iFirstReactionIndex + 2));
    if (pSF == 0 || dElle == 0 || state == MISS) {
        state = MISS;
        //printf("Impact::AssRes() dElle == 0\n");
	    for (int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		    WorkVec.PutCoef(2 * iNumRowsNode + iCnt, -XCurr(iFirstReactionIndex + iCnt));
	    }
        return WorkVec;
    }
    doublereal F_rod = GetF();
    if (F_rod < 0.0 || state == NO_RESISTANCE) {
        state = NO_RESISTANCE;
        //printf("Impact::AssRes() F_rod < 0.0\n");
	    for (int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		    WorkVec.PutCoef(2 * iNumRowsNode + iCnt, -XCurr(iFirstReactionIndex + iCnt));
	    }
        return WorkVec;
    }
    const StructNode* pStructNode1(dynamic_cast<const StructNode *>(pNode1));
    const StructNode* pStructNode2(dynamic_cast<const StructNode *>(pNode2));
    Vec3 G1(Zero3);
    Vec3 B1(Zero3);
    if (pStructNode1->GetVCurr().Dot() + pStructNode1->GetWCurr().Dot() + pStructNode1->GetXPPCurr().Dot() + pStructNode1->GetWPCurr().Dot() != 0.0) {
        G1 = dynamic_cast<const DynamicStructNode *>(pNode1)->GetGCurr();
        B1 = dynamic_cast<const DynamicStructDispNode *>(pNode1)->GetBCurr();
    }
    Vec3 G2(Zero3);
    Vec3 B2(Zero3);
    if (pStructNode2->GetVCurr().Dot() + pStructNode2->GetWCurr().Dot() + pStructNode2->GetXPPCurr().Dot() + pStructNode2->GetWPCurr().Dot() != 0.0) {
        G2 = dynamic_cast<const DynamicStructNode *>(pNode2)->GetGCurr();
        B2 = dynamic_cast<const DynamicStructDispNode *>(pNode2)->GetBCurr();
    }
    Mat3x3 R1(pStructNode1->GetRCurr());
    Mat3x3 R2(pStructNode2->GetRCurr());
	Vec3 R_Arm2(R2*Arm2);
	Vec3 R_Arm1_calc(pNode2->GetXCurr() + R_Arm2 - pNode1->GetXCurr());
    Vec3 dFt1((R_Arm1_calc.Cross(B1) + G1).Cross(R_Arm1_calc));
    Vec3 dFt2((R_Arm2.Cross(B2) + G2).Cross(R_Arm2));
    Vec3 Ft((dFt2 - dFt1) / dCoef);
    Vec3 V(pNode2->GetVCurr() - R_Arm2.Cross(pStructNode2->GetWCurr()) - pNode1->GetVCurr() + R_Arm1_calc.Cross(pStructNode1->GetWCurr()));
    Vec3 unit_normal(pNode2->GetXCurr() + R2 * f2 - pNode1->GetXCurr() + R1 * f1);
    unit_normal /= unit_normal.Norm();
    Ft = Ft - unit_normal*Ft.Dot(unit_normal);
    V = V - unit_normal*V.Dot(unit_normal);
    const doublereal mu((*pSF)(V.Norm()));
    /*
    printf("B1: %f, %f, %f\n", B1(1), B1(2), B1(3));
    printf("B2: %f, %f, %f\n", B2(1), B2(2), B2(3));
    printf("G1: %f, %f, %f\n", G1(1), G1(2), G1(3));
    printf("G2: %f, %f, %f\n", G2(1), G2(2), G2(3));
    printf("dFt1: %f, %f, %f\n", dFt1(1), dFt1(2), dFt1(3));
    printf("dFt2: %f, %f, %f\n", dFt2(1), dFt2(2), dFt2(3));
    printf("Ft: %f, %f, %f before\n", Ft(1), Ft(2), Ft(3));
    printf("Ft_old: %f, %f, %f before\n", Ft_old(1), Ft_old(2), Ft_old(3));
    */
    if (state != NO_SLIP && 0 <= Ft.Dot(Ft_old) && ((mu * F_rod < Ft.Norm() && 1e-9 < V.Norm() * dCoef) || state == SLIP)) {
        state = SLIP;
        //printf("Impact::AssRes() SLIP\n");
        Ft = Ft * (mu * F_rod / Ft.Norm());
        //printf("Ft: %f, %f, %f; mu * F_rod: %f\n", Ft(1), Ft(2), Ft(3), mu * F_rod);
        WorkVec.Add(1, Ft);
        WorkVec.Add(4, R_Arm1_calc.Cross(Ft));
        WorkVec.Sub(7, Ft);
        WorkVec.Sub(10, R_Arm2.Cross(Ft));
	    for (int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		    WorkVec.PutCoef(2 * iNumRowsNode + iCnt, -XCurr(iFirstReactionIndex + iCnt));
	    }
    } else {
        F_reaction = Zero3;
        if (state == NO_SLIP) {
            F_reaction.Put(1, XCurr(iFirstReactionIndex+1));
            F_reaction.Put(2, XCurr(iFirstReactionIndex+2));
        } else {
            //printf("Reset XCurr ???\n");
            state = NO_SLIP;
        }
        //printf("Impact::AssRes() else\n");
        Vec3 dp(v / dElle);
        //printf("dp: %f, %f, %f\n", dp(1), dp(2), dp(3));
        Vec3 cz(-dp(2), dp(1), 0.0); //cross product of dp with (0, 0, 1)
        //printf("cz: %f, %f, %f\n", cz(1), cz(2), cz(3));
        Mat3x3 R_Rv(Eye3);
        if (cz.Norm()) {
            cz /= cz.Norm();
            doublereal u(cz(1));
            doublereal v(cz(2)); //doublereal w = cz(3); = 0.0
            //printf("cz after Norm: %f, %f, %f\n", cz(1), cz(2), cz(3));
            //printf("u, v: %f, %f\n", u, v);
            doublereal rcos(dp(3));
            doublereal rsin(sqrt(dp(1) * dp(1) + dp(2) * dp(2)));
            //printf("rcos, rsin:, %f, %f\n", rcos, rsin);
            R_Rv = Mat3x3(
                rcos + u*u*(1-rcos), v*u*(1-rcos),       -v * rsin,
                u*v*(1-rcos),        rcos + v*v*(1-rcos), u * rsin,
                v * rsin,           -u * rsin,            rcos
            );    
        }
        Vec3 FTmp(R_Rv*F_reaction);
        WorkVec.Add(1, FTmp);
        WorkVec.Add(4, R_Arm1_calc.Cross(FTmp));
        WorkVec.Sub(7, FTmp);

        ASSERT(dCoef != 0.);
        WorkVec.PutCoef(2 * iNumRowsNode + 1, - R_Rv.GetVec(1).Dot(V));
        WorkVec.PutCoef(2 * iNumRowsNode + 2, - R_Rv.GetVec(2).Dot(V));
        //printf("WorkVec Reaction: %f, %f\n", WorkVec(2 * iNumRowsNode + 1), WorkVec(2 * iNumRowsNode + 2));
        //printf("V possible: %f, %f. dCoef: %f\n", - R_Rv.GetVec(1).Dot(V), - R_Rv.GetVec(2).Dot(V), dCoef);
    }
    //printf("R_Arm1_calc: %f, %f, %f\n", R_Arm1_calc(1), R_Arm1_calc(2), R_Arm1_calc(3));

    ASSERTMSGBREAK(state != UNDEFINED, "state is UNDEFINED");
    Ft_old = Ft;
    SaveState();
    return WorkVec;
}

void
Impact::AccumulatorReset(void)
{
    //sum_F = Zero3;
    //sum_M = Zero3;
}

std::ostream&
Impact::OutputAppend(std::ostream& out) const {
    //out << " " << sqrt(sum_F.Dot()) << " " << sqrt(sum_M.Dot()) << " " << sum_F << " " << sum_M;
    //ConstitutiveLaw1DOwner::OutputAppend(out);
}

class Universe {
public:
    std::map<unsigned, btCollisionObject*> label_object_map;
    std::map<btCollisionObject*, const StructNode*> object_node_map;
} universe;

// CollisionWorld: begin

CollisionWorld::CollisionWorld(
    unsigned uLabel, const DofOwner *pDO,
    DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
iNumRowsNode(6),
iNumColsNode(6),
iNumDofsImpact(8)
{
    printf("Entering CollisionWorld::CollisionWorld()\n");
    // help
    if (HP.IsKeyWord("help")) {
        silent_cout(
            "\n"
            "Module:     Collision\n"
            "\n"
            "    This element implements a collision world\n"
            "\n"
            "    collision world,\n"
            "       <number_of_collision_pairs>,\n"
            "           <pair> [,...]\n"
            "\n"
            "   <pair> ::= (CollisionObject) <label_1>, (CollisionObject) <label_2>, (ConstitutiveLaw<1D>) <const_law>\n"
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
    N = HP.GetInt();
    ConstLawType::Type VECLType(ConstLawType::VISCOELASTIC);
    ConstLawType::Type VCLType(ConstLawType::VISCOUS);
    std::set<NodeObjectPair> node_object_pairs;
    std::set<btCollisionObject*> objects;
    unsigned k(0);
    for (int i; i < N; i++) {
        btCollisionObject* ob1(universe.label_object_map[HP.GetInt()]);
        btCollisionObject* ob2(universe.label_object_map[HP.GetInt()]);
        objects.insert(ob1);
        objects.insert(ob2);
        ObjectPair object_pair(std::make_pair(ob1, ob2));
        ConstitutiveLaw1D* pCL(HP.GetConstLaw1D(VECLType));
        const BasicScalarFunction* pSF(0);
        doublereal penetration(0.0);
        if (HP.IsKeyWord("friction" "function")) {
            pSF = ParseScalarFunction(HP, pDM);
            penetration = HP.GetReal();
        }
        objectpair_impact_map[object_pair] = new Impact(pDO, pCL, pSF, penetration,
            universe.object_node_map[ob1], universe.object_node_map[ob2]);
        objectpair_index_map[object_pair] = k++;
        node_object_pairs.insert(std::make_pair(universe.object_node_map[ob1], ob1));
        node_object_pairs.insert(std::make_pair(universe.object_node_map[ob2], ob2));
    }
    for (std::set<btCollisionObject*>::iterator it = objects.begin();
        it != objects.end(); it++) {
        world->addCollisionObject((*it));
    }
    k = 0;
    nodes_vector.resize(node_object_pairs.size());
    for (std::set<NodeObjectPair>::iterator it = node_object_pairs.begin();
        it != node_object_pairs.end(); it++) {
        nodes_vector[k] = (*it).first;
        object_index_map[(*it).second] = k++;
    }

    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

CollisionWorld::~CollisionWorld(void)
{
    for (std::map<ObjectPair, Impact*>::const_iterator it = objectpair_impact_map.begin();
        it != objectpair_impact_map.end(); it++) {
        delete it->second;
    }
    objectpair_impact_map.clear();
    delete world;
    delete broadphase;
    delete dispatcher;
    delete configuration;
}

unsigned int
CollisionWorld::iGetNumDof(void) const
{
    return iNumDofsImpact * objectpair_impact_map.size();
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
CollisionWorld::Output(OutputHandler& OH) const
{
    //printf("Entering CollisionWorld::OutputHandler()\n");
    if (fToBeOutput()) {
        if ( OH.UseText(OutputHandler::LOADABLE) ) {
            std::ostream& os = OH.Loadable();
            os << GetLabel();
            for (std::map<ObjectPair, Impact*>::const_iterator it = objectpair_impact_map.begin();
                it != objectpair_impact_map.end(); it++) {
                std::set<ObjectPair>::const_iterator collision_it = collisions.find(it->first);
                if (collision_it != collisions.end()) {
                    it->second->OutputAppend(os);
                } else {
                    for (int i = 0; i < 8; i++) {
                    os << " " << 0.0;
                    }
                }
            }
            os << std::endl;
        }
    }
}

void
CollisionWorld::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    //printf("Entering CollisionWorld::WorkSpaceDim()\n");
    *piNumRows = iNumRowsNode * nodes_vector.size() + iGetNumDof();
    *piNumCols = iNumColsNode * nodes_vector.size() + iGetNumDof();
}

VariableSubMatrixHandler& 
CollisionWorld::AssJac(VariableSubMatrixHandler& WorkMat,
    doublereal dCoef, 
    const VectorHandler& XCurr,
    const VectorHandler& XPrimeCurr)
{
    //printf("Entering CollisionWorld::AssJac()\n");
    //printf("XCurr.iGetSize(), XPrimeCurr.iGetSize()= %d, %d\n", XCurr.iGetSize(), XPrimeCurr.iGetSize());
    DEBUGCOUT("Entering CollisionWorld::AssJac()" << std::endl);
    //WorkMat.SetNullMatrix();
    //return WorkMat;
    integer iNumRows(0);
    integer iNumCols(0);
    WorkSpaceDim(&iNumRows, &iNumCols);
    const integer iFirstReactionIndex(iGetFirstIndex());
    const integer iFirstRowDof(iNumRows - iGetNumDof());
    const integer iFirstColDof(iNumCols - iGetNumDof());
    FullMatrixHandler MyMat = FullMatrixHandler();
    MyMat.ResizeReset(iNumRows, iNumCols);
    for (std::map<ObjectPair, Impact*>::const_iterator it = objectpair_impact_map.begin();
        it != objectpair_impact_map.end(); it++) {
        ObjectPair object_pair(it->first);
        int i1(object_index_map[object_pair.first]);
        int i2(object_index_map[object_pair.second]);
        int k(objectpair_index_map[object_pair]);
        int iRow1(iNumRowsNode * i1);
        int iRow2(iNumRowsNode * i2);
        /*
        printf("MyMat:\n");
        for (int r = 1; r <= 2 * iNumRowsNode + iNumDofsImpact; r++) {
            printf("Row %d:", r);
            for (int c = 1; c <= 2 * iNumColsNode + iNumDofsImpact; c++) {
                printf(", %f", MyMat.dGetCoef(r, c));
            }
            printf("\n");
        }
        */
        for (int i = 0; i < it->second->ContactsSize(); i++) {
            int iRowDof((iNumRows - iGetNumDof()) + (iNumDofsImpact * k) + (2 * i));
            int iColDof((iNumCols - iGetNumDof()) + (iNumDofsImpact * k) + (2 * i));
            if (it->second->SetOffsets(i)) {
                it->second->SetFirstReactionIndex(iFirstReactionIndex + (iNumDofsImpact * k) + (2 * i));
                FullSubMatrixHandler& WM(it->second->AssJac(WorkMat, dCoef, XCurr, XPrimeCurr));
                // Requires iNumRowsNode = 6, otherwise more nested 'for' loops must be written
                for (int r = 1; r <= 6; r++) {
                    for (int c = 1; c <= 6; c++) {
                        MyMat.IncCoef(iRow1 + r, (iNumColsNode * i1) + c, WM.dGetCoef(r, c));
                        MyMat.IncCoef(iRow1 + r, (iNumColsNode * i2) + c, WM.dGetCoef(r, 6 + c));
                        MyMat.IncCoef(iRow2 + r, (iNumColsNode * i1) + c, WM.dGetCoef(6 + r, c));
                        MyMat.IncCoef(iRow2 + r, (iNumColsNode * i2) + c, WM.dGetCoef(6 + r, 6 + c));
                    }
                    for (int c = 1; c <= (iNumColsNode - 6); c++) {
                        MyMat.IncCoef(iRow1 + r, (iNumColsNode * i1) + 6 + c, WM.dGetCoef(r, 12 + c));
                        MyMat.IncCoef(iRow1 + r, (iNumColsNode * i2) + 6 + c, WM.dGetCoef(r, 6 + iNumColsNode + c));
                        MyMat.IncCoef(iRow2 + r, (iNumColsNode * i1) + 6 + c, WM.dGetCoef(6 + r, 12 + c));
                        MyMat.IncCoef(iRow2 + r, (iNumColsNode * i2) + 6 + c, WM.dGetCoef(6 + r, 6 + iNumColsNode + c));
                    }
                    for (int c = 1; c <= 2; c++) {
                        MyMat.IncCoef(iRow1 + r, iColDof + c, WM.dGetCoef(r, 2 * iNumColsNode + c));
                        MyMat.IncCoef(iRow2 + r, iColDof + c, WM.dGetCoef(6 + r, 2 * iNumColsNode + c));
                    }
                }
                for (int r = 1; r <= 2; r++) {
                    for (int c = 1; c <= 6; c++) {
                        MyMat.IncCoef(iRowDof + r, (iNumColsNode * i1) + c, WM.dGetCoef(2 * iNumRowsNode + r, c));
                        MyMat.IncCoef(iRowDof + r, (iNumColsNode * i2) + c, WM.dGetCoef(2 * iNumRowsNode + r, 6 + c));
                    }
                    for (int c = 1; c <= (iNumColsNode - 6); c++) {
                        MyMat.IncCoef(iRowDof + r, (iNumColsNode * i1) + 6 + c, WM.dGetCoef(2 * iNumRowsNode + r, 12 + c));
                        MyMat.IncCoef(iRowDof + r, (iNumColsNode * i2) + 6 + c, WM.dGetCoef(2 * iNumRowsNode + r, 6 + iNumColsNode + c));
                    }
                    for (int c = 1; c <= 2; c++) {
                        MyMat.IncCoef(iRowDof + r, iColDof + c, WM.dGetCoef(2 * iNumRowsNode + r, 2 * iNumColsNode + c));
                    }
                }
            } else {
                for (int rc = 1; rc <= 2; rc++) {
                    MyMat.IncCoef(iRowDof + rc, iColDof + rc, 1.0);
                }
            }
            /*
            printf("MyMat:\n");
            for (int r = 1; r <= 2 * iNumRowsNode + iNumDofsImpact; r++) {
                printf("Row %d:", r);
                for (int c = 1; c <= 2 * iNumColsNode + iNumDofsImpact; c++) {
                    printf(", %f", MyMat.dGetCoef(r, c));
                }
                printf("\n");
            }
            */
        }
        /*
        printf("NoSlip\n");
        for (int i = 0; i < it->second->ContactsSize(); i++) {
            for (int j = 0; j < it->second->ContactsSize(); j++) {
                if (i != j && it->second->NoSlip(i) && it->second->NoSlip(j)) {
                    int iRowDof((iNumRows - iGetNumDof()) + (iNumDofsImpact * k) + (2 * i));
                    int iColDof((iNumCols - iGetNumDof()) + (iNumDofsImpact * k) + (2 * j));
                    Mat3x3 dDi_dDj(- it->second->GetRv(i).Transpose() * it->second->GetRv(j));
                    printf("%d, %d\n", i, j);
                    for (int r = 1; r <=2; r++) {
                        for (int c = 1; c <= 2; c++) {
                            MyMat.IncCoef(iRowDof + r, iColDof + c, dDi_dDj(r, c));
                        }
                    }
                }
            }
        }
        printf("MyMat:\n");
        for (int r = 1; r <= 2 * iNumRowsNode + iNumDofsImpact; r++) {
            printf("Row %d:", r);
            for (int c = 1; c <= 2 * iNumColsNode + iNumDofsImpact; c++) {
                printf(", %f", MyMat.dGetCoef(r, c));
            }
            printf("\n");
        }
        */
        for (int i = it->second->ContactsSize(); i < iNumDofsImpact / 2; i++) {
            int iRowDof((iNumRows - iGetNumDof()) + (iNumDofsImpact * k) + (2 * i));
            int iColDof((iNumCols - iGetNumDof()) + (iNumDofsImpact * k) + (2 * i));
            for (int rc = 1; rc <= 2; rc++) {
                MyMat.IncCoef(iRowDof + rc, iColDof + rc, 1.0);
            }
        }
    }
    WorkMat.SetNullMatrix();
    FullSubMatrixHandler& WM = WorkMat.SetFull();
    WM.ResizeReset(iNumRows, iNumCols);
    for (int i = 0; i < nodes_vector.size(); i++) {
        const integer iFirstPositionIndex(nodes_vector[i]->iGetFirstPositionIndex());
        const integer iFirstMomentumIndex(nodes_vector[i]->iGetFirstMomentumIndex());
        for (int r = 1; r <= 6; r++) {
            WM.PutRowIndex((6 * i) + r, iFirstMomentumIndex + r);
        }
        for (int c = 1; c <= 6; c++) {
            WM.PutColIndex((iNumColsNode * i) + c, iFirstPositionIndex + c);
        }
        for (int c = 1; c <= (iNumColsNode - 6); c++) {
            WM.PutColIndex((iNumColsNode * i) + 6 + c, iFirstMomentumIndex + c);
        }
    }
    for (int r = 1; r <= iNumRows; r++) {
        for (int c = 1; c <= iNumCols; c++) {
            WM.PutCoef(r, c, MyMat(r, c));
        }
    }
    for (int i = 0; i < objectpair_index_map.size(); i++) {
        for (int increment = (iNumDofsImpact * i) + 1; increment <= (iNumDofsImpact * i) + iNumDofsImpact; increment++) {
            WM.PutRowIndex(iFirstRowDof + increment, iFirstReactionIndex + increment);
            WM.PutColIndex(iFirstColDof + increment, iFirstReactionIndex + increment);
        }
    }
    for (std::map<ObjectPair, unsigned int>::const_iterator it = objectpair_index_map.begin();
        it != objectpair_index_map.end(); it++) {
        std::set<ObjectPair>::const_iterator collision_it = collisions.find(it->first);
        if (collision_it == collisions.end()) {
            //printf("CollisionWorld::AssJac() Miss!\n");
            for (int increment = (iNumDofsImpact * (it->second)) + 1; increment <= (iNumDofsImpact * (it->second)) + iNumDofsImpact; increment++) {
                WM.PutCoef(iFirstRowDof + increment, iFirstColDof + increment, 1.0);
            }
        }
    }
    /*
    printf("End:\n");
    for (int r = 1; r <= iNumRows; r++) {
        printf("Row %i:", r);
        for (int c = 1; c <= iNumCols; c++) {
            printf(", %f", WM.dGetCoef(r, c));
        }
        printf("\n");
    }
    */
    return WorkMat;
}

void
CollisionWorld::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
    //printf("Entering CollisionWorld::AfterConvergence()\n");
    for (std::map<ObjectPair, Impact*>::const_iterator it = objectpair_impact_map.begin();
        it != objectpair_impact_map.end(); it++) {
        it->second->ClearContacts();
    }
    world->performDiscreteCollisionDetection();
    for (int i = 0; i < dispatcher->getNumManifolds() ; i++) {
        btPersistentManifold* manifold(dispatcher->getManifoldByIndexInternal(i));
        ObjectPair object_pair = std::make_pair(
            (btCollisionObject*)(manifold->getBody0()), (btCollisionObject*)(manifold->getBody1()));
        bool swapped(false);
        if (objectpair_impact_map.count(object_pair) == 0) {
            std::swap(object_pair.first, object_pair.second);
            swapped = true;
        }
        if (objectpair_impact_map.count(object_pair)) {
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
                objectpair_impact_map[object_pair]->PushBackContact(p1, p2);
                //printf("Contact!!\n");
            }
        }
    }
}

SubVectorHandler& 
CollisionWorld::AssRes(SubVectorHandler& WorkVec,
    doublereal dCoef,
    const VectorHandler& XCurr, 
    const VectorHandler& XPrimeCurr)
{
    //printf("Entering CollisionWorld::AssRes()\n");
    DEBUGCOUT("Entering CollisionWorld::AssRes()" << std::endl);

    // compute contact forces and write them in contribution to residual

    for (std::map<ObjectPair, Impact*>::const_iterator it = objectpair_impact_map.begin();
        it != objectpair_impact_map.end(); it++) {
        it->second->AccumulatorReset();
    }
    integer iNumRows(0);
    integer iNumCols(0);
    WorkSpaceDim(&iNumRows, &iNumCols);
    const integer iFirstReactionIndex(iGetFirstIndex());
    const integer iFirstRowDof(iNumRows - iGetNumDof());
    MyVectorHandler MyVec = MyVectorHandler();
    MyVec.ResizeReset(iNumRows);
    bool ChangeJac(false);
    collisions.clear();
    for (std::map<ObjectPair, Impact*>::const_iterator it = objectpair_impact_map.begin();
        it != objectpair_impact_map.end(); it++) {
        ObjectPair object_pair(it->first);
        int k(objectpair_index_map[object_pair]);
        int iRow1((iNumRowsNode * object_index_map[object_pair.first]));
        int iRow2((iNumRowsNode * object_index_map[object_pair.second]));
        for (int i = 0; i < it->second->ContactsSize(); i++) {
            //printf("i = %d", i);
            int iRowDof((iNumRows - iGetNumDof()) + (iNumDofsImpact * k) + (2 * i));
            if (it->second->SetOffsets(i)) {
                //printf(": YES\n");
                it->second->SetFirstReactionIndex(iFirstReactionIndex + (iNumDofsImpact * k) + (2 * i));
                try {
                    it->second->AssRes(WorkVec, dCoef, XCurr, XPrimeCurr);
                    for (int r = 1; r <= 6; r++) {
                        MyVec.IncCoef(iRow1 + r, WorkVec(r));
                        MyVec.IncCoef(iRow2 + r, WorkVec(6 + r));
                    }
                    for (int r = 1; r <= 2; r++) {
                        MyVec.IncCoef(iRowDof + r, WorkVec(2 * iNumRowsNode + r));
                    }
                } catch (Elem::ChangedEquationStructure) {
                    ChangeJac = true;
                    //printf("Changed Jac\n");
                }
                collisions.insert(object_pair);
                //printf("Collision!!\n");
            } else {
                //printf(": NO\n");
                for (int r = 1; r <= 2; r++) {
                    MyVec.IncCoef(iRowDof + r, - XCurr(iFirstReactionIndex + (iNumDofsImpact * k) + (2 * i) + r));
                }
            }
        }
        for (int i = it->second->ContactsSize(); i < iNumDofsImpact / 2; i++) {
            int iRowDof((iNumRows - iGetNumDof()) + (iNumDofsImpact * k) + (2 * i));
            for (int r = 1; r <= 2; r++) {
                MyVec.IncCoef(iRowDof + r, - XCurr(iFirstReactionIndex + (iNumDofsImpact * k) + (2 * i) + r));
            }
        }
    }
    if (ChangeJac) {
        throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
    }
    WorkVec.ResizeReset(iNumRows);
    for (int r = 1; r <= iNumRows; r++) {
        WorkVec.PutCoef(r, MyVec(r));
    }
    for (int i = 0; i < nodes_vector.size(); i++) {
        const integer iFirstMomentumIndex(nodes_vector[i]->iGetFirstMomentumIndex());
        for (int r = 1; r <= 6; r++) {
            WorkVec.PutRowIndex((6 * i) + r, iFirstMomentumIndex + r);
        }
    }
    for (int i = 0; i < objectpair_index_map.size(); i++) {
        for (int increment = (iNumDofsImpact * i) + 1; increment <= (iNumDofsImpact * i) + iNumDofsImpact; increment++) {
            WorkVec.PutRowIndex(iFirstRowDof + increment, iFirstReactionIndex + increment);
        }
    }
    for (std::map<ObjectPair, unsigned int>::const_iterator it = objectpair_index_map.begin();
        it != objectpair_index_map.end(); it++) {
        std::set<ObjectPair>::const_iterator collision_it = collisions.find(it->first);
        if (collision_it == collisions.end()) {
            //printf("CollisionWorld::AssRes() Miss!\n");
            for (int increment = (iNumDofsImpact * (it->second)) + 1; increment <= (iNumDofsImpact * (it->second)) + iNumDofsImpact; increment++) {
                WorkVec.PutCoef(iFirstRowDof + increment, - XCurr(iFirstReactionIndex + increment));
            }
        }
    }
    return WorkVec;
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

int
CollisionWorld::iGetNumConnectedNodes(void) const
{
    printf("Entering CollisionWorld::iGetNumConnectedNodes()\n");
    // return number of connected nodes
    return nodes_vector.size(); //0;
}

void
CollisionWorld::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
    printf("Entering CollisionWorld::GetConnectedNodes()\n");
    // resize to number of connected nodes
    // fill array with connected nodes
    //connectedNodes.resize(0);
    connectedNodes.resize(iGetNumConnectedNodes());
    for (int i = 0; i < nodes_vector.size(); i++) {
        connectedNodes[i] = nodes_vector[i];
    }
}

void
CollisionWorld::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
    //printf("Entering CollisionWorld::AfterPredict()\n");
    const integer iFirstReactionIndex = iGetFirstIndex();
    for (int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
        X.PutCoef(iFirstReactionIndex + iCnt, 0.0);
        XP.PutCoef(iFirstReactionIndex + iCnt, 0.0);
    }
    AfterConvergence(X, XP);
}

void
CollisionWorld::SetValue(DataManager *pDM,
    VectorHandler& X, VectorHandler& XP,
    SimulationEntity::Hints *ph)
{
    printf("Entering CollisionWorld::SetValue()\n");
    AfterConvergence(X, XP);
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
    // help
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
            "       <shape> [,margin, <margin>]\n"
            "\n"
            "   <shape> ::= {\n"
            "       btBoxShape, <x_half_extent>, <y_half_extent>, <z_half_extent>\n"
            "       | btCapsuleShape, <radius>, <height>\n"
            "       | btConeShape, <radius>, <height>\n"
            "       | btSphereShape, <radius>"
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
    } else {
        silent_cerr("collision object(" << GetLabel() << "): a valid btShape type is expected at line " << HP.GetLineData() << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
    if (HP.IsKeyWord("margin")) {
        btScalar margin(HP.GetReal());
        shape->setMargin(margin);
    }
    printf("getMargin() %f\n", shape->getMargin());
    ob->setCollisionShape(shape);
    universe.label_object_map[uLabel] = ob;
    universe.object_node_map[ob] = pNode;
    
    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

CollisionObject::~CollisionObject(void)
{

    delete ob;
    delete shape;
    NO_OP;
}

void
CollisionObject::Output(OutputHandler& OH) const
{
    if (fToBeOutput()) {
        // std::ostream& out = OH.Loadable();

        // write output as appropriate
    }
}

void
CollisionObject::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    // set *piNumRows to 6 * number_of_nodes
    // set *piNumCols to 1 for residual only
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
    // return number of private data
    return 0;
}

unsigned int
CollisionObject::iGetPrivDataIdx(const char *s) const
{
    // parse string and compute index of requested private data

    // shouldn't get here until private data are defined
    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

doublereal
CollisionObject::dGetPrivData(unsigned int i) const
{
    ASSERT(i > 1 && i <= iGetNumPrivData());

    // compute requested private data

    // shouldn't get here until private data are defined
    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

int
CollisionObject::iGetNumConnectedNodes(void) const
{
    // return number of connected nodes
    return 0;
}

void
CollisionObject::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
    // resize to number of connected nodes
    // fill array with connected nodes
    connectedNodes.resize(0);
    // connectedNodes[0] = m_pNode;
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

