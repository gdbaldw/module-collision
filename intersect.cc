#include "intersect.h"

namespace FCL
{

void
Intersect(const fcl::Sphere* s1, const fcl::Transform3f& tf1, const fcl::Sphere* s2, const fcl::Transform3f& tf2, Vec3f_pairs& Rf_pairs)
{
    const fcl::Vec3f normal(tf2.getTranslation() - tf1.getTranslation());
    if (normal.length() <= s1->radius + s2->radius) {
        Rf_pairs.push_back(std::make_pair(tf1.getTranslation() + normal * (s1->radius / normal.length()), tf2.getTranslation() - normal * (s2->radius / normal.length())));
    }
}

void
Intersect(const fcl::Sphere* s1, const fcl::Transform3f& tf1, const fcl::Plane* s2, const fcl::Transform3f& tf2, Vec3f_pairs& Rf_pairs)
{
    const fcl::Plane new_s2 = fcl::transform(*s2, tf2);
    const fcl::FCL_REAL signed_dist = new_s2.signedDistance(tf1.getTranslation());
    if (std::abs(signed_dist) <= s1->radius) {
       Rf_pairs.push_back(std::make_pair(tf1.getTranslation() - new_s2.n * s1->radius, tf1.getTranslation() - new_s2.n * signed_dist)); 
    }
}

template<typename T_SH1, typename T_SH2>
void
GenFunc(fcl::CollisionObject* pObject1, fcl::CollisionObject* pObject2, Vec3f_pairs& Rf_pairs)
{
    const T_SH1* s1 = static_cast<const T_SH1*>(pObject1->collisionGeometry().get());
    const T_SH2* s2 = static_cast<const T_SH2*>(pObject2->collisionGeometry().get());
    Intersect(s1, pObject1->getTransform(), s2, pObject2->getTransform(), Rf_pairs);
}

FuncMatrix::FuncMatrix(void)
{
    for(int i = 0; i < fcl::NODE_COUNT; i++) {
        for(int j = 0; j < fcl::NODE_COUNT; j++) {
            funcs[i][j] = NULL;
        }
    }
    funcs[fcl::GEOM_SPHERE][fcl::GEOM_SPHERE] = &GenFunc<fcl::Sphere, fcl::Sphere>;;
    funcs[fcl::GEOM_SPHERE][fcl::GEOM_PLANE] = &GenFunc<fcl::Sphere, fcl::Plane>;;
}

Func
FuncMatrix::GetFunc(ObjectPair object_pair) {
    return funcs[object_pair.first->getNodeType()][object_pair.second->getNodeType()];
}  

} // FCL

