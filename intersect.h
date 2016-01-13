#include <boost/foreach.hpp>
#include <fcl/shape/geometric_shapes.h>
#include <fcl/shape/geometric_shapes_utility.h>
#include <fcl/narrowphase/narrowphase.h>
#include <fcl/broadphase/broadphase.h>
#include <fcl/collision.h>

namespace FCL
{
typedef boost::shared_ptr <const fcl::CollisionGeometry> constCollisionGeometryPtr_t;
typedef boost::shared_ptr <fcl::CollisionGeometry> CollisionGeometryPtr_t;
typedef std::pair<fcl::CollisionObject*, fcl::CollisionObject*> ObjectPair;
typedef std::vector<std::pair<fcl::Vec3f, fcl::Vec3f> > Vec3f_pairs;
typedef void (*Func)(fcl::CollisionObject* pObject1, fcl::CollisionObject* pObject2, Vec3f_pairs& Rf_pairs);

class FuncMatrix {
private:
    Func funcs[fcl::NODE_COUNT][fcl::NODE_COUNT];
public:
    FuncMatrix(void);
    Func GetFunc(ObjectPair object_pair);
};

} // FCL

