#pragma once
#include <scene/geometry/geometry.h>

//A square is aligned with the XY plane with a normal aligned with the positive Z axis. Its side length is 1, and it is centered at the origin.
//These attributes can be altered by applying a transformation matrix to the square.
class SquarePlane : public Geometry
{
    Intersection GetIntersection(Ray r);
    virtual glm::vec2 GetUVCoordinates(const glm::vec3 &point);
    virtual glm::vec3 ComputeNormal(const glm::vec3 &P);
    void create();
    Intersection GetSurfaceSample(float u1, float u2, const glm::vec3 &isx_normal);
    virtual void ComputeArea();
    virtual glm::vec3 SetRandomPoint(glm::vec3& local_p,float x,float y);
};
