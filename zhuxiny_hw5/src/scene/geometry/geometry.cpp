#include <scene/geometry/geometry.h>

float Geometry::RayPDF(const Intersection &isx, const Ray &ray)
{
    //TODO
    //The isx passed in was tested ONLY against us (no other scene objects), so we test if NULL
    //rather than if != this.
    if(isx.object_hit == NULL)
    {
        return 0;
    }
    //Add more here
   //float area=isx.object_hit->area;
   float r=glm::distance2(isx.point,ray.origin);
    //r=glm::length(isx.point-ray.origin);
   float cos_theta=glm::dot(-ray.direction,isx.normal);
   //cos_theta=fabs(cos_theta);
   return r/(cos_theta*area);
}
Intersection Geometry::GetSurfaceSample(float u1, float u2, const glm::vec3 &isx_normal)
{
    return Intersection();}
