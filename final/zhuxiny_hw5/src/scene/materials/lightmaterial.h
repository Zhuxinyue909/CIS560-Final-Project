#pragma once
#include <scene/materials/material.h>
#include <scene/geometry/sphere.h>
class LightMaterial : public Material
{
public:
    //Already implemented. Just returns the emitted light color * intensity
    virtual glm::vec3 EvaluateScatteredEnergy(const Intersection &isx, const glm::vec3 &woW, const glm::vec3 &wiW, BxDFType flags = BSDF_ALL) const;

    //Given an intersection with some geometry, generate a point on the geometry to which this material is applied and
    //
    virtual Intersection SampleLight(float x,float y,Geometry *geo) const;
};
