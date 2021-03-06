#pragma once

#include <scene/materials/material.h>
#include <raytracing/intersection.h>
#include <openGL/drawable.h>
#include <raytracing/ray.h>
#include <scene/transform.h>
#include <math.h>
#include <random>
class Material;
class Intersection;


class Geometry : public Drawable
{
public:
//Constructors/destructors
    Geometry() : name("GEOMETRY"), transform()
    {
        material = NULL;
    }
//Functions
    virtual ~Geometry(){}
    virtual Intersection GetIntersection(Ray r) = 0;
    virtual void SetMaterial(Material* m){material = m;}
    virtual glm::vec2 GetUVCoordinates(const glm::vec3 &point) = 0;
    virtual glm::vec3 ComputeNormal(const glm::vec3 &P) = 0;

    //Returns the solid-angle weighted probability density function given a point we're trying to illuminate and
    //a ray going towards the Geometry
    virtual float RayPDF(const Intersection &isx, const Ray &ray);
    virtual Intersection GetSurfaceSample(float u1, float u2, const glm::vec3 &isx_normal);
    //This is called by the XML Reader after it's populated the scene's list of geometry
    //Computes the surface area of the Geometry in world space
    //Remember that a Geometry's Transform's scale is applied before its rotation and translation,
    //so you'll never have a skewed shape
    virtual void ComputeArea() = 0;
    virtual glm::vec3 SetRandomPoint(glm::vec3& local_p,float x,float y) =0;

//Member variables
    QString name;//Mainly used for debugging purposes
    Transform transform;
 //   Voxel * voxels;
    glm::vec3 voxelsize;
    Material* material;
    float area;
};
