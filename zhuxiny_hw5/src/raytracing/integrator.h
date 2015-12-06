#pragma once
#include <la.h>
#include <raytracing/ray.h>
#include <raytracing/intersection.h>
#include <raytracing/intersectionengine.h>
#include <scene/scene.h>
#include <scene/geometry/sphere.h>
#include<random>
#include <chrono>
class Scene;
struct LightNode
{
    glm::vec3 color;
    int lightId;
    Intersection isx;
};

//------------------------------------voxel class----------

//------------------------------------------------------
//The Integrator class recursively evaluates the path a ray takes throughout a scene
//and computes the color a ray becomes as it bounces.
//It samples the materials, probability density functions, and BRDFs of the surfaces the ray hits
//to do this.

class Integrator
{
public:
    Integrator();
    Integrator(Scene *s);
    glm::vec3 Volemetric_Render(Ray r);
    glm::vec3 TraceRay(Ray r, unsigned int depth);
    glm::vec3 direct_light(Ray r,float rand1,float rand2,Ray &new_ray,glm::vec3 &f_brdf_color,float &my_tracing_pdf,float &fcos);
    void SetDepth(unsigned int depth);
    bool shadowTest(Intersection isx, Ray light);
    //Samples taken must be set in the function. It tells the how many ray samples were used to estimate the direct
    glm::vec3 EstimateDirectLighting(const Intersection &isx, unsigned int &samples_taken);
    QList <LightNode> lightnode;
    bool Recurse_Lightray(Ray &r,int depth);
    //--------------------------------------
    //--------------volumetric---------------
    //---------------------------------------

    void SetVoxel_Density();
    float SetVoxel_Light(int num,glm::vec3 voxelpos);
    glm::vec3 Volumetric_Render(Ray r,bool if_first);
    bool first;
    //-----------------------
    Ray BounceRay(Ray r,glm::vec3 &color);
    void ConstructLightPath(Ray r);
    Scene* scene;
    IntersectionEngine* intersection_engine;

protected:
    unsigned int max_depth;//Default value is 5.
};
