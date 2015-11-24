#include <scene/materials/lightmaterial.h>

glm::vec3 LightMaterial::EvaluateScatteredEnergy(const Intersection &isx, const glm::vec3 &woW, const glm::vec3 &wiW, BxDFType flags) const
{
    return glm::dot(wiW, isx.normal) > 0.0f ? (this->base_color * isx.texture_color * this->intensity) : glm::vec3(0.0f);
}
//Given an intersection with some geometry, generate a point on the geometry to which this material is applied and
//


Intersection LightMaterial::SampleLight(float x,float y,Geometry* geo) const//geo=light
{

    Intersection light_isx;
    glm::vec3 local_p;

    //else{
    glm::vec3 rpoint=geo->SetRandomPoint(local_p, x, y);
    light_isx.point=rpoint;

    glm::vec3 normal_on_light = geo->ComputeNormal(rpoint);
    glm::vec4 tm = geo->transform.invTransT()*glm::vec4(normal_on_light,0.0);
    normal_on_light=glm::vec3(tm.x,tm.y,tm.z);
    light_isx.normal=glm::normalize(normal_on_light);

    light_isx.object_hit=geo;

    glm::vec2 light_uv=geo->GetUVCoordinates(glm::vec3(local_p));
    light_isx.texture_color=Material::GetImageColorInterp(light_uv, geo->material->texture);


    return light_isx;
//    }

}
