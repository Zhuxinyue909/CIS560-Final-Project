#include <scene/materials/bxdfs/lambertBxDF.h>

glm::vec3 LambertBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    return this->diffuse_color/PI;

}
glm::vec3 LambertBxDF::EvaluateHemisphereScatteredEnergy(const glm::vec3 &wo, int num_samples, const glm::vec2 *samples) const
{
    //TODO
    return glm::vec3(0);
}
glm::vec3 LambertBxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{
    glm::vec3 normal=glm::vec3(0,0,1);
    float t = 2*PI*rand1;
    float u = rand2-0.5;//[-0.5,+0.5]*/


    float x=u*cos(t);
    float y=u*sin(t);
    float z=sqrt(0.25-pow(x,2.0)-pow(y,2.0));
    wi_ret=glm::normalize(glm::vec3(x,y,z));
    pdf_ret=glm::dot(wi_ret,normal)/PI;
    return EvaluateScatteredEnergy(wo, wi_ret);
}
