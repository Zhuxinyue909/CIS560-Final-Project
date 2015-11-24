#include <scene/materials/bxdfs/specularreflectionBxDF.h>

glm::vec3 SpecularReflectionBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    //TODO

    float EPSILON=0.01;
    glm::vec3 normal=glm::vec3(0,0,1.0);
    glm::vec3 black=glm::vec3(0,0,0);
    glm::vec3 relf=glm::reflect(wi,normal);
    if(glm::dot(relf,wo)<=1.0+EPSILON&&glm::dot(relf,wo)>=1.0-EPSILON) return this->reflection_color;

    else
    return black;

}
//This generates an incoming light direction wi based on rand1 and rand2 and returns the result of EvaluateScatteredEnergy based on wi.
//It "returns" wi by storing it in the supplied reference to wi. Likewise, it "returns" the value of its PDF given wi and wo in the reference to pdf.
glm::vec3 SpecularReflectionBxDF::EvaluateHemisphereScatteredEnergy(const glm::vec3 &wo, int num_samples, const glm::vec2 *samples) const
{
    //TODO

    return glm::vec3(0);
}
glm::vec3 SpecularReflectionBxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{
    glm::vec3 normal=glm::vec3(0,0,1);

    //wi_ret=glm::normalize(glm::reflect(wo,normal));
    wi_ret=glm::normalize(glm::vec3(-wo.x,-wo.y,wo.z));
    pdf_ret=1.0;
    return this->reflection_color;//EvaluateScatteredEnergy(wo, wi_ret);
}
