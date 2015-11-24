#include <scene/materials/bxdfs/blinnmicrofacetbxdf.h>

glm::vec3 BlinnMicrofacetBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi) const
{

    glm::vec3 N=glm::vec3(0,0,1);

    glm::vec3 wh=glm::normalize((wo+wi)/2.f);
    float a=fabs(glm::dot(wh,N));
    float M=fabs(2*glm::dot(N,wh)*glm::dot(N,wo)/glm::dot(wo,wh));
    float S=fabs(2*glm::dot(N,wh)*glm::dot(N,wi)/glm::dot(wo,wh));
    float _G=std::min(glm::min(M,S),1.f);


    //*********************

    float _D=( exponent+2.f)/(2.f*PI)*powf(a, exponent);

    float costheta_o=glm::dot(N,wo);//normalize!!!
    float costheta_i=glm::dot(N,wi);
    float brdf=_D*_G/(4.0*costheta_o*costheta_i);
//glm::clamp(brdf,0.f,1.f);
    return brdf*this->reflection_color;
}
glm::vec3 BlinnMicrofacetBxDF::EvaluateHemisphereScatteredEnergy(const glm::vec3 &wo, int num_samples, const glm::vec2 *samples) const
{
    //TODO
    return glm::vec3(0);
}
//This generates an incoming light direction wi based on rand1 and rand2 and returns the result of EvaluateScatteredEnergy based on wi.
//It "returns" wi by storing it in the supplied reference to wi. Likewise, it "returns" the value of its PDF given wi and wo in the reference to pdf.
bool samehemi(glm::vec3 x, glm::vec3 y)
{
    return x.z*y.z>0.f;
}
glm::vec3 BlinnMicrofacetBxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{


    float costheta=powf(rand1,1.f/(exponent+1.f));
    float sintheta=sqrt(std::max(0.f,1.f-costheta*costheta));
    float phi=rand2*2.f*PI;
    glm::vec3 wh;
    wh.x=cos(phi)*sintheta;
    wh.y=sin(phi)*sintheta;
    wh.z=costheta;
    if(!samehemi(wo, wh))wh=-wh;

    wi_ret=-wo+2.f*glm::dot(wo,wh)*wh;
    glm::vec3 result=EvaluateScatteredEnergy(wo,wi_ret);

    glm::vec3 wm=glm::normalize(wi_ret+wo);
    float cst=fabsf(wm.z);
    pdf_ret=(this->exponent+1.f)*pow(cst,exponent)/(2.f*PI*4.0*glm::dot(wo,wm));

   // if(glm::dot(wo,wm)<0.f)pdf_ret=0.0;
    return result;

}
