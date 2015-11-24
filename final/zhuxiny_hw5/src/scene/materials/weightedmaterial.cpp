#include <scene/materials/weightedmaterial.h>

WeightedMaterial::WeightedMaterial() : Material(){}
WeightedMaterial::WeightedMaterial(const glm::vec3 &color) : Material(color){}

//Given an intersection with some Geometry, evaluate the scattered energy at isx
//given a world-space wo and wi for all BxDFs we contain that match the input flags
glm::vec3 WeightedMaterial::EvaluateScatteredEnergy(const Intersection &isx, const glm::vec3 &woW, const glm::vec3 &wiW, BxDFType flags) const
{

    std::mt19937 mg(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<float> u01(0.f,1.f);
    float rand1=u01(mg);

    glm::mat4 tran_m=glm::mat4(glm::vec4(isx.tangent.x,isx.bitangent.x,isx.normal.x,0),
                               glm::vec4(isx.tangent.y,isx.bitangent.y,isx.normal.y,0),
                               glm::vec4(isx.tangent.z,isx.bitangent.z,isx.normal.z,0),
                               glm::vec4(0,0,0,1));
    glm::vec3 wo=glm::normalize(glm::vec3(tran_m*glm::vec4(woW,0)));
    glm::vec3 wi=glm::normalize(glm::vec3(tran_m*glm::vec4(wiW,0)));


    glm::vec3 result;

    if(this->bxdf_weights.size()!=0)
    {
        Geometry *geo=isx.object_hit;

        for(int i=0;i<bxdf_weights.size();i++)
        {
            if(rand1<=bxdf_weights[i])
            {
                //geo->material->bxdfs[i]->EvaluateScatteredEnergy(woW,wiW)*this->bxdf_weights[i*2];
                return geo->material->bxdfs[i]->EvaluateScatteredEnergy(wo,wi);
            }
            else
            {
                rand1-=bxdf_weights[i];
            }
        }
    }

    return result;
}
 //Given an intersection with some Geometry, generate a world-space wi
//then evaluate the scattered energy along the world-space wo.
glm::vec3 WeightedMaterial::SampleAndEvaluateScatteredEnergy(const Intersection &isx, const glm::vec3 &woW, glm::vec3 &wiW_ret, float &pdf_ret, BxDFType flags) const
{
    //TODO
    glm::vec3 dir;
    glm::vec3 color;
    std::mt19937 mg(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<float> u01(0.f,1.f);
    float rand1=u01(mg);
    float x=u01(mg);
    float y=u01(mg);
    glm::mat4 tran_m=glm::mat4(glm::vec4(isx.tangent.x,isx.bitangent.x,isx.normal.x,0),
                               glm::vec4(isx.tangent.y,isx.bitangent.y,isx.normal.y,0),
                               glm::vec4(isx.tangent.z,isx.bitangent.z,isx.normal.z,0),
                               glm::vec4(0,0,0,1));

    glm::vec3 wo=glm::normalize(glm::vec3(tran_m*glm::vec4(woW,0)));
    glm::vec3 wi=glm::normalize(glm::vec3(tran_m*glm::vec4(wiW_ret,0)));
    //color= bxdfs[0]->SampleAndEvaluateScatteredEnergy(wo,wi,x,y,pdf_ret);

    if(this->bxdf_weights.size()!=0)
    {
        Geometry *geo=isx.object_hit;

        for(int i=0;i<bxdf_weights.size();i++)
        {
            if(rand1<=bxdf_weights[i])
            {
                color=geo->material->bxdfs[i]->SampleAndEvaluateScatteredEnergy(wo,wi,x,y,pdf_ret);
                break;
            }
            else
            {
                rand1-=bxdf_weights[i];
            }
        }
    }
    glm::mat4 tran_m1=glm::transpose(tran_m);
    wiW_ret=glm::normalize(glm::vec3(tran_m1*glm::vec4(wi,0)));
    return color;

}
