#include <raytracing/integrator.h>
#include<vector>
float N = 1;
float EPSILON=0.00001;
Integrator::Integrator():
    max_depth(5)
{
    scene = NULL;
    intersection_engine = NULL;
}

glm::vec3 ComponentMult(const glm::vec3 &a, const glm::vec3 &b)
{
    return glm::vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}


bool fequal(float x,float y)
{
    if(fabs(x-y)<=EPSILON)return true;
    else return false;
}
bool checkshadow(Scene *s, Intersection isx)
{
    for(int i=0;i<s->lights.size();i++)
    {
        if(isx.object_hit==s->lights[i]){return false;}
    }
    return true;
}
glm::vec3 Integrator::direct_light(Ray r,float x,float y,Ray &new_ray,glm::vec3 &f_brdf_color,float &my_tracing_pdf,float &fcos)
{
    glm::vec3 black=glm::vec3(0.0,0.0,0.0);
    Intersection light_isx;



    Intersection isx=  intersection_engine->GetIntersection(r);
    if(isx.t<0)
    {
        return black;
    }

    glm::vec3 surface_color=isx.object_hit->material->base_color;
    if(isx.object_hit->material->is_light_source)
    {
        return surface_color;
    }

    glm::vec3 color;
    glm::vec3 light_color;
    glm::vec3 brdf_color;
    for(int j=0;j<N;j++){
        for(int i=0;i<scene->lights.size();i++)
        {
            light_isx= scene->lights[i]->GetSurfaceSample(x,y,isx.normal);
            Ray light_ray;
            light_ray.origin=isx.point;
            light_ray.direction=glm::normalize(light_isx.point-isx.point);
            light_ray.origin+= light_ray.direction*(float)0.001;
            Intersection shadow=intersection_engine->GetIntersection(light_ray);
            glm::vec3 sample_light_color;
            float light_pdf;
            if(shadow.object_hit==scene->lights[i])
            {
                glm::vec3 radiance=scene->lights[i]->material->EvaluateScatteredEnergy(light_isx,r.direction, -light_ray.direction);

                glm::vec3 light_brdf=isx.object_hit->material->EvaluateScatteredEnergy(isx, -r.direction,light_ray.direction)*surface_color;
                light_pdf=scene->lights[i]->RayPDF(light_isx, light_ray);
                sample_light_color =radiance*light_brdf*(float)fabs(glm::dot(isx.normal,light_ray.direction))/light_pdf;
                light_color=sample_light_color;
             }
                float brdf_pdf;
                new_ray.origin=isx.point;


                glm::vec3 brdf_brdf=isx.object_hit->material->SampleAndEvaluateScatteredEnergy(isx,-r.direction,new_ray.direction,brdf_pdf)*surface_color;
                new_ray.origin+=new_ray.direction*(float)0.001;

                glm::vec3 Ld;
                Intersection isx_new=  intersection_engine->GetIntersection(new_ray);
                if(isx_new.object_hit==NULL)
                {
                    Ld=glm::vec3(0);
                }
                else if(isx_new.object_hit->material->is_light_source)
                {
                    Ld=scene->lights[i]->material->EvaluateScatteredEnergy(isx,r.direction, new_ray.direction);
                }
                else Ld=glm::vec3(0);


               // if(fabs(brdf_pdf-0.f)>0.01){
                //radiance*light_brdf*(float)fabs(glm::dot(isx.normal,light_ray.direction))/light_pdf;
               glm::vec3 sample_brdf_color = Ld*brdf_brdf*(float)fabs(glm::dot(isx.normal,new_ray.direction))/brdf_pdf;
               brdf_color =sample_brdf_color;

               f_brdf_color=brdf_brdf*(float)fabs(glm::dot(isx.normal,new_ray.direction));
               fcos=std::max(std::max(brdf_brdf.x,brdf_brdf.y),brdf_brdf.z)*(float)fabs(glm::dot(isx.normal,new_ray.direction));
               my_tracing_pdf=brdf_pdf;
               //  }
                float w_light=pow(light_pdf,2.0)/(pow(light_pdf,2.0)+pow(brdf_pdf,2.0));
                float w_brdf=pow(brdf_pdf,2.0)/(pow(light_pdf,2.0)+pow(brdf_pdf,2.0));

                glm::vec3 com_color =(light_color*w_light+brdf_color*w_brdf);
                color+=com_color;

        }
    }
     return color/((float)scene->lights.size()*N);
}
glm::vec3 Integrator::TraceRay(Ray r, unsigned int depth)
{
    std::mt19937 mg(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<float> u01(0.f,1.f);
    glm::vec3 temp_color;
    float throughput=1.0;
    Ray Li;
  /* for(int i=0;i<N;i++){
        float x=u01(mg);
        float y=u01(mg);
        temp_color+= direct_light(r,x,y,Li);
    }
    temp_color/=N;
*/
    std::vector<glm::vec3> f_brdf_cos;
    std::vector<glm::vec3> Ld_color;
    Ray in=r;
    Ray out;
    float rass=1.0;
    glm::vec3 result=glm::vec3(0);
    while(depth<max_depth)
    {
        if(depth>0){
            throughput*=rass;
        }
        if(depth>2){

            float x0=u01(mg);
            if(x0>throughput)break;
        }
        glm::vec3 d_color;
        float tracing_pdf;
        glm::vec3 ccolor;
        float x=u01(mg);
        float y=u01(mg);
        //direct_light(Ray r,float x,float y,Ray &new_ray,glm::vec3 &f_brdf_color,float &my_tracing_pdf,float &fcos)
        d_color=direct_light(in,x,y,out,ccolor,tracing_pdf,rass);
        Ld_color.push_back(d_color);
        f_brdf_cos.push_back(ccolor);

        if(tracing_pdf<0)break;
        in=out;
        depth++;
    }
    for(int i=f_brdf_cos.size()-1;i>=0;i--)
    {
        result=result*f_brdf_cos[i]+Ld_color[i];
    }
return result;
}



/*
glm::vec3 Integrator::TraceRay(Ray r, unsigned int depth)
{

    glm::vec3 black=glm::vec3(0.0,0.0,0.0);
    glm::vec3 toatl_light;
    Intersection light_isx;


    Ray Li;
    std::mt19937 mg(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<float> u01(0.f,1.f);


    for(int s=0;s<N;s++){
        glm::vec3 com_color=glm::vec3(0);
         glm::vec3 single_light_color=glm::vec3(0);
        for(int i=0;i<scene->lights.size();i++)
        {

            Li=r;
            depth=4;

            glm::vec3 acc_color=glm::vec3(1);
            while(depth<max_depth)
            {
                Intersection isx=  intersection_engine->GetIntersection(Li);
                if(isx.t<0)
                {
                    acc_color*=black;
                    break;
                }

                glm::vec3 surface_color=isx.object_hit->material->base_color;
                if(isx.object_hit->material->is_light_source)
                {
                    acc_color*=surface_color;
                    break;
                }

                glm::mat4 tran_m=glm::mat4(glm::vec4(isx.tangent,0), glm::vec4(isx.bitangent,0),glm::vec4(isx.normal,0),
                                           glm::vec4(0,0,0,1));
                glm::vec3 local_r_dir=glm::normalize(glm::vec3(tran_m*glm::vec4(Li.direction,0)));

                //direct lighting
                float x=u01(mg);
                float y=u01(mg);

                light_isx= scene->lights[i]->GetSurfaceSample(x,y,isx.normal);
                Ray light_ray;
                light_ray.origin=isx.point;
                light_ray.direction=glm::normalize(light_isx.point-isx.point);
                light_ray.origin+= light_ray.direction*(float)0.001;
                Intersection shadow=intersection_engine->GetIntersection(light_ray);

                if(!checkshadow(scene, shadow))
               {
                glm::vec3 radiance=scene->lights[i]->material->EvaluateScatteredEnergy(light_isx,r.direction, -light_ray.direction);
                glm::vec3 local_light=glm::normalize(glm::vec3(tran_m*glm::vec4(light_ray.direction,0)));
                glm::vec3 light_brdf=isx.object_hit->material->EvaluateScatteredEnergy(isx, -local_r_dir,local_light)*surface_color;
                float light_pdf=scene->lights[i]->RayPDF(light_isx, light_ray);
                glm::vec3 sample_light_color=radiance*light_brdf*(float)fabs(glm::dot(isx.normal,light_ray.direction))/light_pdf;
                return light_brdf;


                float brdf_pdf;
                glm::vec3 brdf_brdf=isx.object_hit->material->SampleAndEvaluateScatteredEnergy(isx,local_r_dir,Li.direction,brdf_pdf)*surface_color;
                glm::vec3 sample_brdf_color=radiance*brdf_brdf*(float)fabs(glm::dot(isx.normal,Li.direction))/brdf_pdf;
                if(brdf_pdf<0.f){
                    break;
                }

                Li.origin=isx.point;
                Li.origin+=Li.direction*(float)0.001;

                float w_light=pow(N*light_pdf,2.0)/(pow(N*light_pdf,2.0)+pow(N*brdf_pdf,2.0));
                float w_brdf=pow(N*brdf_pdf,2.0)/(pow(N*light_pdf,2.0)+pow(N*brdf_pdf,2.0));

                single_light_color+=sample_light_color*w_light;
                acc_color*=w_brdf*sample_brdf_color;
                }
                else acc_color*=black;
                depth++;
            }
            com_color+=(single_light_color+acc_color);

        }
        toatl_light+=(com_color/(float)scene->lights.size());

    }
    return toatl_light/(float)N;
}

*/
/*
glm::vec3 Integrator::TraceRay(Ray r, unsigned int depth)
{

    glm::vec3 black=glm::vec3(0.0,0.0,0.0);
    glm::vec3 local_p;
    glm::vec3 acc_color=glm::vec3(1,1,1);
    glm::vec3 light_color=glm::vec3(0);
    glm::vec3 sample_light_color;
    glm::vec3 brdf_color=glm::vec3(0);

    glm::vec3 com_color;
    Intersection light_isx;

    float rassl;
    Ray Li;
    Li=r;

    float throughput=1.0;
    while(depth<max_depth)
    {
        Intersection isx=  intersection_engine->GetIntersection(Li);
        if(isx.t<0) {
            return black;
        }

        glm::vec3 surface_color=isx.object_hit->material->base_color;
        if(isx.object_hit->material->is_light_source) {
            glm::vec3 radience=isx.object_hit->material->EvaluateScatteredEnergy(light_isx,r.direction, -wj.direction);
            return acc_color*surface_color*;

       }
        std::mt19937 mg(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<float> u01(0.f,1.f);

        glm::mat4 tran_m=glm::mat4(glm::vec4(isx.tangent,0), glm::vec4(isx.bitangent,0),glm::vec4(isx.normal,0),
                                   glm::vec4(0,0,0,1));

        glm::vec3 local_r_dir=glm::normalize(glm::vec3(tran_m*glm::vec4(Li.direction,0)));


        light_color=glm::vec3(0);

        if(depth>0){
            throughput*=rassl;
        }
        if(depth>1){

            float x0=u01(mg);
           // if(x0>rassl)return acc_color;

        }

       com_color=glm::vec3(0);
        for(int j=0;j<N;j++)
        {
            float x=u01(mg);
            float y=u01(mg);
            com_color=glm::vec3(0);

            for(int i=0;i<scene->lights.size();i++)
            {
                light_isx= scene->lights[i]->GetSurfaceSample(x,y,isx.normal);
#define squre_light
#ifdef squre_light
                light_isx=scene->lights[i]->material->SampleLight(x,y,scene->lights[i]);
#endif
                Ray wj;//new_light ray
                wj.origin=isx.point;
                wj.direction=glm::normalize(light_isx.point-isx.point);
                wj.origin+=wj.direction*(float)0.001;

                glm::vec3 local_wj_dir=glm::normalize(glm::vec3(tran_m*glm::vec4(wj.direction,0)));
                Intersection shadow=intersection_engine->GetIntersection(wj);
#ifdef squre_light
                light_isx.normal=shadow.normal;
#endif
                if(shadow.object_hit== scene->lights[i])
                {

                    glm::vec3 light_brdf=isx.object_hit->material->EvaluateScatteredEnergy(isx, local_r_dir,-local_wj_dir)*surface_color;
                    glm::vec3 Ld_light=scene->lights[i]->material->EvaluateScatteredEnergy(light_isx,r.direction, -wj.direction);
                    float light_pdf=scene->lights[i]->RayPDF(light_isx, wj);
                    sample_light_color=Ld_light*light_brdf*(float)fabs(glm::dot(isx.normal,wj.direction))/light_pdf;


                    wjj.origin=isx.point;
                    wjj.origin+=wjj.direction*(float)0.001;
                    float brdf_pdf;
                    glm::vec3 brdf_brdf=isx.object_hit->material->SampleAndEvaluateScatteredEnergy(isx,local_r_dir,wjj.direction,brdf_pdf)*surface_color;
                    glm::vec3 Ld_brdf=scene->lights[i]->material->EvaluateScatteredEnergy(light_isx,r.direction, -wjj.direction);
                    brdf_color=Ld_light*brdf_brdf*(float)fabs(glm::dot(isx.normal,wjj.direction))/brdf_pdf;


                    float w_light=pow(N*light_pdf,2.0)/(pow(N*light_pdf,2.0)+pow(N*brdf_pdf,2.0));
                    float w_brdf=pow(N*brdf_pdf,2.0)/(pow(N*light_pdf,2.0)+pow(N*brdf_pdf,2.0));


                    com_color+=w_light*light_color+w_brdf*brdf_color;
                    rassl=std::max(std::max(light_brdf.x,light_brdf.y),light_brdf.z)*(float)fabs(glm::dot(isx.normal,wj.direction));
                }

            }

            light_color=com_color/(float)scene->lights.size();
        }

        acc_color*=com_color/(float)N;
        depth++;
    }
    return acc_color;
}*/
/*
glm::vec3 Integrator::TraceRay(Ray r, unsigned int depth)
{

    glm::vec3 black=glm::vec3(0.0,0.0,0.0);
    glm::vec3 local_p;
    glm::vec3 radiance;
    glm::vec3 acc_color=glm::vec3(1,1,1);
    glm::vec3 light_color=glm::vec3(0);
    glm::vec3 sample_light_color;
    glm::vec3 brdf_color=glm::vec3(0);
    float light_pdf;
    float brdf_pdf;
    glm::vec3 com_color;
    Intersection light_isx;
    glm::vec3 light_brdf;
    float rassl;
    Ray wjj;
    wjj=r;
    depth=4;
    float throughput=1.0;
    while(depth<max_depth)
    {
        Intersection isx=  intersection_engine->GetIntersection(wjj);
        if(isx.t<0)
            return black;
        glm::vec3 surface_color=isx.object_hit->material->base_color;
        if(isx.object_hit->material->is_light_source)
            return surface_color;

        std::mt19937 mg(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<float> u01(0.f,1.f);

        glm::mat4 tran_m=glm::mat4(glm::vec4(isx.tangent,0), glm::vec4(isx.bitangent,0),glm::vec4(isx.normal,0),
                                   glm::vec4(0,0,0,1));

        glm::vec3 local_r_dir=glm::normalize(glm::vec3(tran_m*glm::vec4(wjj.direction,0)));


        light_color=glm::vec3(0);
        sample_light_color=glm::vec3(0);
        if(depth>0){
            throughput*=rassl;
        }
        if(depth>1){

            float x0=u01(mg);
           // if(x0>rassl)return acc_color;

        }

       com_color=glm::vec3(0);
        for(int j=0;j<N;j++)
        {
            float x=u01(mg);
            float y=u01(mg);
            sample_light_color=glm::vec3(0);

            for(int i=0;i<scene->lights.size();i++)
            {
                light_isx= scene->lights[i]->GetSurfaceSample(x,y,isx.normal);
#define squre_light
#ifdef squre_light
                light_isx=scene->lights[i]->material->SampleLight(x,y,scene->lights[i]);
#endif
                Ray wj;//new_light ray
                wj.origin=isx.point;
                wj.direction=glm::normalize(light_isx.point-isx.point);
                wj.origin+=wj.direction*(float)0.001;

                glm::vec3 local_wj_dir=glm::normalize(glm::vec3(tran_m*glm::vec4(wj.direction,0)));
                Intersection shadow=intersection_engine->GetIntersection(wj);
#ifdef squre_light
                light_isx.normal=shadow.normal;
#endif
                if(shadow.object_hit== scene->lights[i])
                {

                    light_brdf=isx.object_hit->material->EvaluateScatteredEnergy(isx, local_r_dir,-local_wj_dir)*surface_color;

                    radiance=scene->lights[i]->material->EvaluateScatteredEnergy(light_isx,r.direction, -wj.direction);
                    light_pdf=scene->lights[i]->RayPDF(light_isx, wj);
                    sample_light_color+=radiance*light_brdf*(float)fabs(glm::dot(isx.normal,wj.direction))/light_pdf;
                    rassl=std::max(std::max(light_brdf.x,light_brdf.y),light_brdf.z)*(float)fabs(glm::dot(isx.normal,wj.direction));
                }
                else sample_light_color+=black;
            }

            light_color=sample_light_color;//(float)scene->lights.size();



           // wjj.origin=isx.point;
            glm::vec3 brdf_brdf=isx.object_hit->material->SampleAndEvaluateScatteredEnergy(isx,local_r_dir,wjj.direction,brdf_pdf);
            //wjj.origin+=wjj.direction*(float)0.001;
           // Intersection brdf_isx=  intersection_engine->GetIntersection(wjj);
            float brdf_ld=scene->lights[i]->material->EvaluateScatteredEnergy(light_isx,r.direction, -wj.direction);
            brdf_color=brdf_brdf*(float)fabs(glm::dot(isx.normal,wjj.direction))/brdf_pdf;


            float w_light=pow(N*light_pdf,2.0)/(pow(N*light_pdf,2.0)+pow(N*brdf_pdf,2.0));
            float w_brdf=pow(N*brdf_pdf,2.0)/(pow(N*light_pdf,2.0)+pow(N*brdf_pdf,2.0));
            com_color+=w_light*light_color+w_brdf*brdf_color;
        }
        //light_size done:N sample
        com_color/=(float)N;
        acc_color=com_color;
        depth++;

    }
    return brdf_color;
}*/
void Integrator::SetDepth(unsigned int depth)
{
    max_depth = depth;
}
