#include <raytracing/integrator.h>
#include <scene/geometry/cube.h>
#include<vector>
#include<random>
#define Np 0x1000//4096

float N = 5;
float EPSILON=0.00001;
int LigthDepth =2;
//--------------------------------------------------
//------------------volumetric rendering------------
//--------------------------------------------------
static const int PSample=1024;
static const int BM=PSample-1;
static const int voxel_size = 100;
//------------------------------------------------
//--------------perlin noise---------------------
//-------------------------------------------------
#define lerp(a, b,t) (a + t * (b - a))
#define s_curve(t) (t * t * (3.0f - 2.0f * t))
#define setup(i,b0,b1,r0,r1)\
    t = i + N;\
    b0 = ((int)t) & BM;\
    b1 = (b0+1) & BM;\
    r0 = t - (int)t;\
    r1 = r0 - 1.0f;

float noise(glm::vec3 mnoise,bool initial) {

    if (initial) {
            qsrand(1);
        }

    QList<int>p;
    QList<glm::vec3>g3;
    for (int i = 0 ; i < PSample ; i++) {
        p.push_back(i);
        float r1=(float)(qrand() % (2*PSample)-PSample)/(float)PSample;// - PSample) / PSample;
        float r2=(float)(qrand() % (2*PSample)-PSample)/(float)PSample;
        float r3=(float)(qrand() % (2*PSample)-PSample)/(float)PSample;
        glm::vec3 randvec3=glm::vec3(r1,r2,r3);
        randvec3=glm::normalize(randvec3);
        g3.push_back(randvec3);
        }//(-1,1)random number with 1024 precise
    for(int i=PSample-1; i < 0; i--)
    {
        int k=p[i];
        int t=qrand()% (PSample);
        p[i]=p[t];
        p[t]=k;
    }
    for(int i=0;i<PSample+2;i++)
    {
        p.push_back(p[i]);
        g3.push_back(g3[i]);
    }
    //glm::vec3 b0,b1,r0,r1;

    int bx0, bx1, by0, by1, bz0, bz1, b00, b10, b01, b11;
    float rx0, rx1, ry0, ry1, rz0, rz1, sy, sz, a, b, c, d, t, u, v;
    int i, j;
    glm::vec3 q;

    setup(mnoise.x, bx0,bx1, rx0,rx1);//
    setup(mnoise.y, by0,by1, ry0,ry1);
    setup(mnoise.z, bz0,bz1, rz0,rz1);

    i = p[bx0];
    j = p[bx1];

    b00 = p[i + by0];
    b10 = p[j + by0];
    b01 = p[i + by1];
    b11 = p[j + by1];

    t  = s_curve(rx0);
    sy = s_curve(ry0);
    sz = s_curve(rz0);

    #define at3(rx,ry,rz) ( rx * q[0] + ry * q[1] + rz * q[2] )

    q = g3[b00 + bz0] ; u=at3(rx0,ry0,rz0);
    q = g3[b10 + bz0] ; v =at3(rx1,ry0,rz0);
    a = lerp(u, v, t);

    q = g3[b01 + bz0] ; u = at3(rx0,ry1,rz0);
    q = g3[b11 + bz0] ; v = at3(rx1,ry1,rz0);
    b = lerp(u, v, t);

    c = lerp(a, b, sy);

    q = g3[b00 + bz1] ; u = at3(rx0,ry0,rz1);
    q = g3[b10 + bz1] ; v = at3(rx1,ry0,rz1);
    a = lerp(u, v, t);

    q = g3[b01 + bz1] ; u = at3(rx0,ry1,rz1);
    q = g3[b11 + bz1] ; v = at3(rx1,ry1,rz1);
    b = lerp(u, v, t);
    d = lerp(a, b,sy);

    return lerp(c, d, sz);
}
float perlin_noise_3D(glm::vec3 p,bool initial) {
    int Octaves=6;
    float  Frequency=0.125;
    float result = 0.0f;
    int seed=1;
    float Amplitude=2;

    p *= Frequency;
    for(int i=0; i<Octaves; i++ ) {
        result += noise(p,initial)*Amplitude;
        p *= 2.0f;
        Amplitude*=0.5f;
    }

    return result;
}
glm::vec3 World_To_Zero_One(glm::vec3 world_p)
{
    //return glm::vec3(world_p.x/10.0+1.0,world_p.y/10.0+1.0,world_p.z/10.0+1.0);[-5,5]
return glm::vec3(world_p.x,world_p.y,world_p.z);
}

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
bool in_voxel(glm::vec3 p)
{
    float max_abcd=1.f;
    float min_abcd=0.f;
    if(p.x<max_abcd && p.x>=min_abcd &&p.y<max_abcd && p.y>=min_abcd &&p.z<=0 && p.z>-1)
        return true;
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
    //-------------get the random light---
    int light_index;
    float ca=1.0/float(scene->lights.size());
    for(int i=0;i<scene->lights.size();i++)
    {
        if(x<(i+1)*ca){
            light_index=i;
            break;
        }
        else continue;
    }
    //-------------------------------------
    for(int j=0;j<N;j++)
    {
            int i=light_index;
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


               if(fabs(brdf_pdf-0.f)>0.01){
                //radiance*light_brdf*(float)fabs(glm::dot(isx.normal,light_ray.direction))/light_pdf;
                glm::vec3 sample_brdf_color = Ld*brdf_brdf*(float)fabs(glm::dot(isx.normal,new_ray.direction))/brdf_pdf;
                brdf_color =sample_brdf_color;

                f_brdf_color=brdf_brdf*(float)fabs(glm::dot(isx.normal,new_ray.direction));
                fcos=std::max(std::max(brdf_brdf.x,brdf_brdf.y),brdf_brdf.z)*(float)fabs(glm::dot(isx.normal,new_ray.direction));
                my_tracing_pdf=brdf_pdf;
                 }
                float w_light=pow(light_pdf,2.0)/(pow(light_pdf,2.0)+pow(brdf_pdf,2.0));
                float w_brdf=pow(brdf_pdf,2.0)/(pow(light_pdf,2.0)+pow(brdf_pdf,2.0));

                glm::vec3 com_color =(light_color*w_light+brdf_color*w_brdf);
                color+=com_color;

        }

     return color/(float)N;
}
glm::vec3 Integrator::TraceRay(Ray r, unsigned int depth)
{
    if(depth>=max_depth)return glm::vec3(0,0,0);
    glm::vec3 Emitted_Light(0,0,0);
    glm::vec3 Direct_light(0,0,0);
    glm::vec3 Reflect_light(0,0,0);


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

//--------------------------------
//-------------bidirectional-----------
//-------------------------------

Ray Integrator::BounceRay(Ray r,glm::vec3 &color)
{
    Ray next_ray;
    Intersection isx=  intersection_engine->GetIntersection(r);
    glm::vec3 surface_color=isx.object_hit->material->base_color;

    next_ray.origin=isx.point;
    float brdf_pdf;
    glm::vec3 brdf_brdf=isx.object_hit->material->SampleAndEvaluateScatteredEnergy(isx,-r.direction,next_ray.direction,brdf_pdf)*surface_color;
    color=brdf_brdf;
    next_ray.origin+=next_ray.direction*(float)0.001;
    return next_ray;

}
bool Integrator::Recurse_Lightray(Ray &r,int depth)
{
    if(depth>LigthDepth)
    {
        return true;
    }
    glm::vec3 color;
    Ray rNew=BounceRay(r,color);

    LightNode ln;
    ln.isx=intersection_engine->GetIntersection(rNew);
    ln.color=color;
    lightnode.push_back(ln);
    r=rNew;
    return false;

}
void Integrator::ConstructLightPath(Ray r)
{
    std::mt19937 mg(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<float> u01(0.f,1.f);
    int LightSize=scene->lights.size();
    glm::vec3 color;
    Ray current;
    for (int i=0;i<LightSize;i++)
    {
        LightNode ln;
        float random1=u01(mg);
        float random2=u01(mg);
        Intersection isx= intersection_engine->GetIntersection(r);

        ln.color=scene->lights[i]->material->base_color;
        ln.isx=scene->lights[i]->GetSurfaceSample(random1, random2,isx.normal);

        current.origin=ln.isx.point;
        current.direction=glm::normalize(ln.isx.point-scene->lights[i]->transform.position());
        ln.lightId=i;
        lightnode.push_back(ln);//starter node:depth=0
    }

    //depth=1
    for(int i=0;i<LigthDepth;i++)
    {

        if(Recurse_Lightray( current,i))
        {
            break;
        }
    }
}

//--------------------------------
//-------------volumetric-----------
//-------------------------------

float Integrator::SetVoxel_Light(int num,glm::vec3 voxelpos)
{
    if(scene->voxels[num].light>0)return scene->voxels[num].light;
    float step=0.005;
    glm::vec3 ppos;
    int w=100;
    int h=100;
    int l=100;
    int new_num;
    float light=1;
    float Qlight=1;
    float density=0.f;

    for(int i=0;i<scene->lights.size();i++)
    {
        glm::vec3 lpos=scene->lights[i]->transform.position();
        glm::vec3 light_dir=lpos-voxelpos;
        //ray marching
        for(int j=0;;j++)
        {
            if(j!=0)voxelpos+=step*light_dir;
            if(in_voxel(voxelpos))
            {
                ppos=voxelpos*glm::vec3(w,h,-l);
                new_num=ppos.x+ppos.y*h+ppos.z*w*h;
                density=scene->voxels[new_num].density*w;
                if(!fequal(density,0))
                {
                    light=glm::exp(step*density);
                    Qlight*=light;
                }

            }
            else break;
        }
    }
    scene->voxels[num].light=+Qlight;
    return Qlight;


}
void Integrator ::SetVoxel_Density()
{
   // int *x=new int [1000000];
    bool initial=true;
    for(int i=0;i<voxel_size*voxel_size*voxel_size;i++)
    {
        Voxel v;
        v.density=0;
        v.light=-1;
        scene->voxels.push_back(v);
    }

    for(int k=0;k<voxel_size;k++)//zlength
        for(int j=0;j<voxel_size;j++)//ylength
            for(int i=0;i< voxel_size; i++)//xlength
            {
                //---------------------------------------------
                //---------initilize the voxel of the scene------
                //-------------------------------------------------
                float ori_dens=0;
                float noise_dens=0;
                float final_dens=0;
                //-----------------------------------------------------------
                //--calculate distance between each voxel and objects--------
                //--both range from (0-1)--set the density to dens/voxel_size;
                //-----------------------------------------------------------
                for(int obj=0;obj<scene->objects.size();obj++)
                {
                    glm::vec3 voxel_local_pos=glm::vec3((float)i/float(voxel_size),(float)j/float(voxel_size),(float)k/float(voxel_size));
                    glm::vec3 obj_world_pos = scene->objects[obj]->transform.position();
                    glm::vec3 obj_local_pos = World_To_Zero_One(obj_world_pos);

                    float dis=glm::distance(voxel_local_pos,obj_local_pos);//0-1

                    glm::vec3 scl=scene->objects[obj]->transform.getScale();
                    glm::vec3 local_scl=World_To_Zero_One(scl);//0-1
                    //--!!!!!only consider sphere first//
                    if(dis<local_scl.x)
                    {
                        ori_dens+= dis/local_scl.x;
                    }
                    glm::clamp(ori_dens,0.f,1.f);

                    //-------------------------------------
                    //------------set cloud noise--------
                    if(scene->objects[obj]->voxelsize.x>0.f)
                    {
                        glm::vec3 p=voxel_local_pos-obj_local_pos;
                        float p_noise=perlin_noise_3D(p, initial);
                        initial=false;
                        if(p_noise<0){p_noise=0;}
                        noise_dens+=p_noise;
                    }
                    //clamp...??
                    glm::clamp(noise_dens,0.f,1.f);
                    //-----------------------------------------

                }
                final_dens=noise_dens+ori_dens;
                glm::clamp(final_dens,0.f,1.0f);
                scene->voxels[i+j*voxel_size+k*voxel_size*voxel_size].density=final_dens;
                //if(!fequal(final_dens,0.0))std::cout<<final_dens;
            }

}
//point light
bool cub_intersect(Ray r, glm::vec3 Bl, glm::vec3 Bh)//in:rdir,bl,bh//out:if intersect
    {
    float a;
    glm::vec3 R=r.direction;
    glm::vec3 Eyep=r.origin;
    float Tnear,Tfar;
    ////////////////////YOZ
    if(R.x==0&&(Eyep.x<Bl.x||Eyep.x>Bh.x))return false;
    else if(R.y==0&&(Eyep.y<Bl.y||Eyep.y>Bh.y))return false;
    else if(R.z==0&&(Eyep.z<Bl.z||Eyep.z>Bh.z))return false;
    else{
        float t1, t2;
        Tnear = -100000;
        Tfar = 100000;
    ///////////Z
        t1 = (Bl.z -Eyep.z)/R.z;
        t2 = (Bh.z -Eyep.z)/R.z;
        if(t1>t2){a = t1;t1 = t2;t2 = a;}
        if(t1>Tnear) Tnear = t1;
        if(t2<Tfar)Tfar = t2;
                ////////////x
        t1 = (Bl.x - Eyep.x)/R.x;
        t2 = (Bh.x -Eyep.x)/R.x;
        if(t1>t2){a = t1;t1 = t2;t2 = a;}
        if(t1>Tnear)Tnear = t1;
        if(t2<Tfar) Tfar = t2;
        ////////////////////Y
        t1 = (Bl.y - Eyep.y)/R.y;
        t2 = (Bh.y- Eyep.y)/R.y;
        if(t1>t2){a = t1;t1 = t2;t2 = a;}
        if(t1>Tnear) Tnear = t1;
        if(t2<Tfar) Tfar = t2;
        if(Tnear<Tfar)
            return true;
        else return false;
    }
}
glm::vec3 Integrator::Volumetric_Render(Ray r,bool if_first)//ray marching
{
    glm::vec3 blue=glm::vec3(0,0,1);
    glm::vec3 red=glm::vec3(1,0,0);
    glm::vec3 Bl=glm::vec3(0,0,0);
    glm::vec3 Bh=glm::vec3(1,1,1);

    //12-5=7;
    float tNear=0.1;
    glm::vec3 ray_begin=tNear*r.direction+r.origin;

   // float maxstep=10-tNear;
    glm::vec3 acc_color=glm::vec3(0,0,0);

    bool if_intersect=cub_intersect(r, Bl, Bh);
    if(if_intersect){

        for(float i=0;i<1000;i++)
        {
            ray_begin+=step*r.direction;
            //std::cout<<ray_begin.x<<","<<ray_begin.y<<","<<ray_begin.z<<","<<std::endl;
            if(in_voxel(ray_begin))
            {
                return red;
                std::cout<<"ss";
                //------world to voxel--------
               // glm::vec3 lenth=
                glm::vec3 local=glm::vec3(ray_begin-Bl)/glm::vec3(Bh-Bl);
               // glm::vec3 longray=ray_begin*glm::vec3(voxel_size);
                int index=local.x+local.y*voxel_size+ local.z*voxel_size*voxel_size;
                float rouxi=scene->voxels[index].density;
                if(!fequal(rouxi,0))std::cout<<"-:"<<rouxi;
                //if(!fequal(rouxi,0))std::cout<<index<<":"<<rouxi<<std::endl;
                /*if(rouxi!=0){

                    deltT=glm::exp(-i*rouxi);
                    T*=deltT;
                    acc_color+=(1-deltT)*surface_color*T*SetVoxel_Light(index,ray_begin);
                }*/
                acc_color+=rouxi;

            }
            else break;
        }
        return acc_color;

    }
    else return blue;
}
  /*  Cube c;
    Intersection isx=  intersection_engine->GetIntersection(r);
    if(isx.t<0){return black;}

    glm::vec3 surface_color=isx.object_hit->material->base_color;
    //use point light:so the light does not need to be shown



    float step=0.005;
    float tNear=scene->camera.near_clip;
    glm::vec3 ray_begin=tNear*r.direction+r.origin;
    float maxstep=scene->camera.far_clip-tNear;
    glm::vec3 acc_color;

    float deltT=0,rouxi=0,T=1;
    for(float i=0;i<maxstep;i+=step)
    {
        ray_begin+=i*r.direction;
        if(in_voxel(ray_begin))
        {
            glm::vec3 longray=ray_begin*glm::vec3(w,l,-h);
            int index=longray.x+longray.y*w+ longray.z*w*w;
            rouxi=scene->voxels[index].density*w;
            if(rouxi!=0){

                deltT=glm::exp(-i*rouxi);
                T*=deltT;
                acc_color+=(1-deltT)*surface_color*T*SetVoxel_Light(index,ray_begin);
            }

        }
        else break;
    }
    acc_color+=T*surface_color;
    return acc_color;*/



void Integrator::SetDepth(unsigned int depth)
{
    max_depth = depth;
}

