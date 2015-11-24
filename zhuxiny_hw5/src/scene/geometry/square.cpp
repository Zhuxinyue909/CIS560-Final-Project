#include <scene/geometry/square.h>

void SquarePlane::ComputeArea()
{
    //TODO
    area=1.0*this->transform.getScale().x*this->transform.getScale().y;

}
Intersection SquarePlane::GetSurfaceSample(float u1, float u2, const glm::vec3 &isx_normal)
{
    Intersection result;
    glm::vec3 local_p;
    result.point=SetRandomPoint(local_p, u1, u2);
    glm::vec4 Lnormal=glm::vec4(ComputeNormal(result.point),0);
    result.normal=glm::normalize(glm::vec3(transform.invTransT()*Lnormal));
    glm::vec2 uv = this->GetUVCoordinates(glm::vec3(local_p));
    glm::vec3 color = Material::GetImageColor(uv, this->material->texture);
    result.object_hit=this;
    result.texture_color=color;


    return result;
}
glm::vec3 SquarePlane:: SetRandomPoint(glm::vec3& local_p,float x,float y)
{
    /*std::random_device rd;
    std::mt19937 mg(rd());
    std::uniform_real_distribution<float> u01(0.f,1.f);

    float x = u01(mg)-0.5;
    float y = u01(mg)-0.5;//[-0.5,+0.5]*/
    glm::vec3 randomp=glm::vec3(x-0.5,y-0.5,0.0);

    glm::vec4 point=transform.T()*glm::vec4(randomp,1.0);
    local_p=randomp;
    return glm::vec3(point.x,point.y,point.z);


}
Intersection SquarePlane::GetIntersection(Ray r)
{
    //Transform the ray
    Ray r_loc = r.GetTransformedCopy(transform.invT());
    Intersection result;

    //Ray-plane intersection
    float t = glm::dot(glm::vec3(0,0,1), (glm::vec3(0.5f, 0.5f, 0) - r_loc.origin)) / glm::dot(glm::vec3(0,0,1), r_loc.direction);
    glm::vec4 P = glm::vec4(t * r_loc.direction + r_loc.origin, 1);
    //Check that P is within the bounds of the square
    if(t > 0 && P.x >= -0.5f && P.x <= 0.5f && P.y >= -0.5f && P.y <= 0.5f)
    {
        result.point = glm::vec3(transform.T() * P);
        result.normal = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(ComputeNormal(glm::vec3(P)), 0)));
        result.object_hit = this;
        result.t = glm::distance(result.point, r.origin);
        result.texture_color = Material::GetImageColorInterp(GetUVCoordinates(glm::vec3(P)), material->texture);
        //TODO: Store the tangent and bitangent
        result.tangent = glm::normalize(glm::cross(glm::vec3(0,1,0),result.normal));
        result.bitangent = glm::normalize(glm::cross(result.normal,result.tangent));
        return result;
    }
    return result;
}


glm::vec2 SquarePlane::GetUVCoordinates(const glm::vec3 &point)
{
    return glm::vec2(point.x + 0.5f, point.y + 0.5f);
}

glm::vec3 SquarePlane::ComputeNormal(const glm::vec3 &P)
{
        return glm::vec3(0,0,1);
}

void SquarePlane::create()
{
    GLuint cub_idx[6];
    glm::vec3 cub_vert_pos[4];
    glm::vec3 cub_vert_nor[4];
    glm::vec3 cub_vert_col[4];

    cub_vert_pos[0] = glm::vec3(-0.5f, 0.5f, 0);  cub_vert_nor[0] = glm::vec3(0, 0, 1); cub_vert_col[0] = material->base_color;
    cub_vert_pos[1] = glm::vec3(-0.5f, -0.5f, 0); cub_vert_nor[1] = glm::vec3(0, 0, 1); cub_vert_col[1] = material->base_color;
    cub_vert_pos[2] = glm::vec3(0.5f, -0.5f, 0);  cub_vert_nor[2] = glm::vec3(0, 0, 1); cub_vert_col[2] = material->base_color;
    cub_vert_pos[3] = glm::vec3(0.5f, 0.5f, 0);   cub_vert_nor[3] = glm::vec3(0, 0, 1); cub_vert_col[3] = material->base_color;

    cub_idx[0] = 0; cub_idx[1] = 1; cub_idx[2] = 2;
    cub_idx[3] = 0; cub_idx[4] = 2; cub_idx[5] = 3;

    count = 6;

    bufIdx.create();
    bufIdx.bind();
    bufIdx.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufIdx.allocate(cub_idx, 6 * sizeof(GLuint));

    bufPos.create();
    bufPos.bind();
    bufPos.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufPos.allocate(cub_vert_pos, 4 * sizeof(glm::vec3));

    bufNor.create();
    bufNor.bind();
    bufNor.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufNor.allocate(cub_vert_nor, 4 * sizeof(glm::vec3));

    bufCol.create();
    bufCol.bind();
    bufCol.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufCol.allocate(cub_vert_col, 4 * sizeof(glm::vec3));
}
