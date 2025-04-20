#include "shapes.hpp"
Mesh unitCube(int m, int n, int o)
{
    //unit cube with m,n,o divisions across x,y,z axis
    std::vector<MeshVertex> vertices;
    std::vector<MeshFace> faces;
    std::vector<HalfEdge*> halfEdges;
    //faces except those along z axis

    //vertices
    for(int i=0;i<(o+1);i++)
    {
        vec3 position = vec3(0.0,0.0,(float)i/(float)o);
        for(int j=0;j<2*m+2*n;j++)
        {
            vertices.emplace_back(MeshVertex(nullptr,position));
            if(j<m) position.x+=1.0/(float)m;
            else if(j<m+n) position.y+=1.0/(float)n;
            else if(j<2*m+n) position.x-=1.0/(float)m;
            else position.y-=1.0/(float)n;
            // std::cout<<position.x<<" "<<position.y<<" "<<position.z<<std::endl;
        }
    }
    //faces
    for(int i=0;i<o;i++)
    {
        for(int j=0;j<2*m+2*n;j++)
        {
            int l1=i*(2*m+2*n)+j,l2=l1+1,l4=l1+(2*m+2*n),l3=l4+1;
            if(j==(2*m+2*n-1)) {l2=i*(2*m+2*n);l3=l2+(2*m+2*n);}
            HalfEdge *e1 = new HalfEdge(l1);
            HalfEdge *e2 = new HalfEdge(l2);
            HalfEdge *e3 = new HalfEdge(l3);
            HalfEdge *e4 = new HalfEdge(l4);
            e1->next=e2;
            e2->next=e3;
            e3->next=e4;
            e4->next=e1;
            vertices[l1].edge=e1;
            vertices[l2].edge=e2;
            vertices[l3].edge=e3;
            vertices[l4].edge=e4;
            halfEdges.emplace_back(e1);
            halfEdges.emplace_back(e2);
            halfEdges.emplace_back(e3);
            halfEdges.emplace_back(e4);

            MeshFace currFace = MeshFace(e1,vec3(0.0,0.0,1.0));
            if(j<m) currFace.normal=vec3(0.0,-1.0,0.0);
            else if (j<m+n) currFace.normal=vec3(1.0,0.0,0.0);
            else if (j<2*m+n) currFace.normal=vec3(0.0,1.0,0.0);
            else currFace.normal=vec3(-1.0,0.0,0.0);
            faces.emplace_back(currFace);

            e1->face=&faces[faces.size()-1];
            e2->face=&faces[faces.size()-1];
            e3->face=&faces[faces.size()-1];
            e4->face=&faces[faces.size()-1];

            if(i>0)
            {
                faces[(i-1)*(2*m+2*n)+j].edge->next->next->twin=e1;
                e1->twin=faces[(i-1)*(2*m+2*n)+j].edge->next->next;
            }
            if(j>0)
            {
                faces[i*(2*m+2*n)+j-1].edge->next->twin=e4;
                e4->twin=faces[i*(2*m+2*n)+j-1].edge->next;
            }
        }
        faces[i*(2*m+2*n)].edge->next->next->next->twin=faces[i*(2*m+2*n)+(2*m+2*n)-1].edge->next;
        faces[i*(2*m+2*n)+(2*m+2*n)-1].edge->next->twin=faces[i*(2*m+2*n)].edge->next->next->next;
    }

    //bottom face
    int sz=vertices.size();
    vec3 position = vec3(1.0/(float)m,1.0/(float)n,0.0);
    for(int i=0;i<n-1;i++)
    {
        for(int j=0;j<m-1;j++)
        {
            // std::cout<<position.x<<" "<<position.y<<" "<<position.z<<std::endl;
            vertices.emplace_back(MeshVertex(nullptr,position));
            position.x+=1.0/(float)m;
        }
        position.y+=1.0/(float)n;
        position.x=1.0/(float)m;
    }
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            int l1,l2,l3,l4;
            //l1
            if(i==0) l1=j;
            else if (j==0) l1=2*m+2*n-i;
            else l1=sz+(i-1)*(m-1)+j-1;
            //l2
            if(i==0) l2=j+1;
            else if (j==m-1) l2=m+i;
            else l2=sz+(i-1)*(m-1)+j;
            //l3
            if(i==n-1) l3=2*m+n-j-1;
            else if (j==m-1) l3=m+i+1;
            else l3=sz+i*(m-1)+j;
            //l4
            if(i==n-1) l4=2*m+n-j;
            else if (j==0) l4=2*m+2*n-i-1;
            else l4=sz+i*(m-1)+j-1;

            if(l1<0 || l2<0 || l3<0 || l4<0) std::cout<<l1<<" "<<l2<<" "<<l3<<" "<<l4<<std::endl;
            if(l1>=vertices.size() || l2>=vertices.size() || l3>=vertices.size() || l4>=vertices.size()) std::cout<<l1<<" "<<l2<<" "<<l3<<" "<<l4<<std::endl;

            HalfEdge *e1 = new HalfEdge(l1);
            HalfEdge *e2 = new HalfEdge(l2);
            HalfEdge *e3 = new HalfEdge(l3);
            HalfEdge *e4 = new HalfEdge(l4);
            e1->next=e4;
            e4->next=e3;
            e3->next=e2;
            e2->next=e1;
            MeshFace currFace(e1,vec3(0.0,0.0,-1.0));
            halfEdges.emplace_back(e1);
            halfEdges.emplace_back(e2);
            halfEdges.emplace_back(e3);
            halfEdges.emplace_back(e4);
            faces.emplace_back(currFace);
            if(l1>=sz) vertices[l1].edge=e1;
            if(l2>=sz) vertices[l2].edge=e2;
            if(l3>=sz) vertices[l3].edge=e3;
            if(l4>=sz) vertices[l4].edge=e4;

            e1->face=&faces[faces.size()-1];
            e2->face=&faces[faces.size()-1];
            e3->face=&faces[faces.size()-1];
            e4->face=&faces[faces.size()-1];

            int currFaceIndex = faces.size()-1;
            if(j>0)
            {
                faces[currFaceIndex-1].edge->next->next->twin=e1;
                e1->twin=faces[currFaceIndex-1].edge->next->next;
            }
            if (j==0)
            {
                // if(faces[2*m+2*n-i-1].edge->twin) std::cout<<"error"<<std::endl;
                faces[2*m+2*n-i-1].edge->twin=e1;
                e1->twin=faces[2*m+2*n-i-1].edge;
            }
            if (j==m-1)
            {
                faces[m+i].edge->twin=e3;
                e3->twin=faces[m+i].edge;
            }
            if(i>0)
            {
                faces[currFaceIndex-m].edge->next->twin=e2;
                e2->twin=faces[currFaceIndex-m].edge->next;
            }
            if (i==0)
            {
                faces[j].edge->twin=e2;
                e2->twin=faces[j].edge;
            }
            if (i==n-1)
            {
                faces[2*m+n-j-1].edge->twin=e4;
                e4->twin=faces[2*m+n-j-1].edge;
            }
        }
    }

    //top face
    sz=vertices.size();
    position = vec3(1.0/(float)m,1.0/(float)n,1.0);
    for(int i=0;i<n-1;i++)
    {
        for(int j=0;j<m-1;j++)
        {
            // std::cout<<position.x<<" "<<position.y<<" "<<position.z<<std::endl;
            vertices.emplace_back(MeshVertex(nullptr,position));
            position.x+=1.0/(float)m;
        }
        position.y+=1.0/(float)n;
        position.x=1.0/(float)m;
    }
    int offset=o*(2*m+2*n), face_offset=(o-1)*(2*m+2*n);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            int l1,l2,l3,l4;
            //l1
            if(i==0) l1=offset+j;
            else if (j==0) l1=offset+2*m+2*n-i;
            else l1=sz+(i-1)*(m-1)+j-1;
            //l2
            if(i==0) l2=offset+j+1;
            else if (j==m-1) l2=offset+m+i;
            else l2=sz+(i-1)*(m-1)+j;
            //l3
            if(i==n-1) l3=offset+2*m+n-j-1;
            else if (j==m-1) l3=offset+m+i+1;
            else l3=sz+i*(m-1)+j;
            //l4
            if(i==n-1) l4=offset+2*m+n-j;
            else if (j==0) l4=offset+2*m+2*n-i-1;
            else l4=sz+i*(m-1)+j-1;

            if(l1<0 || l2<0 || l3<0 || l4<0) std::cout<<l1<<" "<<l2<<" "<<l3<<" "<<l4<<std::endl;
            if(l1>=vertices.size() || l2>=vertices.size() || l3>=vertices.size() || l4>=vertices.size()) std::cout<<l1<<" "<<l2<<" "<<l3<<" "<<l4<<std::endl;

            HalfEdge *e1 = new HalfEdge(l1);
            HalfEdge *e2 = new HalfEdge(l2);
            HalfEdge *e3 = new HalfEdge(l3);
            HalfEdge *e4 = new HalfEdge(l4);
            e1->next=e2;
            e2->next=e3;
            e3->next=e4;
            e4->next=e1;
            MeshFace currFace(e1,vec3(0.0,0.0,1.0));
            halfEdges.emplace_back(e1);
            halfEdges.emplace_back(e2);
            halfEdges.emplace_back(e3);
            halfEdges.emplace_back(e4);
            faces.emplace_back(currFace);
            if(l1>=sz) vertices[l1].edge=e1;
            if(l2>=sz) vertices[l2].edge=e2;
            if(l3>=sz) vertices[l3].edge=e3;
            if(l4>=sz) vertices[l4].edge=e4;

            e1->face=&faces[faces.size()-1];
            e2->face=&faces[faces.size()-1];
            e3->face=&faces[faces.size()-1];
            e4->face=&faces[faces.size()-1];

            int currFaceIndex = faces.size()-1;
            if(j>0)
            {
                faces[currFaceIndex-1].edge->next->twin=e4;
                e4->twin=faces[currFaceIndex-1].edge->next;
            }
            if (j==0)
            {
                // if(faces[2*m+2*n-i-1].edge->twin) std::cout<<"error"<<std::endl;
                // if(face_offset+2*m+2*n-i-1>=faces.size()) std::cout<<offset+2*m+2*n-i-1<<std::endl;
                faces[face_offset+2*m+2*n-i-1].edge->next->next->twin=e4;
                e4->twin=faces[face_offset+2*m+2*n-i-1].edge->next->next;
            }
            if (j==m-1)
            {
                faces[face_offset+m+i].edge->next->next->twin=e2;
                e2->twin=faces[face_offset+m+i].edge->next->next;
            }
            if(i>0)
            {
                faces[currFaceIndex-m].edge->next->next->twin=e1;
                e1->twin=faces[currFaceIndex-m].edge->next->next;
            }
            if (i==0)
            {
                faces[face_offset+j].edge->next->next->twin=e1;
                e1->twin=faces[face_offset+j].edge->next->next;
            }
            if (i==n-1)
            {
                faces[face_offset+2*m+n-j-1].edge->next->next->twin=e3;
                e3->twin=faces[face_offset+2*m+n-j-1].edge->next->next;
            }
        }
    }

    //make the cube origin centered
    for(int i=0;i<vertices.size();i++) vertices[i].position-=vec3(0.5,0.5,0.5);

    return Mesh(halfEdges,vertices,faces);
}


Mesh unitSquare(int m, int n)
{
    std::vector<MeshVertex> vertices;
    std::vector<MeshFace> faces;
    std::vector<HalfEdge*> halfEdges;

    for(int i=0;i<(m+1)*(n+1);i++)
    {
        int x = i/(n+1), y = i%(n+1);
        float x_cord = (float)y/(float)n;
        float y_cord = (float)x/(float)m;
        vertices.emplace_back(MeshVertex(nullptr,vec3(x_cord,y_cord,0.0)));
    }

    for(int i=0;i<m*n;i++)
    {
        int x = i/n, y = i%n;

        int l1=x*(n+1)+y, l2=x*(n+1)+y+1, l3=l2+(n+1), l4=l1+(n+1);

        HalfEdge * e1=new HalfEdge(l1);
        HalfEdge * e2=new HalfEdge(l2);
        HalfEdge * e3=new HalfEdge(l3);
        HalfEdge * e4=new HalfEdge(l4);

        vertices[l1].edge=e1;
        vertices[l2].edge=e2;
        vertices[l3].edge=e3;
        vertices[l4].edge=e4;

        e1->next=e2;
        e2->next=e3;
        e3->next=e4;
        e4->next=e1;
        MeshFace currFace = MeshFace(e1,vec3(0.0,0.0,1.0));
        faces.emplace_back(currFace);
        e1->face=&faces[faces.size()-1];
        e2->face=&faces[faces.size()-1];
        e3->face=&faces[faces.size()-1];
        e4->face=&faces[faces.size()-1];
        if(y>0)
        {
            faces[i-1].edge->next->twin=e4;
            e4->twin=faces[i-1].edge->next;
        }
        if(x>0)
        {
            faces[i-n].edge->next->next->twin=e1;
            e1->twin=faces[i-n].edge->next->next;
        }
        halfEdges.emplace_back(e1);
        halfEdges.emplace_back(e2);
        halfEdges.emplace_back(e3);
        halfEdges.emplace_back(e4);
    }

    return Mesh(halfEdges,vertices,faces);
}

vec3 euclideanPos(vec3 polar_position)
{
    polar_position.y=radians(polar_position.y);
    polar_position.z=radians(polar_position.z);
    return vec3(
        polar_position.x * sin(polar_position.y) * cos(polar_position.z),
        polar_position.x * sin(polar_position.y) * sin(polar_position.z),
        polar_position.x * cos(polar_position.y)
    );
}

Mesh unitSphere(int m, int n)
{
    std::vector<MeshVertex> vertices;
    std::vector<MeshFace> faces;
    std::vector<HalfEdge*> halfEdges;

    float alpha = float(180)/float(n);
    float beta = float(360)/float(m);

    //quadrilateral faces
    for(int i=0;i<n-1;i++)
    {
        // std::cout<<"printing positions "<<std::endl;
        for(int j=0;j<m;j++)
        {
            float theta = 180-alpha*(i+1);
            float phi = beta*j;
            vec3 position = euclideanPos(vec3(1.0,theta,phi));
            // std::cout<<position.x<<" "<<position.y<<" "<<position.z<<std::endl;
            vertices.emplace_back(MeshVertex(nullptr,position));
        }
        // std::cout<<"done printing positions"<<std::endl;
    }
    for(int i=0;i<(n-2);i++)
    {
        for(int j=0;j<m;j++)
        {
            int l1,l2,l3,l4;
            if(j<m-1)
            {
                l1=i*m+j;l2=i*m+j+1;l3=l2+m;l4=l1+m;
            }
            else
            {
                l1=i*m+j;l2=i*m;l3=l2+m;l4=l1+m;
                // std::cout<<l1<<" "<<l2<<" "<<l3<<" "<<l4<<std::endl;
            }
            // std::cout<<l1<<" "<<l2<<" "<<l3<<" "<<l4<<std::endl;
            HalfEdge *e1 = new HalfEdge(l1);
            HalfEdge *e2 = new HalfEdge(l2);
            HalfEdge *e3 = new HalfEdge(l3);
            HalfEdge *e4 = new HalfEdge(l4);
            e1->next=e2;
            e2->next=e3;
            e3->next=e4;
            e4->next=e1;
            vertices[l1].edge=e1;
            vertices[l2].edge=e2;
            vertices[l3].edge=e3;
            vertices[l4].edge=e4;
            float theta = 180-alpha*(i+1)-alpha/(float)(2);
            float phi = beta*(j+0.5);
            vec3 currNormal=euclideanPos(vec3(1.0,theta,phi));
            currNormal=normalize(currNormal);
            MeshFace currFace = MeshFace(e1,currNormal);
            faces.emplace_back(currFace);
            halfEdges.emplace_back(e1);
            halfEdges.emplace_back(e2);
            halfEdges.emplace_back(e3);
            halfEdges.emplace_back(e4);

            e1->face=&faces[faces.size()-1];
            e2->face=&faces[faces.size()-1];
            e3->face=&faces[faces.size()-1];
            e4->face=&faces[faces.size()-1];

            if(j>0)
            {
                faces[i*m+j-1].edge->next->twin=e4;
                e4->twin=faces[i*m+j-1].edge->next;
            }
            if(i>0)
            {
                // std::cout<<"here"<<std::endl;
                faces[(i-1)*(m)+j].edge->next->next->twin=e1;
                e1->twin=faces[(i-1)*m+j].edge->next->next;
            }
        }
        faces[i*m].edge->next->next->next->twin=faces[i*m+m-1].edge->next;
        faces[i*m+m-1].edge->next->twin=faces[i*m].edge->next->next->next;
    }

    //pole
    vec3 position = vec3(0.0,0.0,-1.0);
    vertices.emplace_back(MeshVertex(nullptr,position));
    position = vec3(0.0,0.0,1.0);
    vertices.emplace_back(MeshVertex(nullptr,position));

    int sz=faces.size();

    for(int i=0;i<m;i++)
    {
        //bottom pole
        int l1,l2,l3;
        if(i<m-1)
        {
            l1=i;l2=i+1;l3=vertices.size()-2;
        }
        else
        {
            l1=i;l2=0;l3=vertices.size()-2;
        }
        HalfEdge *e1 = new HalfEdge(l1);
        HalfEdge *e2 = new HalfEdge(l2);
        HalfEdge *e3 = new HalfEdge(l3);
        e2->next=e1;
        e1->next=e3;
        e3->next=e2;
        vertices[l3].edge=e3;
        float theta = 180 - alpha/(float)(2);
        float phi = beta*(i+0.5);
        vec3 currNormal=euclideanPos(vec3(1.0,theta,phi));
        currNormal=normalize(currNormal);
        MeshFace currFace = MeshFace(e1,currNormal);
        faces.emplace_back(currFace);
        halfEdges.emplace_back(e1);
        halfEdges.emplace_back(e2);
        halfEdges.emplace_back(e3);

        e1->face=&faces[faces.size()-1];
        e2->face=&faces[faces.size()-1];
        e3->face=&faces[faces.size()-1];

        if(i>0)
        {
            faces[sz+i-1].edge->next->twin=e1;
            e3->twin=faces[sz+i-1].edge->next;
        }
        if(i==m-1)
        {
            faces[sz].edge->twin=e3;
            e3->twin=faces[sz].edge;
        }
        if(n>2)
        {
            faces[i].edge->twin=e2;
            e2->twin=faces[i].edge;
        }
    }
    sz=faces.size();
    for(int i=0;i<m;i++)
    {
        //top pole
        int l1,l2,l3;
        if(i<m-1)
        {
            l1=(n-2)*m+i;l2=l1+1;l3=vertices.size()-1;
        }
        else
        {
            l1=(n-2)*m+i;l2=(n-2)*m;l3=vertices.size()-1;
        }
        HalfEdge *e1 = new HalfEdge(l1);
        HalfEdge *e2 = new HalfEdge(l2);
        HalfEdge *e3 = new HalfEdge(l3);
        e1->next=e2;
        e2->next=e3;
        e3->next=e1;
        vertices[l3].edge=e3;
        float theta = alpha/(float)(2);
        float phi = beta*(i+0.5);
        vec3 currNormal=euclideanPos(vec3(1.0,theta,phi));
        currNormal=normalize(currNormal);
        MeshFace currFace = MeshFace(e1,currNormal);
        faces.emplace_back(currFace);
        halfEdges.emplace_back(e1);
        halfEdges.emplace_back(e2);
        halfEdges.emplace_back(e3);

        e1->face=&faces[faces.size()-1];
        e2->face=&faces[faces.size()-1];
        e3->face=&faces[faces.size()-1];

        if(i>0)
        {
            faces[sz+i-1].edge->next->twin=e3;
            e3->twin=faces[sz+i-1].edge->next;
        }
        if(i==m-1)
        {
            faces[sz].edge->next->next->twin=e2;
            e2->twin=faces[sz].edge->next->next;
        }
        if(n>2)
        {
            faces[(n-3)*m+i].edge->next->next->twin=e1;
            e1->twin=faces[(n-3)*m+i].edge->next->next;
        }
        else
        {
            if(n<2)
            {
                std::cout<<"Need at least two stacks"<<std::endl;
            }
            else
            {
                faces[i].edge->next->next->twin=e1;
                e1->twin=faces[i].edge->next->next;
            }
        }
    }

    return Mesh(halfEdges,vertices,faces);
}
