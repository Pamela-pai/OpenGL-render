//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  ------------------------
//  FakeGL.cpp
//  ------------------------
//  
//  A unit for implementing OpenGL workalike calls
//  
///////////////////////////////////////////////////

#include "FakeGL.h"
#include <math.h>

//-------------------------------------------------//
//                                                 //
// CONSTRUCTOR / DESTRUCTOR                        //
//                                                 //
//-------------------------------------------------//

// constructor
FakeGL::FakeGL()
{ // constructor
} // constructor

// destructor
FakeGL::~FakeGL()
{ // destructor
} // destructor

//-------------------------------------------------//
//                                                 //
// GEOMETRIC PRIMITIVE ROUTINES                    //
//                                                 //
//-------------------------------------------------//

// starts a sequence of geometric primitives
void FakeGL::Begin(unsigned int PrimitiveType)
{ // Begin()
    // set the primitive type
    this->currentPrimitive = PrimitiveType;
} // Begin()

// ends a sequence of geometric primitives
void FakeGL::End()
{ // End()
    this->currentPrimitive = -1;
} // End()

// sets the size of a point for drawing
void FakeGL::PointSize(float size)
{ // PointSize()
    //set the point size
    this->pointSize = size;
} // PointSize()

// sets the width of a line for drawing purposes
void FakeGL::LineWidth(float width)
{ // LineWidth()
    this->lineWidth = width;
} // LineWidth()

//-------------------------------------------------//
//                                                 //
// MATRIX MANIPULATION ROUTINES                    //
//                                                 //
//-------------------------------------------------//

// set the matrix mode (i.e. which one we change)   
void FakeGL::MatrixMode(unsigned int whichMatrix)
{ // MatrixMode()
    this->currentMatMode = whichMatrix;
} // MatrixMode()

// pushes a matrix on the stack
void FakeGL::PushMatrix()
{ // PushMatrix()
    if(this->currentMatMode==FAKEGL_MODELVIEW){
        this->matStack.push(this->modelViewMat);
    }else if(this->currentMatMode == FAKEGL_PROJECTION){
        this->matStack.push(this->projectionMat);
    }
} // PushMatrix()

// pops a matrix off the stack
void FakeGL::PopMatrix()
{ // PopMatrix()
    if(this->currentMatMode==FAKEGL_MODELVIEW){
        this->modelViewMat = this->matStack.top();
        this->matStack.pop();
    }else if(this->currentMatMode == FAKEGL_PROJECTION){
        this->projectionMat = this->matStack.top();
        this->matStack.pop();
    }
} // PopMatrix()

// load the identity matrix
void FakeGL::LoadIdentity()
{ // LoadIdentity()
    if(this->currentMatMode==FAKEGL_MODELVIEW){
        this->modelViewMat.SetIdentity();
    }else if(this->currentMatMode == FAKEGL_PROJECTION){
        this->projectionMat.SetIdentity();
    }
} // LoadIdentity()

// multiply by a known matrix in column-major format
void FakeGL::MultMatrixf(const float *columnMajorCoordinates)
{ // MultMatrixf()
    Matrix4 mat4;
    for(int i=0;i<16;i++){
        int x = i%4;
        int y = i/4;
        mat4[x][y] = columnMajorCoordinates[i];
    }
    if(this->currentMatMode==FAKEGL_MODELVIEW){
        this->modelViewMat =  this->modelViewMat *mat4;
    }else if(this->currentMatMode == FAKEGL_PROJECTION){
        this->projectionMat = this->projectionMat * mat4;
    }



} // MultMatrixf()

// sets up a perspective projection matrix
void FakeGL::Frustum(float left, float right, float bottom, float top, float zNear, float zFar)
{ // Frustum()
    Matrix4 mat4;
    mat4.SetZero();
    mat4[0][0] = (2.f*zNear)/(right-left);
    mat4[1][1] = (2.f*zNear)/(top-bottom);
    mat4[0][2] = (right+left)/(right-left);
    mat4[1][2] = (top+bottom)/(top-bottom);
    mat4[2][2] = -(zFar+zNear)/(zFar-zNear);
    mat4[2][3] = -(2*zFar*zNear)/(zFar-zNear);
    mat4[3][2] = -1;

    if(this->currentMatMode == FAKEGL_MODELVIEW){
        this->modelViewMat = this->modelViewMat * mat4;
    }else if(this->currentMatMode == FAKEGL_PROJECTION){
        this->projectionMat = this->projectionMat * mat4;
    }



} // Frustum()

// sets an orthographic projection matrix
void FakeGL::Ortho(float left, float right, float bottom, float top, float zNear, float zFar)
{ // Ortho()
    Matrix4 mat4;
    mat4[0][3] = -(right+left)/(right-left);
    mat4[1][3] = -(top+bottom)/(top-bottom);
    mat4[2][3] = -(zFar+zNear)/(zFar-zNear);
    mat4[0][0] = 2.f/(right-left);
    mat4[1][1] = 2/(top-bottom);
    mat4[2][2] = -2/(zFar-zNear);
    mat4[3][3] = 1;

    if(this->currentMatMode == FAKEGL_MODELVIEW){
        this->modelViewMat =  this->modelViewMat * mat4;
    }else if(this->currentMatMode == FAKEGL_PROJECTION){
        this->projectionMat = this->projectionMat * mat4;
    }


} // Ortho()

// rotate the matrix
void FakeGL::Rotatef(float angle, float axisX, float axisY, float axisZ)
{ // Rotatef()
    Matrix4 mat4;
    float theta = angle * 3.14 / 180.f;
    mat4.SetRotation(Cartesian3(axisX, axisY, axisZ),theta);

    if(this->currentMatMode==FAKEGL_MODELVIEW){
        this->modelViewMat =  this->modelViewMat * mat4;
    }else if(this->currentMatMode == FAKEGL_PROJECTION){
        this->projectionMat = this->projectionMat * mat4;
    }

} // Rotatef()

// scale the matrix
void FakeGL::Scalef(float xScale, float yScale, float zScale)
{ // Scalef()
    Matrix4 mat4;
    mat4.SetScale(xScale, yScale, zScale);

    if(this->currentMatMode==FAKEGL_MODELVIEW){
        this->modelViewMat =  this->modelViewMat *mat4;
    }else if(this->currentMatMode == FAKEGL_PROJECTION){
        this->projectionMat =  this->projectionMat *mat4;
    }


} // Scalef()

// translate the matrix
void FakeGL::Translatef(float xTranslate, float yTranslate, float zTranslate)
{ // Translatef()
    Matrix4 mat4;
    mat4.SetTranslation(Cartesian3(xTranslate, yTranslate, zTranslate));

    if(this->currentMatMode==FAKEGL_MODELVIEW){
        this->modelViewMat =   this->modelViewMat * mat4 ;
    }else if(this->currentMatMode == FAKEGL_PROJECTION){
        this->projectionMat =  this->projectionMat * mat4;
    }
} // Translatef()

// sets the viewport
void FakeGL::Viewport(int x, int y, int width, int height)
{ // Viewport()
    viewPortMat.SetZero();
    viewPortMat[0][3] = (x+width/2.f);
    viewPortMat[1][3] = (y+height/2.f);
    viewPortMat[2][3] = 1;
    viewPortMat[0][0] = width/2.f;
    viewPortMat[1][1] = height/2.f;
    viewPortMat[2][2] = 1;
    viewPortMat[3][3] = 1;
} // Viewport()

auto FakeGL::reflect(const Cartesian3 & vec,const Cartesian3 & normal) -> Cartesian3
    {
        float dn = 2 * vec.dot(normal);
        return vec - normal * dn;
    }

//-------------------------------------------------//
//                                                 //
// VERTEX ATTRIBUTE ROUTINES                       //
//                                                 //
//-------------------------------------------------//

// sets colour with floating point
void FakeGL::Color3f(float red, float green, float blue)
{ // Color3f()
    this->colorf.red = red*255;
    this->colorf.green = green*255;
    this->colorf.blue = blue*255;


} // Color3f()

// sets material properties
void FakeGL::Materialf(unsigned int parameterName, const float parameterValue)
{ // Materialf()
    if(parameterName & FAKEGL_SHININESS){
        shinessM = parameterValue;
    }
} // Materialf()

void FakeGL::Materialfv(unsigned int parameterName, const float *parameterValues)
{ // Materialfv()
    if(parameterName & FAKEGL_EMISSION){
        this->emissionM[0] = parameterValues[0];
        this->emissionM[1] = parameterValues[1];
        this->emissionM[2] = parameterValues[2];
        this->emissionM[3] = parameterValues[3];
    }
    if(parameterName&FAKEGL_DIFFUSE){
        this->diffuseM[0] = parameterValues[0];
        this->diffuseM[1] = parameterValues[1];
        this->diffuseM[2] = parameterValues[2];
        this->diffuseM[3] = parameterValues[3];
    }
    if(parameterName&FAKEGL_SPECULAR){
        this->specularM[0] = parameterValues[0];
        this->specularM[1] = parameterValues[1];
        this->specularM[2] = parameterValues[2];
        this->specularM[3] = parameterValues[3];
    }
    if(parameterName & FAKEGL_AMBIENT){
        this->ambientM[0] = parameterValues[0];
        this->ambientM[1] = parameterValues[1];
        this->ambientM[2] = parameterValues[2];
        this->ambientM[3] = parameterValues[3];
    }

} // Materialfv()

// sets the normal vector
void FakeGL::Normal3f(float x, float y, float z)
{ // Normal3f()
    this->normal = Homogeneous4(x,y,z,0);
} // Normal3f()

// sets the texture coordinates
void FakeGL::TexCoord2f(float u, float v)
{ // TexCoord2f()

    this->textureU = u;
    this->textureV = v;

} // TexCoord2f()

// sets the vertex & launches it down the pipeline
void FakeGL::Vertex3f(float x, float y, float z)
{ // Vertex3f()
    vertexWithAttributes vertex(x,y,z);
    vertex.u = this->textureU;
    vertex.v = this->textureV;
    vertex.normal = this->normal;

    std::copy(std::begin(this->ambientM), std::end(this->ambientM), std::begin(vertex.ambientM));
    std::copy(std::begin(this->emissionM), std::end(this->emissionM), std::begin(vertex.emissionM));
    std::copy(std::begin(this->specularM), std::end(this->specularM), std::begin(vertex.specularM));
    std::copy(std::begin(this->diffuseM), std::end(this->diffuseM), std::begin(vertex.diffuseM));
    vertex.shinessM = this->shinessM;

    this->vertexQueue.push_back(vertex);
    TransformVertex();
} // Vertex3f()

//-------------------------------------------------//
//                                                 //
// STATE VARIABLE ROUTINES                         //
//                                                 //
//-------------------------------------------------//

// disables a specific flag in the library
void FakeGL::Disable(unsigned int property)
{ // Disable()
    if(property == FAKEGL_DEPTH_TEST) {this->enable_depth_test = false;}
    if(property == FAKEGL_LIGHTING) {this->enable_lighting = false;}
    if(property == FAKEGL_TEXTURE_2D) {this->enable_texture_2D = false;}
    if(property == FAKEGL_PHONG_SHADING) {this->enable_phong_shading = false;}

} // Disable()

// enables a specific flag in the library
void FakeGL::Enable(unsigned int property)
{ // Enable()
    if(property == FAKEGL_DEPTH_TEST) {this->enable_depth_test = true;
        this->depthBuffer.Resize(frameBuffer.width, frameBuffer.height);}
    if(property ==  FAKEGL_LIGHTING) {this->enable_lighting = true;}
    if(property == FAKEGL_TEXTURE_2D) {this->enable_texture_2D = true;}
    if(property == FAKEGL_PHONG_SHADING) {this->enable_phong_shading = true;}

} // Enable()

//-------------------------------------------------//
//                                                 //
// LIGHTING STATE ROUTINES                         //
//                                                 //
//-------------------------------------------------//

// sets properties for the one and only light
void FakeGL::Light(int parameterName, const float *parameterValues)
{ // Light()
    if(parameterName & FAKEGL_AMBIENT){
        ambientL[0] = parameterValues[0];
        ambientL[1] = parameterValues[1];
        ambientL[2] = parameterValues[2];
        ambientL[3] = parameterValues[3];

    }
    if(parameterName & FAKEGL_DIFFUSE){
        diffuseL[0]=parameterValues[0];
        diffuseL[1]=parameterValues[1];
        diffuseL[2]=parameterValues[2];
        diffuseL[3]=parameterValues[3];

    }
    if(parameterName & FAKEGL_SPECULAR){
        specularL[0] = parameterValues[0];
        specularL[1] = parameterValues[1];
        specularL[2] = parameterValues[2];
        specularL[3] = parameterValues[3];

    }
    if(parameterName & FAKEGL_POSITION){
        positionL[0] = parameterValues[0];
        positionL[1] = parameterValues[1];
        positionL[2] = parameterValues[2];
        positionL[3] = parameterValues[3];

    }
} // Light()

//-------------------------------------------------//
//                                                 //
// TEXTURE PROCESSING ROUTINES                     //
//                                                 //
// Note that we only allow one texture             //
// so glGenTexture & glBindTexture aren't needed   //
//                                                 //
//-------------------------------------------------//

// sets whether textures replace or modulate
void FakeGL::TexEnvMode(unsigned int textureMode)
{ // TexEnvMode()
    if(textureMode == FAKEGL_REPLACE){
        this->textureMode = FAKEGL_REPLACE;
    }
    if(textureMode == FAKEGL_MODULATE){
        this->textureMode = FAKEGL_MODULATE;
    }
} // TexEnvMode()

// sets the texture image that corresponds to a given ID
void FakeGL::TexImage2D(const RGBAImage &textureImage)
{ // TexImage2D()
    this->textureImg = textureImage;
} // TexImage2D()

//-------------------------------------------------//
//                                                 //
// FRAME BUFFER ROUTINES                           //
//                                                 //
//-------------------------------------------------//

// clears the frame buffer
void FakeGL::Clear(unsigned int mask)
{ // Clear()
    int width = this->frameBuffer.width;
    int height = this->frameBuffer.height;
    for(int i=0;i<height;i++){
        for(int j=0;j<width;j++){
            //clear frame & depth buffer
            this->frameBuffer[i][j] = backGroundColor;
            if(this->enable_depth_test) {this->depthBuffer[i][j].alpha = 255;}
        }
    }

} // Clear()

// sets the clear colour for the frame buffer
void FakeGL::ClearColor(float red, float green, float blue, float alpha)
{ // ClearColor()
    this->backGroundColor=RGBAValue(red*255,green*255,blue*255, alpha*255);
} // ClearColor()

//-------------------------------------------------//
//                                                 //
// MAJOR PROCESSING ROUTINES                       //
//                                                 //
//-------------------------------------------------//

//thread transform
// transform one vertex & shift to the raster queue
void FakeGL::TransformVertex()
{ // TransformVertex()
    auto vertex = this->vertexQueue.front();
    vertexQueue.pop_front();

    Homogeneous4 hg4(vertex.position.x, vertex.position.y, vertex.position.z);
    // model view transformation
    auto temp =  this->modelViewMat* hg4;
    // projection transformation
    Homogeneous4 screenResult = this->projectionMat* temp;
    // perspective division
    Cartesian3 ndcs = screenResult.Point();
    // viewport mapping
    Homogeneous4 ndcs4 = Homogeneous4(ndcs.x, ndcs.y, ndcs.z);
    screenResult = this->viewPortMat * ndcs4;
    screenVertexWithAttributes  screenVertex(screenResult.x, screenResult.y, screenResult.z);

    screenVertex.colour = this->colorf;
    Homogeneous4 normal4 = this->modelViewMat * vertex.normal;
    screenVertex.normal = normal4;

    if(this->enable_texture_2D){
        screenVertex.u = this->textureU;
        screenVertex.v = this->textureV;
    }
    if(this->enable_lighting){
        std::copy(std::begin(vertex.ambientM), std::end(vertex.ambientM), std::begin(screenVertex.ambientM));
        std::copy(std::begin(vertex.diffuseM), std::end(vertex.diffuseM), std::begin(screenVertex.diffuseM));
        std::copy(std::begin(vertex.specularM), std::end(vertex.specularM), std::begin(screenVertex.specularM));
        std::copy(std::begin(vertex.emissionM), std::end(vertex.emissionM), std::begin(screenVertex.emissionM));
        screenVertex.shinessM = vertex.shinessM;

    }
    // start rasterise
    this->rasterQueue.push_back(screenVertex);
    RasterisePrimitive();

} // TransformVertex()

// thread rasterise
// rasterise a single primitive if there are enough vertices on the queue
bool FakeGL::RasterisePrimitive()
{ // RasterisePrimitive()
    switch (this->currentPrimitive) {
        case 0:{
            // primitive is a point
            if(this->rasterQueue.size()>=1){
                auto vertex = rasterQueue.front();
                rasterQueue.pop_front();
                //rasterise the point
                RasterisePoint(vertex);
                return true;
            }else
                return false;
        }
        case 1:{
            // primitive is a line
            if(this->rasterQueue.size()>=2){
                auto vertex1 = rasterQueue.front();rasterQueue.pop_front();
                auto vertex2 = rasterQueue.front();rasterQueue.pop_front();
                //rasterise the line
                this->RasteriseLineSegment(vertex1, vertex2);
                return true;
            }else
                return false;
        }
        case 2:{
            // primitive is a triangle
            if(this->rasterQueue.size()>=3){
                auto vertex1 = rasterQueue.front();rasterQueue.pop_front();
                auto vertex2 = rasterQueue.front();rasterQueue.pop_front();
                auto vertex3 = rasterQueue.front();rasterQueue.pop_front();
                //rasterise the triangle
                RasteriseTriangle(vertex1, vertex2, vertex3);
                return true;
            }else
                return false;
        }
    }
    // it doesn't have enough vertices to rasterise this primitive.
    return false;
} // RasterisePrimitive()

// depth test
bool FakeGL::depthTest(int x, int y, float z){
    if(this->depthBuffer[x][y].alpha >= z/2){
        this->depthBuffer[x][y].alpha = z/2;
        return true;
    }else
        return false;
}

// rasterises a single point
void FakeGL::RasterisePoint(screenVertexWithAttributes &vertex0)
{ // RasterisePoint()
    int x = vertex0.position.x;
    int y = vertex0.position.y;
//    int z = vertex0.position.z;
    int pointSize = this->pointSize;
    if(pointSize==0) {return;}

    for(int i=x-(pointSize/2);i<=x+(pointSize/2);i++){
        for(int j=y-(pointSize/2);j<=y+(pointSize/2);j++){
            fragmentWithAttributes vertex(j,i,vertex0.colour);
            this->fragmentQueue.push_back(vertex);
        }
    }
    ProcessFragment();

} // RasterisePoint()

// rasterises a single line segment
void FakeGL::RasteriseLineSegment(screenVertexWithAttributes &vertex0, screenVertexWithAttributes &vertex1)
{ // RasteriseLineSegment()

    int x0 = vertex0.position.x;
    int y0 = vertex0.position.y;
    int x1 = vertex1.position.x;
    int y1 = vertex1.position.y;

    bool flag = false;
    if (std::abs(x0-x1)<std::abs(y0-y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        flag = true;
    }
    if (x0>x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int disX = x1-x0;
    int disY = y1-y0;
    int dis2Y = std::abs(disY)*2;
    int acc = 0;
    int y = y0;
    for (int x=x0; x<=x1; x++) {
        if (flag) {
            fragmentWithAttributes fragmentVertex(x, y, colorf);
            this->fragmentQueue.push_back(fragmentVertex);
        } else {
            fragmentWithAttributes fragmentVertex(y, x, colorf);
            this->fragmentQueue.push_back(fragmentVertex);
        }
        acc += dis2Y;
        if (acc > disX) {
            y += (y1>y0?1:-1);
            acc -= disX*2;
        }
    }
    this->ProcessFragment();
} // RasteriseLineSegment()

// rasterises a single triangle
void FakeGL::RasteriseTriangle(screenVertexWithAttributes &vertex0, screenVertexWithAttributes &vertex1, screenVertexWithAttributes &vertex2)
{ // RasteriseTriangle()



    // compute a bounding box that starts inverted to frame size
    // clipping will happen in the raster loop proper
    float minX = frameBuffer.width, maxX = 0.0;
    float minY = frameBuffer.height, maxY = 0.0;

    // test against all vertices
    if (vertex0.position.x < minX) minX = vertex0.position.x;
    if (vertex0.position.x > maxX) maxX = vertex0.position.x;
    if (vertex0.position.y < minY) minY = vertex0.position.y;
    if (vertex0.position.y > maxY) maxY = vertex0.position.y;

    if (vertex1.position.x < minX) minX = vertex1.position.x;
    if (vertex1.position.x > maxX) maxX = vertex1.position.x;
    if (vertex1.position.y < minY) minY = vertex1.position.y;
    if (vertex1.position.y > maxY) maxY = vertex1.position.y;

    if (vertex2.position.x < minX) minX = vertex2.position.x;
    if (vertex2.position.x > maxX) maxX = vertex2.position.x;
    if (vertex2.position.y < minY) minY = vertex2.position.y;
    if (vertex2.position.y > maxY) maxY = vertex2.position.y;

    // now for each side of the triangle, compute the line vectors
    Cartesian3 vector01 = vertex1.position - vertex0.position;
    Cartesian3 vector12 = vertex2.position - vertex1.position;
    Cartesian3 vector20 = vertex0.position - vertex2.position;

    // now compute the line normal vectors
    Cartesian3 normal01(-vector01.y, vector01.x, 0.0);
    Cartesian3 normal12(-vector12.y, vector12.x, 0.0);
    Cartesian3 normal20(-vector20.y, vector20.x, 0.0);

    // we don't need to normalise them, because the square roots will cancel out in the barycentric coordinates
    float lineConstant01 = normal01.dot(vertex0.position);
    float lineConstant12 = normal12.dot(vertex1.position);
    float lineConstant20 = normal20.dot(vertex2.position);

    // and compute the distance of each vertex from the opposing side
    float distance0 = normal12.dot(vertex0.position) - lineConstant12;
    float distance1 = normal20.dot(vertex1.position) - lineConstant20;
    float distance2 = normal01.dot(vertex2.position) - lineConstant01;

    // if any of these are zero, we will have a divide by zero error
    // but notice that if they are zero, the vertices are collinear in projection and the triangle is edge on
    // we can render that as a line, but the better solution is to render nothing.  In a surface, the adjacent
    // triangles will eventually take care of it
    if ((distance0 == 0) || (distance1 == 0) || (distance2 == 0))
        return;

    // create a fragment for reuse
    fragmentWithAttributes rasterFragment;

    int textureWidth = textureImg.width;
    int textureHeight = textureImg.height;

    float vertex0Light[3]={0,0,0}, vertex1Light[3]={0,0,0}, vertex2Light[3]={0,0,0};
    Cartesian3 v0Normal(vertex0.normal.x, vertex0.normal.y, vertex0.normal.z);
    Cartesian3 v1Normal(vertex1.normal.x, vertex1.normal.y, vertex1.normal.z);
    Cartesian3 v2Normal(vertex2.normal.x, vertex2.normal.y, vertex2.normal.z);


    Homogeneous4 eye{0,0,0,1};
    auto eyeDir = this->modelViewMat * eye;
    Cartesian3 eyeDir_(eyeDir.x, eyeDir.y, eyeDir.z);
//    eyeDir_ = eyeDir_.unit();



    v0Normal = v0Normal.unit();
    v1Normal = v1Normal.unit();
    v2Normal = v2Normal.unit();
    Cartesian3 lightPosition(this->positionL[0], this->positionL[1], this->positionL[2]);
    lightPosition = lightPosition.unit();

    auto v0reflectDir = reflect(lightPosition, v0Normal);
    auto v1reflectDir = reflect(lightPosition, v1Normal);
    auto v2reflectDir = reflect(lightPosition, v2Normal);



    for(int i=0;i<3;i++){


        auto v0ambient = ambientL[i]*vertex0.ambientM[i];
        auto v1ambient = ambientL[i]*vertex1.ambientM[i];
        auto v2ambient = ambientL[i]*vertex2.ambientM[i];

        auto v0diffuse = diffuseL[i]*vertex0.diffuseM[i]* std::max(v0Normal.dot(lightPosition), 0.0f);
        auto v1diffuse = diffuseL[i]*vertex1.diffuseM[i]* std::max(v1Normal.dot(lightPosition), 0.0f);
        auto v2diffuse = diffuseL[i]*vertex2.diffuseM[i]* std::max(v2Normal.dot(lightPosition), 0.0f);

        auto v0specular = specularL[i]*vertex0.specularM[i]*std::pow(std::max(v0reflectDir.dot(eyeDir_), 0.0f), this->shinessM);
        auto v1specular = specularL[i]*vertex1.specularM[i]*std::pow(std::max(v1reflectDir.dot(eyeDir_), 0.0f), this->shinessM);
        auto v2specular = specularL[i]*vertex2.specularM[i]*std::pow(std::max(v2reflectDir.dot(eyeDir_), 0.0f), this->shinessM);

        vertex0Light[i] += v0ambient + v0diffuse + v0specular+ this->emissionM[i];

        vertex1Light[i] += v1ambient + v1diffuse + v1specular+ this->emissionM[i];

        vertex2Light[i] += v2ambient + v2diffuse + v2specular + this->emissionM[i];

    }





    // loop through the pixels in the bounding box
    for (rasterFragment.row = minY; rasterFragment.row <= maxY; rasterFragment.row++)
    { // per row
        // this is here so that clipping works correctly
        if (rasterFragment.row < 0) continue;
        if (rasterFragment.row >= frameBuffer.height) continue;
        for (rasterFragment.col = minX; rasterFragment.col <= maxX; rasterFragment.col++)
        { // per pixel
            // this is also for correct clipping
            if (rasterFragment.col < 0) continue;
            if (rasterFragment.col >= frameBuffer.width) continue;

            // the pixel in cartesian format
            Cartesian3 pixel(rasterFragment.col, rasterFragment.row, 0.0);

            // right - we have a pixel inside the frame buffer AND the bounding box
            // note we *COULD* compute gamma = 1.0 - alpha - beta instead
            float alpha = (normal12.dot(pixel) - lineConstant12) / distance0;
            float beta = (normal20.dot(pixel) - lineConstant20) / distance1;
            float gamma = (normal01.dot(pixel) - lineConstant01) / distance2;

            // now perform the half-plane test
            if ((alpha < 0.0) || (beta < 0.0) || (gamma < 0.0))
                continue;

            float z = vertex0.position.z * alpha+ vertex1.position.z * beta + vertex2.position.z * gamma;


            rasterFragment.colour = alpha * vertex0.colour + beta * vertex1.colour + gamma * vertex2.colour;

            if(this->enable_depth_test){
                if(depthTest(rasterFragment.row, rasterFragment.col,z*255)){
                    // now we add it to the queue for fragment processing
                }else
                    continue;
            }



            if(this->enable_texture_2D){

                int interpU = (alpha * vertex0.u + beta * vertex1.u + gamma * vertex2.u)* textureWidth;
                int interpV = (alpha * vertex0.v + beta * vertex1.v + gamma * vertex2.v)* textureHeight;

                if(this->textureMode == FAKEGL_REPLACE){
                    rasterFragment.colour = this->textureImg[interpV][interpU];

                }else if( this->textureMode == FAKEGL_MODULATE){
                    rasterFragment.colour.modulate(this->textureImg[interpV][interpU]);
                }
            }

            if(this->enable_lighting){

                if(this->enable_phong_shading){

                    Cartesian3 interplNormal;
                    interplNormal = alpha * v0Normal + beta*v1Normal + gamma*v2Normal;

                    float interplambientM[3];
                    float interpldiffuseM[3];
                    float interplspecularMastrial[3];
                    for(int i=0;i<3;i++){
                        interplambientM[i] = alpha*vertex0.ambientM[i] + beta* vertex1.ambientM[i] + gamma*vertex2.ambientM[i];
                        interpldiffuseM[i] = alpha*vertex0.diffuseM[i] + beta * vertex1.diffuseM[i] + gamma*vertex2.diffuseM[i];
                        interplspecularMastrial[i] = alpha*vertex0.specularM[i] + beta *vertex1.specularM[i] + gamma*vertex2.specularM[i];
                    }

                    auto lightReflect = reflect(lightPosition, interplNormal);

                    rasterFragment.colour.red *= ambientL[0]*interplambientM[0] +
                                                 diffuseL[0]*interplspecularMastrial[0]*std::max(interplNormal.dot(lightPosition), 0.0f)+
                                                 specularL[0]*interplspecularMastrial[0]*std::pow(std::max(lightReflect.dot(eyeDir_), 0.0f), this->shinessM)+
                                                 this->emissionM[0];

                    rasterFragment.colour.green *=  ambientL[1]*interplambientM[1] +
                                                    diffuseL[1]*interplspecularMastrial[1]*std::max(interplNormal.dot(lightPosition), 0.0f)+
                                                    specularL[1]*interplspecularMastrial[1]*std::pow(std::max(lightReflect.dot(eyeDir_), 0.0f), this->shinessM)+
                                                    this->emissionM[1];

                    rasterFragment.colour.blue *=  ambientL[2]*interplambientM[2] +
                                                   diffuseL[2]*interplspecularMastrial[2]*std::max(interplNormal.dot(lightPosition), 0.0f)+
                                                   specularL[2]*interplspecularMastrial[2]*std::pow(std::max(lightReflect.dot(eyeDir_), 0.0f), this->shinessM)+
                                                   this->emissionM[2];








                }else{ //ground shading

                    float red = (alpha*vertex0Light[0] + beta*vertex1Light[0] + gamma * vertex2Light[0]);
                    float green = (alpha*vertex0Light[1] + beta*vertex1Light[1] + gamma * vertex2Light[1]);
                    float blue = (alpha*vertex0Light[2] + beta*vertex1Light[2] + gamma * vertex2Light[2]);
                    rasterFragment.colour.red *=red;
                    rasterFragment.colour.green *=green;
                    rasterFragment.colour.blue *= blue;
                }
            }







            fragmentQueue.push_back(rasterFragment);

        } // per pixel

    } // per row
    ProcessFragment();

} // RasteriseTriangle()

// thread fragment
// process a single fragment
void FakeGL::ProcessFragment()
{ // ProcessFragment()
    while(!this->fragmentQueue.empty()){
        auto fragment = this->fragmentQueue.front();
        this->fragmentQueue.pop_front();
        int x = fragment.col;
        int y = fragment.row;
        this->frameBuffer[y][x] = fragment.colour;
    }
} // ProcessFragment()

// standard routine for dumping the entire FakeGL context (except for texture / image)
std::ostream &operator << (std::ostream &outStream, FakeGL &fakeGL)
    { // operator <<
    outStream << "=========================" << std::endl;
    outStream << "Dumping FakeGL Context   " << std::endl;
    outStream << "=========================" << std::endl;


    outStream << "-------------------------" << std::endl;
    outStream << "Vertex Queue:            " << std::endl;
    outStream << "-------------------------" << std::endl;
    for (auto vertex = fakeGL.vertexQueue.begin(); vertex < fakeGL.vertexQueue.end(); vertex++)
        { // per matrix
        outStream << "Vertex " << vertex - fakeGL.vertexQueue.begin() << std::endl;
        outStream << *vertex;
        } // per matrix


    outStream << "-------------------------" << std::endl;
    outStream << "Raster Queue:            " << std::endl;
    outStream << "-------------------------" << std::endl;
    for (auto vertex = fakeGL.rasterQueue.begin(); vertex < fakeGL.rasterQueue.end(); vertex++)
        { // per matrix
        outStream << "Vertex " << vertex - fakeGL.rasterQueue.begin() << std::endl;
        outStream << *vertex;
        } // per matrix


    outStream << "-------------------------" << std::endl;
    outStream << "Fragment Queue:          " << std::endl;
    outStream << "-------------------------" << std::endl;
    for (auto fragment = fakeGL.fragmentQueue.begin(); fragment < fakeGL.fragmentQueue.end(); fragment++)
        { // per matrix
        outStream << "Fragment " << fragment - fakeGL.fragmentQueue.begin() << std::endl;
        outStream << *fragment;
        } // per matrix


    return outStream;
    } // operator <<

// subroutines for other classes
std::ostream &operator << (std::ostream &outStream, vertexWithAttributes &vertex)
    { // operator <<
    std::cout << "Vertex With Attributes" << std::endl;
    std::cout << "Position:   " << vertex.position << std::endl;
    std::cout << "Colour:     " << vertex.colour << std::endl;

	// you

    return outStream;
    } // operator <<

std::ostream &operator << (std::ostream &outStream, screenVertexWithAttributes &vertex) 
    { // operator <<
    std::cout << "Screen Vertex With Attributes" << std::endl;
    std::cout << "Position:   " << vertex.position << std::endl;
    std::cout << "Colour:     " << vertex.colour << std::endl;

    return outStream;
    } // operator <<

std::ostream &operator << (std::ostream &outStream, fragmentWithAttributes &fragment)
    { // operator <<
    std::cout << "Fragment With Attributes" << std::endl;
    std::cout << "Row:        " << fragment.row << std::endl;
    std::cout << "Col:        " << fragment.col << std::endl;
    std::cout << "Colour:     " << fragment.colour << std::endl;

    return outStream;
    } // operator <<


    
    