#include "continuousviewer.h"
#include <assert.h>


ContinuousViewer *newContinuousViewer(){
    ContinuousViewer *v = malloc(sizeof(ContinuousViewer));


    v->modelViewport.x = -1;
    v->modelViewport.y = -1;
    v->modelViewport.w = 2;
    v->modelViewport.h = 2;


    v->maxVal = 1;
    v->minVal = 0;

    //vertex shader

    const char *vertSource = "#version 130\n\
    in mediump vec2 point;\n\
    out mediump vec2 pos;\n\
    uniform mat3 mv;\n\
    void main()\n\
    {\n\
      gl_Position = vec4(mv*vec3(point, 1), 1);\n\
      pos = vec3(mv*vec3(point, 1)).xy;\n\
    }";

    v->_vert = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(v->_vert, 1, &vertSource, NULL);
    glCompileShader(v->_vert);



    //frag shader
    const char *fragSource = "#version 130\n\
    in mediump vec2 pos;\n\
    out mediump vec4 fragColor;\n\
    uniform vec2 valRangeCoefs;\n\
    uniform vec3 ctrlPts[512];\n\
    uniform int nPt;\n\
    void main()\n\
    {\n\
        float value = 0; float sum = 0;float minDist=10000000.0;\n\
        for(int i = 0; i < nPt; i++){\n\
            float w = 1.0 / (pow(10000*dot(ctrlPts[i].xy - pos, ctrlPts[i].xy - pos), 2) + 1);\n\
            sum += w;\n\
            value += w*ctrlPts[i].z;\n\
            minDist = min(minDist, dot(ctrlPts[i].xy - pos, ctrlPts[i].xy - pos));\n\
        }\n\
        value /= sum;\n\
        float val = valRangeCoefs.x*value + valRangeCoefs.y; \n\
// //old palette\n\
//        if(val < 0.25) \n\
//            fragColor = vec4(mix(vec3(72, 142, 202) / 255.0, vec3(72, 187, 70)/255.0, max(val, 0)*4.0), 1);\n\
//        else if(val < 0.5) \n\
//            fragColor = vec4(mix(vec3(72, 187, 70) / 255.0, vec3(240, 120, 41)/255.0, (val-0.25)*4.0), 1);\n\
//        else if(val < 0.75) \n\
//            fragColor = vec4(mix(vec3(240, 120, 41) / 255.0, vec3(210, 31, 40)/255.0, (val-0.5)*4.0), 1);\n\
//        else \n\
//            fragColor = vec4(mix(vec3(210, 31, 40)/255, vec3(165, 47, 90)/255.0, min((val-0.75)*4.0, 1.0)), 1);\n\
////std palette \n\
        fragColor.xyz = vec3(1.5) - 4.0 * abs(val - vec3(0.75, 0.5, 0.25));\n\
        fragColor.a = 1.0-smoothstep(0.8, 1.0, sqrt(minDist)*2.0);\n\
//        fragColor = vec4(vec3(sqrt(minDist)*10.0), 1);\n\
    }";
   
   /* //frag shader with interpolation for missing part
    const char *fragSource = "#version 130\n\
    in mediump vec2 pos;\n\
    out mediump vec4 fragColor;\n\
    uniform vec2 valRangeCoefs;\n\
    uniform vec3 ctrlPts[512];\n\
    uniform int nPt;\n\
    void main()\n\
    {\n\
        float value = 0; float sum = 0;\n\
        for(int i = 0; i < nPt; i++){\n\
            float w = 1.0 / (pow(10000*dot(ctrlPts[i].xy - pos, ctrlPts[i].xy - pos), 2) + 1);\n\
            sum += w;\n\
            value += w*ctrlPts[i].z;\n\
//            minDist = min(minDist, dot(ctrlPts[i].xy - pos, ctrlPts[i].xy - pos));\n\
        }\n\
        value /= sum;\n\
        float val = valRangeCoefs.x*value + valRangeCoefs.y; \n\
// //old palette\n\
//        if(val < 0.25) \n\
//            fragColor = vec4(mix(vec3(72, 142, 202) / 255.0, vec3(72, 187, 70)/255.0, max(val, 0)*4.0), 1);\n\
//        else if(val < 0.5) \n\
//            fragColor = vec4(mix(vec3(72, 187, 70) / 255.0, vec3(240, 120, 41)/255.0, (val-0.25)*4.0), 1);\n\
//        else if(val < 0.75) \n\
//            fragColor = vec4(mix(vec3(240, 120, 41) / 255.0, vec3(210, 31, 40)/255.0, (val-0.5)*4.0), 1);\n\
//        else \n\
//            fragColor = vec4(mix(vec3(210, 31, 40)/255, vec3(165, 47, 90)/255.0, min((val-0.75)*4.0, 1.0)), 1);\n\
////std palette \n\
        fragColor.xyz = vec3(1.5) - 4.0 * abs(val - vec3(0.75, 0.5, 0.25));\n\
        fragColor.a = 1.0;//-smoothstep(0.8, 1.0, sqrt(minDist)*20.0);\n\
//        fragColor = vec4(vec3(sqrt(minDist)*10.0), 1);\n\
    }";*/


    v->_frag = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(v->_frag, 1, &fragSource, NULL);
    glCompileShader(v->_frag);

    //getting compilation error
    GLint erreurCompilation = 0;
    glGetShaderiv(v->_frag, GL_COMPILE_STATUS, &erreurCompilation);


    // S'il y a eu une erreur

    if(erreurCompilation != GL_TRUE)
    {
        printf("frag compile err");
        exit(1);
    }



    //link
    v->_prog = glCreateProgram();
    glAttachShader(v->_prog, v->_vert);
    glAttachShader(v->_prog, v->_frag);
    glLinkProgram(v->_prog);

    GLint erreurLink=0;
    glGetProgramiv(v->_prog, GL_LINK_STATUS, &erreurLink);
    if(erreurLink != GL_TRUE){
        printf("link error\n");
        exit(1);

    }




    glGenVertexArrays(1, &v->_vao);
    glBindVertexArray(v->_vao);

    glGenBuffers(1, &v->_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, v->_vbo);

    float pts[] = {0, 0, 0, 1, 1, 0, 1, 1};
    glBufferData(GL_ARRAY_BUFFER, sizeof(pts), pts, GL_STATIC_DRAW);

    glGenBuffers(1, &v->_ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, v->_ibo);
    unsigned int ids[] = {0, 1, 2, 1, 2, 3 };
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(ids) , ids, GL_STATIC_DRAW);

    glVertexAttribPointer(glGetAttribLocation(v->_prog, "point"), 2, GL_FLOAT, GL_FALSE, 0, (void *)0);

    glEnableVertexAttribArray(0);


    return v;
}

void freeContinuousViewer(ContinuousViewer *v)
{
    glDisableVertexAttribArray(0);

    glDeleteBuffers(1, &v->_ibo);

    glDeleteBuffers(1, &v->_vbo);

    glDeleteVertexArrays(1, &v->_vao);


    glDeleteProgram(v->_prog);
    glDeleteShader(v->_vert);
    glDeleteShader(v->_frag);

    free(v);
}

static int couldAffect(cvPoint pt, cvRect r){
    float k = 6.f;

    //rect contains
    return pt.x >= r.x - r.w  *k && pt.y >= r.y - r.h  *k && pt.x <= r.x+r.w + r.w  *k&& pt.y <= r.y+r.h + r.h  *k;
}


static void drawParticulesContinuousGrid(ContinuousViewer *viewer, const cvPoint *pts, size_t nPt, cvRect rect, int depth)
{

    //filter all the pts that are far
    cvPoint relevant[nPt];
    size_t nRelevantPt=0;
    for(int i = 0; i < nPt; i++){
        if(couldAffect(pts[i], rect))
            relevant[nRelevantPt++] = pts[i];
    }


    //if we have too much pts and we did not recursed forever, split the problem
    if(nRelevantPt > 32 && depth < 5){
        rect.w *= 0.5;
        rect.h *= 0.5;


        drawParticulesContinuousGrid(viewer, relevant, nRelevantPt, rect, depth+1);
        rect.x += rect.w;
        drawParticulesContinuousGrid(viewer, relevant, nRelevantPt, rect, depth+1);
        rect.y += rect.h;
        drawParticulesContinuousGrid(viewer, relevant, nRelevantPt, rect, depth+1);
        rect.x -= rect.w;
        drawParticulesContinuousGrid(viewer, relevant, nRelevantPt, rect, depth+1);

        return;
    }

    assert(nRelevantPt < 512);




  //the viewport matrix
  float mvMat[9] = {rect.w, 0, 0, 0, rect.h, 0, rect.x, rect.y, 1};
  glUniformMatrix3fv(glGetUniformLocation(viewer->_prog, "mv"), 1, GL_FALSE, mvMat);


  //the viewport matrix
  glUniform2f(glGetUniformLocation(viewer->_prog, "valRangeCoefs"), 1.f / (viewer->maxVal - viewer->minVal), - viewer->minVal/ (viewer->maxVal - viewer->minVal));

  //the points
  glUniform3fv(glGetUniformLocation(viewer->_prog, "ctrlPts"), nRelevantPt, relevant);


  glUniform1i(glGetUniformLocation(viewer->_prog, "nPt"), nRelevantPt);


  glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, (void *)0);




}


void drawParticulesContinuous(ContinuousViewer *viewer, cvPoint *pts, size_t nPt){

    cvRect viewport;
    viewport.x = -0.257;
    viewport.y = -0.333;
    viewport.w = 0.257+0.289;
    viewport.h = 0.333+0.384;
    


//      float mvMat[9] = {, 0, 0, 0, 2.f/viewer->viewport.h, 0, , -2.f*viewer->viewport.y/viewer->viewport.h - 1, 1};

      float scaleX = 2.f/viewer->modelViewport.w;
      float scaleY = 2.f/viewer->modelViewport.h;

      float offX = -2.f*viewer->modelViewport.x/viewer->modelViewport.w - 1;
      float offY = -2.f*viewer->modelViewport.y/viewer->modelViewport.h - 1;

    //filter all the pts that are far
    cvPoint transformed[nPt];
    for(int i = 0; i < nPt; i++){
        transformed[i].val = pts[i].val;
        transformed[i].x= pts[i].x * scaleX + offX;
        transformed[i].y= pts[i].y* scaleY + offY;
    }

    glDisable(GL_CULL_FACE);
    glBindVertexArray(viewer->_vao);
    glUseProgram(viewer->_prog);

    drawParticulesContinuousGrid(viewer, transformed, nPt, viewport, 0);
    glBindVertexArray(0);
}




//old, with delaunay
//ContinuousViewer *newContinuousViewer(){
//    ContinuousViewer *v = malloc(sizeof(ContinuousViewer));


//    v->viewport.x = -1;
//    v->viewport.y = -1;
//    v->viewport.w = 2;
//    v->viewport.h = 2;


//    v->maxVal = 1;
//    v->minVal = 0;

//    //vertex shader

//    const char *vertSource = "#version 130\n\
//    in mediump vec2 point;\n\
//    in mediump float texcoord;\n\
//    out mediump float value;\n\
//    uniform mat3 mv;\n\
//    void main()\n\
//    {\n\
//      gl_Position = vec4(mv*vec3(point, 1), 1);\n\
//      value = texcoord;\n\
//    }";

//    v->_vert = glCreateShader(GL_VERTEX_SHADER);
//    glShaderSource(v->_vert, 1, &vertSource, NULL);
//    glCompileShader(v->_vert);



//    //frag shader
//    const char *fragSource = "#version 130\n\
//    in mediump float value;\n\
//    out mediump vec4 fragColor;\n\
//    uniform vec2 valRangeCoefs;\n\
//    void main()\n\
//    {\n\
//        float val = valRangeCoefs.x*value + valRangeCoefs.y; \n\
//        if(val < 0.25) \n\
//            fragColor = vec4(mix(vec3(72, 142, 202) / 255.0, vec3(72, 187, 70)/255.0, max(val, 0)*4.0), 1);\n\
//        else if(val < 0.5) \n\
//            fragColor = vec4(mix(vec3(72, 187, 70) / 255.0, vec3(240, 120, 41)/255.0, (val-0.25)*4.0), 1);\n\
//        else if(val < 0.75) \n\
//            fragColor = vec4(mix(vec3(240, 120, 41) / 255.0, vec3(210, 31, 40)/255.0, (val-0.5)*4.0), 1);\n\
//        else \n\
//            fragColor = vec4(mix(vec3(210, 31, 40)/255, vec3(165, 47, 90)/255.0, min((val-0.75)*4.0, 1.0)), 1);\n\
////        fragColor = vec4(vec3(val), 1);\n\
//    }";


//    v->_frag = glCreateShader(GL_FRAGMENT_SHADER);
//    glShaderSource(v->_frag, 1, &fragSource, NULL);
//    glCompileShader(v->_frag);


//    //link
//    v->_prog = glCreateProgram();
//    glAttachShader(v->_prog, v->_vert);
//    glAttachShader(v->_prog, v->_frag);
//    glLinkProgram(v->_prog);


//    glGenVertexArrays(1, &v->_vao);
//    glBindVertexArray(v->_vao);

//    glGenBuffers(1, &v->_vbo);
//    glBindBuffer(GL_ARRAY_BUFFER, v->_vbo);

//    glGenBuffers(1, &v->_ibo);
//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, v->_ibo);

//    glVertexAttribPointer(glGetAttribLocation(v->_prog, "point"), 2, GL_FLOAT, GL_FALSE, sizeof(cvPoint), (void *)0);
//    glVertexAttribPointer(glGetAttribLocation(v->_prog, "texcoord"), 1, GL_FLOAT, GL_FALSE, sizeof(cvPoint), (void *)(2 * sizeof(float)));

//    glEnableVertexAttribArray(0);
//    glEnableVertexAttribArray(1);


//    return v;
//}

//void freeContinuousViewer(ContinuousViewer *v)
//{
//    glDisableVertexAttribArray(1);
//    glDisableVertexAttribArray(0);

//    glDeleteBuffers(1, &v->_ibo);

//    glDeleteBuffers(1, &v->_vbo);

//    glDeleteVertexArrays(1, &v->_vao);


//    glDeleteProgram(v->_prog);
//    glDeleteShader(v->_vert);
//    glDeleteShader(v->_frag);

//    free(v);
//}

//void drawParticulesContinuous(ContinuousViewer *viewer, cvPoint *pts, size_t nPt)
//{

//    double positions[nPt*2];

//    for(int i = 0; i < nPt; i++){
//        positions[2*i + 0] = pts[i].x;
//        positions[2*i + 1] = pts[i].y;
//    }

//    Delaunator *d = make_delau(positions, nPt*2);



//    //update vertices
//    glBindBuffer(GL_ARRAY_BUFFER, viewer->_vbo);
//    glBufferData(GL_ARRAY_BUFFER, sizeof(cvPoint) * nPt, pts, GL_DYNAMIC_DRAW);
////    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_DYNAMIC_DRAW);

//    //update the indices
//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, viewer->_ibo);
//    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * d->triangles.size, d->triangles.v, GL_DYNAMIC_DRAW);


//    glDisable(GL_CULL_FACE);



//    glBindVertexArray(viewer->_vao);

//  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
//  glClear(GL_COLOR_BUFFER_BIT);
//  glUseProgram(viewer->_prog);

//  //the viewport matrix
//  float mvMat[9] = {2.f/viewer->viewport.w, 0, 0, 0, 2.f/viewer->viewport.h, 0, -2.f*viewer->viewport.x/viewer->viewport.w - 1, -2.f*viewer->viewport.y/viewer->viewport.h - 1, 1};
//  glUniformMatrix3fv(glGetUniformLocation(viewer->_prog, "mv"), 1, GL_FALSE, mvMat);


//  //the viewport matrix
//  glUniform2f(glGetUniformLocation(viewer->_prog, "valRangeCoefs"), 1.f / (viewer->maxVal - viewer->minVal), - viewer->minVal/ (viewer->maxVal - viewer->minVal));


//  glDrawElements(GL_TRIANGLES, d->triangles.size, GL_UNSIGNED_INT, (void *)0);

//  glBindVertexArray(0);

//  free_delau(d);


//}
