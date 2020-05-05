#pragma once

#include "BOV.h"


typedef struct cvRect{
    float x, y, w, h;
} cvRect;

/**
 * @brief The ContinuousViewer struct allows to draw particules as a continuoum field
 *
 * @details
 *
 * example usage :
 *
 * ContinuousViewer *viewer = newContinuousViewer();
 *
 * cvPoint somePoints[nPt];
 * drawParticulesContinuous(viewer, somePoints, nPt);
 *
 * freeContinuousViewer(viewer);
 *
 *
 *
 * @note coefficients of interest :  couldAffect function's  k :
 *
 *
 */
typedef struct ContinuousViewer{
    ///region (paticule in this rect are visible) visible. by default, its [-1 1]^2  (x, y = -1   w h = 2)
    cvRect modelViewport;

//    ///screen that is used . by default, its [-1 1]^2  (x, y = -1   w h = 2)
//    cvRect viewViewport;

    /// the grad start at min and finishes at max.  values outside of range get clamped.
    float minVal, maxVal;


    /// private stuff
    GLuint _vao, _vbo, _ibo;
    GLuint _prog, _vert, _frag;
}ContinuousViewer;





ContinuousViewer *newContinuousViewer();
void freeContinuousViewer(ContinuousViewer *v);


typedef struct cvPoint{
    float x, y, val;
} cvPoint;

/**
 * @brief drawParticulesContinuous draw points as a continuous field.
 * @param viewer    the viewer @see newContinuousViewer
 * @param pts       array of points
 * @param nPt       number of points
 *
 * @note to personalise view range, you should edit the viewer's viewport and value range before calling this function
 * @note Values outside range get range get clamped
 *
 * @note Implement a quad tree, splitting if nb pt > M   (here, 32) and depth < 5
 * @note points of the k neighbouging cell are included in a cell as 'affecting it'. Small k (function couldAffect) is an optimization, but can lead to line artefacts
 */

void drawParticulesContinuous(ContinuousViewer *viewer, cvPoint *pts, size_t nPt);
