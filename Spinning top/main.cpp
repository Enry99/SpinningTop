//#define DATA_OUTPUT
#define SLEEP
#ifndef SLEEP
bool goDraw = true;
#endif
#if defined(WIN32)|| defined(_WIN32)
#define _USE_MATH_DEFINES
#endif

#include "Quaternion.h"
#include <FL/gl.h>
#include <FL/glu.h>
#include <FL/glut.H>
#include <FL/x.H>
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Multiline_Output.H>
#include <FL/Fl_Check_Button.H>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <thread>
#include <time.h>
#include <cstdarg>

#include <exception>
#include <stdexcept>

//extern functions/variables (defined in window.cpp)
extern int CreateMyWindow(int argc, char** argv);
extern void activateButtons();
extern Fl_Check_Button* enable_gravity;
extern Fl_Check_Button* enable_file_output;
extern Fl_Check_Button* persistent_trail;
extern double FPS_display_width, FPS_display_height;
extern double
cameraPos[3],
cameraUp[3],
cameraFront[3];

//variables declarations
constexpr double cyl_radius_ratio = 0.1;
constexpr double opengl_top_scale_factor = 25; //temporary: axis lenght, floor and camera distance should be adjusted to be 1:1
constexpr double trail_duration = 60;
constexpr double trail_step = 0.001;
GLfloat amb_light[4]{ 0.1, 0.1, 0.1, 1.0 };
GLfloat diffuse[4]{ 0.9, 0.9, 0.9, 1.0 };
GLfloat specular[4]{ 0.7, 0.7, 0.7, 1.0 };
float lightAngle = 0.0, lightHeight = 10;
float L_floor = 3;
/* Variable controlling various rendering modes. */
int stencilReflection = 1, stencilShadow = 1, offsetShadow = 1;
int renderShadow = 1, renderReflection = 1;
int directionalLight = 0;

double xRotated, a1, a2, a3; //rotation angle and axis components
double values[15]{}; //values tied to sliders
double height_cone1;
double height_cone2;
double cones_radius;
double cyl_height;
double cyl_radius;
bool run_animation = false; //run starts the animation when the first frame is ready
bool pause_animation = false; //pause is controlled by button
int light_state = 0;
GLfloat myFloorPlane[4] = { 0.,0.,1.,0. };
GLfloat floorShadow[4][4]{};
GLfloat lightPosition[4]{};
#define DARK_TOP true
bool
DRAW_SUN = true,
DRAW_FLOOR = true,
SHOW_FPS = true,
DRAW_TRAIL = true,
DRAW_AXIS = true,
DRAW_TOP = true,
DRAW_VECTORS = false,
hide_vectors, //hides vectors when animation is stopped
rainbow = false;
std::chrono::steady_clock::time_point 
idle_prev_time = std::chrono::steady_clock::now(),
previous_time, 
FPS_previous_time; //clock for fps count, is updated every second
int frames_counter = 0, FPS = 0;
std::vector<std::array<double, 4>> rotation_data;
std::vector<std::array<double, 3>> trail_points;
std::vector<std::array<double, 4>> W_vectors;
std::vector<std::array<double, 4>> L_vectors;
std::vector<std::array<double, 4>> T_vectors;
std::vector<std::array<double, 3>> rcm_data;
std::vector<std::array<double, 3>> L_lab_data;
std::vector<std::array<double, 3>> L_top_data;
std::vector<std::array<double, 3>> w_lab_data;
std::vector<std::array<double, 3>> w_top_data;
std::vector<std::array<double, 3>> theta_phi_phidot_data;
std::vector<double> Energy_data;
std::ofstream
L_stream,
w_stream,
Energy_stream,
theta_phi_phidot_stream,
Rcm_stream;



//Runge-Kutta internal variables and functions
double
w0[3]{},
wt[3]{},
field1[3]{},
field2[3]{},
field3[3]{},
eval0[3]{},
eval1[3]{},
eval2[3]{},
eval3[3]{},
FORCE[3]{},
FORCE_top[3]{},
TORQUE[3]{},
Rcm_top[3]{},
Rcm[3]{}, // Rcm lab frame
mass,
init_angle,
t,
tf,
delta_t,
I1,
I2,
I3;
Quaternion q;
Quaternion Rcm_top_quat;
double accumulator = 0;
//vector field components
double fx(const double* const& vect) { return (TORQUE[0] - vect[1] * vect[2] * (I3 - I2)) / I1; }
double fy(const double* const& vect) { return (TORQUE[1] - vect[2] * vect[0] * (I1 - I3)) / I2; }
double fz(const double* const& vect) { return (TORQUE[2] - vect[0] * vect[1] * (I2 - I1)) / I3; }
double (*F[3])(const double* const&) { fx, fy, fz };


//functions declarations
void keyboardFunction(unsigned char, int, int);
void myCone(double, double, double, int, bool = false);
void shadowMatrix(GLfloat[4][4], GLfloat[4], GLfloat[4]);
void drawAxis();
void drawTrail();
void drawVectors();
void drawFPS(std::string);
void drawTop(bool = false);
void drawFloor();
void displaySpinningTop();
void setInitialConditions();
void evolve();
void startAlgorithm();
void animationLoop();

template<class X, typename std::enable_if<std::is_arithmetic<X>::value, void>* = nullptr>
X vector_modulus(const X vect[3])
{
    return sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2]);
}

template<class X, typename std::enable_if<std::is_arithmetic<X>::value, void>* = nullptr>
void cross_product(const X vector1[3], const X vector2[3], X result[3])
{
    result[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1];
    result[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2];
    result[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0];
}

template<class X, typename std::enable_if<std::is_arithmetic<X>::value, void>* = nullptr>
X dot_product(const X * vector1, const X * vector2, int size = 3)
{
    X result{};
    for (int i = 0; i < size; ++i) result += vector1[i] * vector2[i];
    return result;
}


/*
void printv(va_list args, const char* format) 
{
    constexpr int LEN = 8192;
    char buf[LEN];
    char* ch = buf;
    vsnprintf(buf, LEN, format, args);
    while (*ch) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *ch++);
}

void print(const char* format, ...)
{
    va_list args;
    va_start(args, format);
    printv(args, format);
    va_end(args);
}
*/



void keyboardFunction(unsigned char key, int, int)
{

    switch (key)
    {
    case '+': light_state = 1;    break;
    case '-': light_state = -1;   break;
    case 't': DRAW_TRAIL = !DRAW_TRAIL;    break;
    case 'r': rainbow = !rainbow;  break;
    case 'v': DRAW_VECTORS = !DRAW_VECTORS;    break;
    case 'a': DRAW_AXIS  = !DRAW_AXIS; break;
    case 'p': pause_animation = !pause_animation; break;
    case 'f': DRAW_FLOOR = !DRAW_FLOOR; break;
    case 's': DRAW_SUN = !DRAW_SUN; break;
    case ' ': SHOW_FPS = !SHOW_FPS; break;
    case 'h': DRAW_TOP = !DRAW_TOP;  break;
    case 27: exit(0);   break;// exit program when [ESC] key presseed
    default:   break;
    }

}

void myCone(double base_radius, double top_radius, double height, int slices, bool dark) {

    glFrontFace(GL_CW); //gl_quad_strip has clockwise winding 
    glShadeModel(GL_FLAT);
    glBegin(GL_QUAD_STRIP);
    float step = 2 * M_PI / slices;
    float diagonal = sqrt((base_radius - top_radius) * (base_radius - top_radius) + height * height);
    float sine = (base_radius - top_radius) / diagonal;
    float cosine = height / diagonal;
    bool color = false; //alternates slices color
    for (float k = 0; k <= 2. * M_PI; k += step)
    {
        dark ? glColor4f(0., 0., 0., 0.5) : color ? glColor3f(0.2, 0.2, 1) : glColor3f(1.0, 0.0, 0.0);
        glNormal3f(cosine * cos(k), cosine * sin(k), sine);
        glVertex3f(base_radius * cos(k), base_radius * sin(k), 0);
        glVertex3f(top_radius * cos(k), top_radius * sin(k), height);
        color = !color; //inverts color for the next slice
    }
    glEnd();
    glShadeModel(GL_SMOOTH);
    glFrontFace(GL_CCW);
}

void shadowMatrix(GLfloat shadowMat[4][4], GLfloat groundplane[4], GLfloat lightpos[4])
{
    /* Create a matrix that will project the desired shadow.
    projects the light in lightpos on groundplane*/

    /* Find dot product between light position vector and ground plane normal. */
    GLfloat dot = dot_product(groundplane, lightpos, 4);

    enum { X, Y, Z, W };

    shadowMat[0][0] = dot - lightpos[X] * groundplane[X];
    shadowMat[1][0] = 0.f - lightpos[X] * groundplane[Y];
    shadowMat[2][0] = 0.f - lightpos[X] * groundplane[Z];
    shadowMat[3][0] = 0.f - lightpos[X] * groundplane[W];

    shadowMat[X][1] = 0.f - lightpos[Y] * groundplane[X];
    shadowMat[1][1] = dot - lightpos[Y] * groundplane[Y];
    shadowMat[2][1] = 0.f - lightpos[Y] * groundplane[Z];
    shadowMat[3][1] = 0.f - lightpos[Y] * groundplane[W];

    shadowMat[X][2] = 0.f - lightpos[Z] * groundplane[X];
    shadowMat[1][2] = 0.f - lightpos[Z] * groundplane[Y];
    shadowMat[2][2] = dot - lightpos[Z] * groundplane[Z];
    shadowMat[3][2] = 0.f - lightpos[Z] * groundplane[W];

    shadowMat[X][3] = 0.f - lightpos[W] * groundplane[X];
    shadowMat[1][3] = 0.f - lightpos[W] * groundplane[Y];
    shadowMat[2][3] = 0.f - lightpos[W] * groundplane[Z];
    shadowMat[3][3] = dot - lightpos[W] * groundplane[W];
}

void drawAxis(){
    
    if (DRAW_AXIS)
    {
        double length = 3;  //3 metri
        glLineWidth(3);
        glBegin(GL_LINES);
        glNormal3f(0.4, 0, 1);
        glColor3f(1.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(length, 0.0, 0.0);

        glColor3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, length, 0.0);

        glColor3f(0.0, 0.0, 1.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, length);

        glEnd();
        glLineWidth(1);

        glColor3f(1.0, 1.0, 1.0);
        glRasterPos3d(1.1 * length, 0, 0);
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'X');
        glRasterPos3d(0, 1.1 * length, 0);
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'Y');
        glRasterPos3d(0, 0, 1.1 * length);
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'Z');
    }
}

void drawTrail()
{
    if(trail_points.size()>=2 && DRAW_TRAIL)
    {
        int R=255, G = 0, B = 0;
        glLineWidth(2);
        glBegin(GL_LINE_STRIP);
        glNormal3f(0.4, 0, 1);
        if(!rainbow) glColor3f(0.5*255./255., 0.5*94./255., 0.5*19./255.);
        int size = trail_points.size();
        int increment = int(trail_step / delta_t) > 0 ? int(trail_step / delta_t) : 1;
        if (rainbow && !persistent_trail->value() && int(size - trail_duration / delta_t) > 0) 
        {   //correction for color mismatch when trail is not persistent
            for (int i = 0; i < size - trail_duration / delta_t; i += increment)
            {
                if (R > 0 && B == 0) ++G, --R;
                if (G > 0 && R == 0) ++B, --G;
                if (B > 0 && G == 0) ++R, --B;
            }
        }
        for(int i = persistent_trail->value() || int(size - trail_duration/delta_t)< 0  ? 0 : size - trail_duration / delta_t; i<size; i+=increment)
        {
            if(rainbow)
            {
                glColor3f(R/255., G/255., B/255.);
                if(R>0 && B == 0) ++G, --R;
                if(G>0 && R == 0) ++B, --G;
                if(B>0 && G == 0) ++R, --B;
            }

            glVertex3f(trail_points[i][0], trail_points[i][1], trail_points[i][2]);
        }
        glEnd();
        glLineWidth(1);

    }
}

void drawVectors()
{
       
    if(!W_vectors.empty() && DRAW_VECTORS && !hide_vectors)//DRAW_VECTORS is from keyboard, hide_vectors is when idle
    {
        GLUquadricObj* vector_W_body=gluNewQuadric();
        GLUquadricObj* vector_L_body=gluNewQuadric();
        GLUquadricObj* vector_T_body=gluNewQuadric();

        GLdouble arrow_radius = 0.03, arrow_height = 0.1, stick_radius = 0.006;
        GLint vectorSlices = 15, vectorStacks = 15;

        glPushMatrix(); //PUSH
        glRotated(W_vectors.back()[0], W_vectors.back()[1], W_vectors.back()[2], W_vectors.back()[3]);
        glColor3f(0.0, 1., 0.1);
        gluCylinder(vector_W_body, stick_radius, stick_radius, opengl_top_scale_factor * (height_cone1 + height_cone2 + 1.1*cyl_height),vectorSlices,vectorStacks); //arrow body
        glTranslated(0.0, 0.0, opengl_top_scale_factor * (height_cone1 + height_cone2 + 1.1*cyl_height));
        glutSolidCone(arrow_radius, arrow_height, vectorSlices, vectorStacks); //arrow tip
        glPopMatrix(); //POP

        glPushMatrix(); //PUSH
        glRotated(L_vectors.back()[0], L_vectors.back()[1], L_vectors.back()[2], L_vectors.back()[3]);
        glColor3f(0.0, 0.0, 1.0);
        gluCylinder(vector_L_body, stick_radius*1.1, stick_radius*1.1, opengl_top_scale_factor * (1. + arrow_height + 0.006) * (height_cone1 + height_cone2 + cyl_height),vectorSlices,vectorStacks);
        glTranslated(0.0, 0.0, opengl_top_scale_factor * (1. + arrow_height + 0.007) * (height_cone1 + height_cone2 + cyl_height));
        glutSolidCone(arrow_radius, arrow_height, vectorSlices, vectorStacks);
        glTranslated(0.0, 0.0, arrow_height);
        glRotated(L_vectors.back()[0], -L_vectors.back()[1], -L_vectors.back()[2], -L_vectors.back()[3]);
        glRotated(T_vectors.back()[0], T_vectors.back()[1], T_vectors.back()[2], T_vectors.back()[3]);
        glColor3f(1, 1, 0);
        gluCylinder(vector_T_body, stick_radius, stick_radius, opengl_top_scale_factor / 5 * (height_cone1 + height_cone2 + cyl_height), vectorSlices, vectorStacks);
        glTranslated(0.0, 0.0, opengl_top_scale_factor / 5 * (height_cone1 + height_cone2 + cyl_height));
        glutSolidCone(arrow_radius, arrow_height, vectorSlices, vectorStacks);
        glPopMatrix(); //POP

    }
}

void drawFPS(std::string text)
{

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.0, FPS_display_width, FPS_display_height, 0.0);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glColor3f(0, 1, 0);
    glRasterPos2f(0.007, 0.015);
    for (auto i : text) glutBitmapCharacter(GLUT_BITMAP_9_BY_15, i);

    glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
}

void drawTop(bool dark)
{
    GLint slices = 30, stacks = 30;
    height_cone1 = values[10];
    height_cone2 = values[11];
    cones_radius = values[12];
    cyl_height = values[13];
    cyl_radius = cyl_radius_ratio * cones_radius;

    glPushMatrix(); //initial push
    dark ? glColor4f(0., 0., 0., 0.5) : glColor3f(0.8,0.2,0.1);
    //upper cone
    glTranslated(0.0, 0.0, opengl_top_scale_factor * height_cone1);       
    myCone(opengl_top_scale_factor * cones_radius, opengl_top_scale_factor *cyl_radius, opengl_top_scale_factor * height_cone2, slices, dark);
    
    glPushMatrix(); //push cylinder
     glTranslated(0.0, 0.0, opengl_top_scale_factor * height_cone2);
     GLUquadricObj* cylinder = gluNewQuadric();
     GLUquadricObj* circle = gluNewQuadric();
     gluCylinder(cylinder, opengl_top_scale_factor*cyl_radius, opengl_top_scale_factor*cyl_radius, opengl_top_scale_factor * cyl_height, slices, stacks);
     glTranslated(0.,0., opengl_top_scale_factor *cyl_height);
     gluDisk(circle,0, opengl_top_scale_factor*cyl_radius,slices,1);
    glPopMatrix(); //pop cylinder

    glPushMatrix(); //push lower cone
     glScaled(1.0, 1.0, -1.0); //cone points downwards
     glFrontFace(GL_CW); //because of scaling -1 on z
     glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
     glutSolidCone(opengl_top_scale_factor *cones_radius, opengl_top_scale_factor *height_cone1,slices,stacks);
     glFrontFace(GL_CCW);
    glPopMatrix(); //pop lower cone

    glPopMatrix(); //pop initial push
}

void drawFloor()
{
    glClearColor(0.,0.,0.,0.);

#ifdef monocromatic //no reflections
    glDisable(GL_LIGHTING);
    glColor3f(170./255, 170./255,170./255);
    glBegin(GL_QUADS);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glNormal3f(0, 0, 1);
    glVertex3f(2*L_floor, 2*L_floor, 0); glVertex3f(2*L_floor, -2*L_floor, 0);
    glVertex3f(-2*L_floor, -2*L_floor, 0); glVertex3f(-2*L_floor, 2*L_floor, 0);
    glEnd();
#endif

    GLfloat floorVertices[4][3] = {
    { -2 * L_floor, -2 * L_floor, 0 },
    { 2 * L_floor,  -2 * L_floor, 0 },
    { 2 * L_floor,  2 * L_floor, 0},
    { -2 * L_floor,2 * L_floor, 0 }
    };

    //if(amb_light[0]>0.09) glDisable(GL_LIGHTING);
    glBegin(GL_QUADS);
      glNormal3f(0, 0, 1);
      //glTexCoord2f(0.0, 0.0);
      glVertex3fv(floorVertices[0]);
      //glTexCoord2f(0.0, 16.0);
      glVertex3fv(floorVertices[1]);
      //glTexCoord2f(16.0, 16.0);
      glVertex3fv(floorVertices[2]);
      //glTexCoord2f(16.0, 0.0);
      glVertex3fv(floorVertices[3]);
    glEnd();
    //if (amb_light[0] > 0.09) glEnable(GL_LIGHTING);

    glClearColor(0., 0., 0., 0.);
}

void displaySpinningTop()
{
    if (!frames_counter) FPS_previous_time = std::chrono::steady_clock::now();
    
    // clear the drawing buffer.
    if ((stencilReflection && renderReflection) || (stencilShadow && renderShadow)) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    } else {
    /* Avoid clearing stencil when not using it. */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }
 
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(cameraPos[0], cameraPos[1], cameraPos[2]+L_floor, 0, 0, L_floor, cameraUp[0], cameraUp[1], cameraUp[2]);


    if (light_state != 0) { for (int i = 0; i<3 ; ++i) amb_light[i] += 0.1*light_state; light_state=0; }
    glEnable(GL_LIGHTING);
    glLightfv(GL_LIGHT0, GL_AMBIENT, amb_light);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
    //glEnable(GL_COLOR_MATERIAL);
    //glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    glEnable(GL_LIGHT0);
    /* Reposition the light source. */
    lightPosition[0] = 6*cos(lightAngle);
    lightPosition[1] = 6*sin(lightAngle);
    lightPosition[2] = lightHeight;
    if (directionalLight) lightPosition[3] = 0.0;
    else lightPosition[3] = 1.0;
    /* Tell GL new light source position. */
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

    if (DRAW_FLOOR) 
    {
        shadowMatrix(floorShadow, myFloorPlane, lightPosition); //builds floorShadow
        if (renderReflection && DRAW_TOP) 
        {
            if (stencilReflection) 
            {
                /* We can eliminate the visual "artifact" of seeing the "flipped"
                top underneath the floor by using stencil.  The idea is
                draw the floor without color or depth update but so that
                a stencil value of one is where the floor will be.  Later when
                rendering the spinning top reflection, we will only update pixels
                with a stencil value of 1 to make sure the reflection only
                lives on the floor, not below the floor. */

                /* Don't update color or depth. */
                glDisable(GL_DEPTH_TEST);
                glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

                /* Draw 1 into the stencil buffer. */
                glEnable(GL_STENCIL_TEST);
                glStencilOp(GL_REPLACE, GL_REPLACE, GL_REPLACE);
                glStencilFunc(GL_ALWAYS, 1, 0xffffffff);

                /* Now render floor; floor pixels just get their stencil set to 1. */
                drawFloor();

                /* Re-enable update of color and depth. */
                glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
                glEnable(GL_DEPTH_TEST);

                /* Now, only render where stencil is set to 1. */
                glStencilFunc(GL_EQUAL, 1, 0xffffffff);  /* draw if ==1 */
                glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
            }

            glPushMatrix(); //push reflection

            /* The critical reflection step: Reflect top  through the floor
                (the Z=0 plane) to make a relection. */
            glScalef(1.0, 1.0, -1.0);

            /* Reflect the light position. */
            glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
            //When glLight() is called its position or direction is transformed by the current modelview matrix

            /* To avoid our normals getting reversed and hence botched lighting
            on the reflection, turn on normalize.  */
            glEnable(GL_NORMALIZE);
            glCullFace(GL_FRONT);

            /* Draw the reflected top. */
            glPushMatrix();
             glTranslated(0, 0, L_floor);
             glRotated(xRotated, a1, a2, a3);
             drawTop();
            glPopMatrix();

            /* Disable normalize again and re-enable back face culling. */
            glDisable(GL_NORMALIZE);
            glCullFace(GL_BACK);

            glPopMatrix(); //pop reflection

            /* Switch back to the unreflected light position. */
            glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
            if (stencilReflection) glDisable(GL_STENCIL_TEST);
        }

        /* Back face culling will get used to only draw either the top or the
            bottom floor.  This let's us get a floor with two distinct
            appearances.  The top floor surface is reflective and kind of red.
            The bottom floor surface is not reflective and blue. */

        /* Draw "bottom" of floor in blue. */
        glFrontFace(GL_CW);  /* Switch face orientation. */
        glColor4f(0.1, 0.1, 0.7, 1.0);
        drawFloor();
        glFrontFace(GL_CCW);

        if (renderShadow) 
        {
            if (stencilShadow) 
            {
                /* Draw the floor with stencil value 3.  This helps us only
                    draw the shadow once per floor pixel (and only on the
                    floor pixels). */
                glEnable(GL_STENCIL_TEST);
                glStencilFunc(GL_ALWAYS, 3, 0xffffffff);
                glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
            }
        }

        /* Draw "top" of floor.  Use blending to blend in reflection. */
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(0.7, 0.0, 0.0, 0.3);
        glColor4f(1.0, 1.0, 1.0, 0.3);
        drawFloor();
        glDisable(GL_BLEND);

        if (renderShadow && DRAW_TOP) 
        {

            /* Render the projected shadow. */

            if (stencilShadow) {

                /* Now, only render where stencil is set above 2 (ie, 3 where
                the top floor is).  Update stencil with 2 where the shadow
                gets drawn so we don't redraw (and accidently reblend) the
                shadow). */
                glStencilFunc(GL_LESS, 2, 0xffffffff);  /*   draw if ==1 */
                glStencilOp(GL_REPLACE, GL_REPLACE, GL_REPLACE);
            }


            /* To eliminate depth buffer artifacts, we use polygon offset
            to raise the depth of the projected shadow slightly so
            that it does not depth buffer alias with the floor. */
            if (offsetShadow) glEnable(GL_POLYGON_OFFSET_FILL);


            /* Render 50% black shadow color on top of whatever the
                floor appareance is. */
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glDisable(GL_LIGHTING);  /* Force the 50% black. */
            glColor4f(0.0, 0.0, 0.0, 0.5);

            /* Project the shadow. */
            glPushMatrix(); //push shadow       
             glMultMatrixf((GLfloat*)floorShadow);
             glTranslated(0, 0, L_floor);
             glRotated(xRotated, a1, a2, a3);
             drawTop(DARK_TOP);
            glPopMatrix(); //pop shadow

            glDisable(GL_BLEND);
            glEnable(GL_LIGHTING);

            if (offsetShadow) glDisable(GL_POLYGON_OFFSET_FILL);
            if (stencilShadow) glDisable(GL_STENCIL_TEST);
        }
    }

    /* Draw "actual" top, not its reflection. */
    glPushMatrix(); //push1
    glTranslated(0,0,L_floor);
    drawAxis();
    drawTrail();
    glPushMatrix(); //push2
    glRotated(xRotated,a1,a2,a3);
    drawVectors(); 
    if( DRAW_TOP) drawTop();
    glPopMatrix(); //pop2

    glPopMatrix(); //pop1
        

    glPushMatrix(); //push light
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glColor3f(1.0, 1.0, 0.0);
    if (directionalLight) 
    {
        /* Draw an arrowhead. */
        //this sections is not really functional
        glDisable(GL_CULL_FACE);
        glTranslatef(lightPosition[0], lightPosition[1], lightPosition[2]);
        glRotatef(lightAngle * -180.0 / M_PI, 0, 1, 0);
        glRotatef(atan(lightHeight/6) * 180.0 / M_PI, 1, 0, 0);
        glBegin(GL_TRIANGLE_FAN);
            glVertex3f(0, 0, 0);
            glVertex3f(2, 1, 1);
            glVertex3f(2, -1, 1);
            glVertex3f(2, -1, -1);
            glVertex3f(2, 1, -1);
            glVertex3f(2, 1, 1);
        glEnd();
        /* Draw a white line from light direction. */
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);
            glVertex3f(0, 0, 0);
            glVertex3f(5, 0, 0);
        glEnd();
        glEnable(GL_CULL_FACE);
        //glEnable(GL_COLOR_MATERIAL);
    } else {
        /* Draw a yellow ball at the light source. */
        glTranslatef(lightPosition[0], lightPosition[1], lightPosition[2]);
        if(DRAW_SUN) glutSolidSphere(.8, 20, 20);
    }
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glPopMatrix();//pop light
    //glPopMatrix();


    ++frames_counter;
    if((std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - FPS_previous_time)).count()*1.e-6 >=1) //updated every second
    {
        FPS = frames_counter;
        frames_counter = 0;
    }
    if(SHOW_FPS) drawFPS( std::string("FPS: ") + std::to_string( FPS ) );

    glutSwapBuffers();
    //glFlush();
}

void setInitialConditions()
{
    mass = values[0],
    delta_t = values[1],
    tf = values[2],
    init_angle = values[3],
    w0[0] = values[4],
    w0[1] = values[5],
    w0[2] = values[6],
    height_cone1 = values[10],
    height_cone2 = values[11],
    cones_radius = values[12],
    cyl_height = values[13],
    cyl_radius = cyl_radius_ratio * cones_radius;

    double V1 = M_PI * pow(cones_radius, 2) * height_cone1 / 3;
    double V2 = M_PI * (pow(cones_radius, 2) + pow(cyl_radius, 2) + cones_radius * cyl_radius) * height_cone2 / 3;
    double V3 = M_PI * pow(cyl_radius, 2) * cyl_height;
    double V = V1 + V2 + V3;
    double D = mass / V;
    double mass1 = D * V1;
    double mass2 = D * V2;
    double mass3 = D * V3;
    double Iz1 = D * M_PI / 10 * pow(cones_radius, 4) * height_cone1;
    double Iz2 = D * M_PI / 10 * (pow(cones_radius, 5) - pow(cyl_radius, 5)) / (cones_radius - cyl_radius) * height_cone2;
    double Iz3 = D * M_PI / 2 * pow(cyl_radius, 4) * cyl_height;
    double Iy1 = Iz1 / 2 + D * M_PI / 5 * pow(cones_radius, 2) * pow(height_cone1, 3);
    double Iy2 = Iz2 / 2 + D * M_PI * (
        (pow(height_cone1 + height_cone2, 5) - pow(height_cone1, 5)) / 5 * pow((cones_radius - cyl_radius) / height_cone2, 2)
        - (pow(height_cone1 + height_cone2, 4) - pow(height_cone1, 4)) / 2 * (cones_radius * (cones_radius - cyl_radius) / height_cone2 + pow((cones_radius - cyl_radius) / height_cone2, 2) * height_cone1)
        + (pow(height_cone1 + height_cone2, 3) - pow(height_cone1, 3)) / 3 * pow(cones_radius + (cones_radius - cyl_radius) / height_cone2 * height_cone1, 2));
    double Iy3 = Iz3 / 2 + D * M_PI * pow(cyl_radius, 2) * (pow(height_cone1 + height_cone2 + cyl_height, 3) - pow(height_cone1 + height_cone2, 3)) / 3;

    I3 = Iz1 + Iz2 + Iz3;
    I1 = I2 = Iy1 + Iy2 + Iy3;

    double h_cm_2 = (pow(cones_radius, 2) + 2 * cones_radius * cyl_radius + 3 * pow(cyl_radius, 2)) * height_cone2 / (4 * (pow(cones_radius, 2) + cones_radius * cyl_radius + pow(cyl_radius, 2)));

    Rcm_top[0] = 0;
    Rcm_top[1] = 0;
    Rcm_top[2] = (mass1 * 3. / 4. * height_cone1 + mass2 * (height_cone1 + h_cm_2) + mass3 * (height_cone1 + height_cone2 + cyl_height / 2)) / mass;

    q = Quaternion(cos(init_angle / 2), sin(init_angle / 2), 0, 0);
    Rcm_top_quat = Quaternion(Rcm_top);

    if(enable_file_output->value())
    {
        std::ofstream initial_values_stream("initial_values.txt");
        const char* slider_names[15] = {
            "mass = ",
            "dt = ",
            "tf = ",
            "theta0 = ",
            "w1 = ",
            "w2 = ",
            "w3 =",
            "Fx = ",
            "Fy = ",
            "Fz = ",
            "h1 = ",
            "h2 = ",
            "cone radius = ",
            "cyl height = ",
            "FPS max"
        };
        for (int i = 0; i < 15; ++i)
        {
            initial_values_stream << slider_names[i] << values[i] << '\n';
        }
        initial_values_stream << "l_cm = " << Rcm_top[2] << '\n';
        initial_values_stream << "I1=I2 = : " << I1 << '\n';
        initial_values_stream << "I3 = " << I3 << '\n';
        initial_values_stream << "Density = " << D << '\n';
        initial_values_stream.clear();
        initial_values_stream.close();
        L_stream.open("Ltop_Llab.txt");
        w_stream.open("wtop_wlab.txt");
        Rcm_stream.open("Rcm.txt");
        Energy_stream.open("Energy.txt");
        theta_phi_phidot_stream.open("theta_phi_phidot.txt");
    }

    rotation_data.reserve(tf/delta_t);
    trail_points.reserve(tf / delta_t);
    W_vectors.reserve(tf / delta_t);
    L_vectors.reserve(tf / delta_t);
    T_vectors.reserve(tf / delta_t);
    rcm_data.reserve(tf / delta_t);
    L_lab_data.reserve(tf / delta_t);
    L_top_data.reserve(tf / delta_t);
    w_lab_data.reserve(tf / delta_t);
    w_top_data.reserve(tf / delta_t);
    theta_phi_phidot_data.reserve(tf / delta_t);
}

void evolve()
{
    double L[3]; 
    if (t <= tf && run_animation)
    {
        //RUNGE-KUTTA LOOP
        auto new_time = std::chrono::steady_clock::now();
        double frame_time = std::chrono::duration_cast<std::chrono::microseconds>(new_time - previous_time).count() * 1.e-6;
        previous_time = new_time;
        //frame_time e' il tempo tra un frame e l'altro. l'algoritmo avanza del numero
        //intero di timestep più vicino a frame_time (in difetto)

        if (!pause_animation) accumulator += frame_time;
        //accumulator accumula il tempo che la simulazione deve calcolare.
        //a ogni step viene decrementato di dt. se e' rimasto un resto
        //questo aumenta fino a poter costutire un intero delta_t
        //quindi ogni tanto viene calcolato un timestep in più tra un frame
        //e il successivo, ma e' r_top_quatstanza trascurabile.
        //in pratica ogni volta viene mandato fuori un fotogramma,
        //poi si guarda quanto tempo e' passato, si evolve fino al nuovo tempo e si manda fuori il nuovo fotogramma



        while (accumulator >= delta_t && t<=tf)
        {
            FORCE[0] = values[7];
            FORCE[1] = values[8];
            FORCE[2] = values[9];
            if (enable_gravity->value()) FORCE[2] -= mass * 9.81;
            //FORCE[0] += cos(2*t);
            //FORCE[1] += sin(t);
            

            //S = lab frame, S' = spinning top frame
            //GET FORCE COMPONENTS IN S'*************************
            Quaternion FORCE_quat(FORCE);
            Quaternion force_top_quat = inverse(q) * FORCE_quat * q;
            for (int i = 0; i < 3; i++) { FORCE_top[i] = force_top_quat[i + 1]; }
            //****************************************************


            //CALCULATE TORQUE COMPONENTS IN S'******************
            TORQUE[0] = -Rcm_top[2] * FORCE_top[1];
            TORQUE[1] = Rcm_top[2] * FORCE_top[0];
            TORQUE[2] = 0; 
            //***************************************************


            //Obtain w(t+dt) from da w(t) e tau(t) IN S'****************************************
            //FORMULA PRINCIPALE DELL'ALGORITMO.
            for (int i = 0; i < 3; ++i) eval0[i] = F[i](w0);
            for (int i = 0; i < 3; ++i) field1[i] = w0[i] + delta_t / 2 * eval0[i];
            for (int i = 0; i < 3; ++i) eval1[i] = F[i](field1);
            for (int i = 0; i < 3; ++i) field2[i] = w0[i] + delta_t / 2 * eval1[i]; 
            for (int i = 0; i < 3; ++i) eval2[i] = F[i](field2);
            for (int i = 0; i < 3; ++i) field3[i] = w0[i] + delta_t * eval2[i]; 
            for (int i = 0; i < 3; ++i) eval3[i] = F[i](field3); 

            for (int i = 0; i < 3; ++i) w0[i] += delta_t / 6 * (eval0[i] + 2 * eval1[i] + 2 * eval2[i] + eval3[i]);
            //*******************************************************************************************            
            //end of Runge-Kutta--------------------


            //Get wt components in S*********************************************
            Quaternion w_top_quat(w0);
            Quaternion w_quat = q * w_top_quat * inverse(q);
            for (int i = 0; i < 3; i++) { wt[i] = (w_quat[i + 1] * delta_t); }
            //*******************************************************************


            //Calculate RCM IN S***************************************************
            Quaternion Rcm_quat = q * Rcm_top_quat * inverse(q);
            for (int i = 0; i < 3; i++) { Rcm[i] = Rcm_quat[i + 1]; }
            //*********************************************************************


            //TRAIL DATA GENERATION***********************************************-
            double trace_lenght = opengl_top_scale_factor * (height_cone1 + height_cone2 + 1.01 * cyl_height) / vector_modulus(Rcm);
            trail_points.push_back({ Rcm[0] * trace_lenght, Rcm[1] * trace_lenght, Rcm[2] * trace_lenght });
            //*********************************************************************


            //UPDATE q, ADDING NEW INFINITESIMAL ROTATION GIVEN BY wt****************
            q = quaternionFromVector(wt) * q;  //normalize_this(q);
            //***********************************************************************

            L[0] = w0[0] * I1, L[1] = w0[1] * I2, L[2] = w0[2] * I3;
            Quaternion L_top_quat(L); double L_lab[3];
            Quaternion L_lab_quat = q * L_top_quat * inverse(q);
            for (int i = 0; i < 3; i++) { L_lab[i] = L_lab_quat[i + 1]; }
            double theta = acos(Rcm[2] / Rcm_top[2]); //arccos(z/l_cm)
            double phi = atan2(Rcm[1], Rcm[0]);
            double phidot = (L_lab[2] - L[2] * cos(theta)) / (I1 * sin(theta) * sin(theta));
            double E =(I1*w0[0]*w0[0]+I2*w0[1]*w0[1]+I3*w0[2]*w0[2])/2. + mass*9.81*Rcm[2];

            rcm_data.push_back({ Rcm [0], Rcm[1], Rcm[2] });
            L_lab_data.push_back({ L_lab[0], L_lab[1], L_lab[2] });
            L_top_data.push_back({ L[0], L[1], L[2] });
            w_lab_data.push_back({ w_quat[1], w_quat[2], w_quat[3] });
            w_top_data.push_back({ w0[0], w0[1], w0[2] });
            Energy_data.push_back(E);
            theta_phi_phidot_data.push_back({ theta, phi, phidot });

            if (enable_file_output->value())
            {
                L_stream << t << '\t'
                    << L[0] << '\t' << L[1] << '\t' << L[2] << '\t'
                    << L_lab[0] << '\t' << L_lab[1] << '\t' << L_lab[2] << '\n';

                w_stream << t << '\t'
                    << w0[0] << '\t' << w0[1] << '\t' << w0[2] << '\t'
                    << w_quat[1] << '\t' << w_quat[2] << '\t' << w_quat[3] << '\n';


                Rcm_stream << t << '\t'
                    << Rcm[0] << '\t' << Rcm[1] << '\t' << Rcm[2] << '\n';

                Energy_stream << t << '\t' << E << '\n';

                theta_phi_phidot_stream << t << '\t'
                    << theta << '\t' << phi << '\t' << phidot << '\n';
            }

            accumulator -= delta_t;
            t += delta_t;
        }

        double rotation_angle = 2 * atan2(sqrt(q[1] * q[1] + q[2] * q[2] + q[3] * q[3]), q[0]);
        xRotated = 180. / M_PI * rotation_angle,
        a1 = q[1],
        a2 = q[2],
        a3 = q[3];
        rotation_data.push_back({ rotation_angle, q[1], q[2], q[3] });
        double rotation_for_W[3];
        cross_product(Rcm_top, w0, rotation_for_W);
        //double W_modulus = vector_modulus(rotation_for_W);
        W_vectors.push_back({ 180. / M_PI * acos(dot_product(Rcm_top, w0) / (vector_modulus(Rcm_top) * vector_modulus(w0))), rotation_for_W[0], rotation_for_W[1], rotation_for_W[2] });
        double rotation_for_L[3];
        cross_product(Rcm_top, L, rotation_for_L);
        //double L_modulus = vector_modulus(rotation_for_L);
        L_vectors.push_back({ 180. / M_PI * acos(dot_product(Rcm_top, L) / (vector_modulus(Rcm_top) * vector_modulus(L))), rotation_for_L[0], rotation_for_L[1], rotation_for_L[2] });
        double rotation_for_T[3];
        cross_product(Rcm_top, TORQUE, rotation_for_T);
        double T_modulus = vector_modulus(rotation_for_T);
        T_vectors.push_back({ 180. / M_PI * acos(dot_product(Rcm_top, TORQUE) / (vector_modulus(Rcm_top) * vector_modulus(TORQUE))), rotation_for_T[0]/ T_modulus, rotation_for_T[1]/ T_modulus, rotation_for_T[2]/ T_modulus });

    }
    else
    {
        std::cout << "\nSimulazione terminata.\n";

        if (enable_file_output->value())
        {
            L_stream.clear();
            w_stream.clear();
            Rcm_stream.clear();
            Energy_stream.clear();
            theta_phi_phidot_stream.clear();

            L_stream.close();
            w_stream.close();
            Rcm_stream.close();
            Energy_stream.close();
            theta_phi_phidot_stream.close();

        }

        run_animation = false; hide_vectors = true;
        activateButtons();
    }

}

void startAlgorithm()
{
    t = 0; accumulator = 0;
    run_animation = true; hide_vectors = false;
    previous_time = std::chrono::steady_clock::now();
}

void animationLoop()
{
#ifdef SLEEP
    
    std::this_thread::sleep_for(std::chrono::microseconds((int)(1'000'000 / values[14]) - (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - idle_prev_time)).count()));
    idle_prev_time = std::chrono::steady_clock::now();

    if (run_animation) evolve();
    else
    {
        xRotated = 180. / M_PI * (values[3]),
        a1 = 1, a2 = 0, a3 = 0;
    }

#else
    if ((std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - idle_prev_time)).count() * 1.e-6 >= 1. / values[14])
    {
        idle_prev_time = std::chrono::steady_clock::now();
        
        if (!run_animation)
        {
            xRotated = 180. / M_PI * (values[3]);
            a1 = 1, a2 = 0, a3 = 0;
        }
        else evolve();
        goDraw = true;
    }
#endif
         
}

int main (int argc, char **argv) {

        CreateMyWindow(argc, argv);
}



/*LISTA INVERSIONI DEL WINDING:
(predefinito di opengl: CCW)

inversione per sottopavimento.
inversione per tronco di cono interna alla funzione (glQuadStrips sono CW)
inversione per il cono sotto perche è scalato -1 su z.


LISTA INVERSIONE DI CULL:
cull face front per trottola riflessa sul pavimento
*/

/*


POSSIBILI GRAFICI:
theta(t)
phi(t)
phidot(t)
phi(phidot)
phidot&theta (t)
w1(t) & w2(t)
w1(w2)
wx & wy & wz (t)
Lx & Ly & Lz (t)
Energy(t)
*/
