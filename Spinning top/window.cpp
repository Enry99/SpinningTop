//#define AXIS_BOXES //currently not working properly. when entering values sometimes the spinning top misteriously disappears.
//#define CHANGE_SYSTEM_TIMER_RESOLUTION
#define SLEEP
#ifndef SLEEP
extern bool goDraw;
#endif
#if defined(_WIN32) || defined(WIN32)
#define _USE_MATH_DEFINES
#endif


#include "slider_input.h"
#include "AxisRangeInput.h"
#include <FL/gl.h>
#include <FL/glu.h>
#include <FL/glut.H>
#include <FL/x.H>
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Toggle_Button.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Float_Input.H>
#include <FL/Fl_Choice.H>
#include <vector>
#include <iostream>
#include <fstream>
#include <array>
#include <math.h>
#include <stdio.h>
#include <string.h>


//extern functions/variables (declared in main.cpp)
extern void setInitialConditions();
extern void startAlgorithm();
extern void displaySpinningTop();
extern void animationLoop();
extern void keyboardFunction(unsigned char, int, int);
extern void setAxisRange();
extern void drawGraph();
extern GLfloat amb_light[4];
extern GLfloat diffuse[4];
extern GLfloat specular[4];
extern std::vector<std::array<double, 4>> rotation_data;
extern std::vector<std::array<double, 3>> trail_points;
extern std::vector<std::array<double, 4>> W_vectors;
extern std::vector<std::array<double, 4>> L_vectors;
extern std::vector<std::array<double, 4>> T_vectors;
extern std::vector<std::array<double, 3>> rcm_data;
extern std::vector<std::array<double, 3>> L_lab_data;
extern std::vector<std::array<double, 3>> L_top_data;
extern std::vector<std::array<double, 3>> w_lab_data;
extern std::vector<std::array<double, 3>> w_top_data;
extern std::vector<std::array<double, 3>> theta_phi_phidot_data;
extern std::vector<double> Energy_data;
extern std::ofstream L_stream;
extern std::ofstream w_stream;
extern std::ofstream Energy_stream;
extern std::ofstream theta_phi_phidot_stream;
extern std::ofstream Rcm_stream;
extern double values[15];
extern bool run_animation;
extern bool pause_animation;
extern bool PERSISTENT_TRAIL;
extern float xmin_graph;
extern float xmax_graph;
extern float ymin_graph;
extern float ymax_graph;
extern int last_index;
extern int graphID;
extern const char* menu_labels[17];


//Function declarations
void idleSpinningTop(void*);
void idleGraph(void*);
void mouseFunction(int, int, int, int);
void mouseWheelFunction(int);
void mouseMotionCallback(int, int);
void graphmouseFunction(int, int, int, int);
void graphmouseWheelFunction(int);
void graphmouseMotionCallback(int, int);
void createMenu(int);
void activateButtons();
void deactivateButtons();
void start_callback(Fl_Widget*);
void pause_callback(Fl_Widget*);
void stop_callback(Fl_Widget*);
void main_window_cb(Fl_Widget*, void*);
int CreateMyWindow(int, char**);


//Variables declarations
double camera_distance = 10.;
double cameraPos[3] { camera_distance * sin(M_PI / 3.) * cos(M_PI / 4.), camera_distance * sin(M_PI / 3.) * sin(M_PI / 4.), camera_distance * cos(M_PI / 3.) };
double cameraUp[3] { 0.0,0.0,1.0 };
double cameraFront[3] { -3.0,-3.0,-3.0 };
double theta = 60.0f;
double phi = 45.0f;
double zoom_scale_factor = 1;
int lastX = 0;
int lastY = 0;
int graphlastX = 0;
int graphlastY = 0;
bool mousePressed = false;
double FPS_display_width, FPS_display_height;

//sliders default values
constexpr int FPS_default = 75;
constexpr double
initial_angular_velocity_def[3] = { 0., 0., 80. },
additional_force_def[3] = { 0., 0., 0. },
mass_def = 0.1,
time_increment_def = 0.001,
total_time_def = 60.,
initial_angle_def = 0.3,
low_height_def = 0.025,
high_height_def = 0.02,
cone_radius_def = 0.04,
cylinder_height_def = 0.015;
constexpr bool enable_gravity_def = true;
constexpr bool enable_file_output_def = false;
constexpr bool persistent_trail_def = false;

//right-click menu labels
char plus_label[] = "Increase brightness [+]";
char minus_label[] = "Decrease brightness [-]";
char t_label[] = "Toggle Trail [t]";
char r_label[] = "Toggle Rainbow Trail [r]";
char v_label[] = "Toggle Vectors [v]";
char a_label[] = "Toggle axis [a]";
char f_label[] = "Toggle floor [f]";
char s_label[] = "Toggle Sun [s]";
char x_label[] = "Toggle FPS [space]";
char h_label[] = "Hide spinning top [h]";


//interface elements
class MyGlutWindow;
class MyGlutWindow2;
Fl_Window * main_window;
MyGlutWindow* spinning_top_window;
MyGlutWindow2* graph_window;
Fl_Button * start_button, * stop_button, * pause_button, * reset_button;
Fl_Text_Display * vector_legend_box;
SliderInput 
* mass_slider,
* initial_angular_velocity[3], 
* additional_force[3],
* time_increment, 
* total_time, 
* initial_angle, 
* low_height, 
* high_height, 
* cone_radius, 
* cylinder_height, 
* FPS_slider;
Fl_Check_Button * enable_gravity;
Fl_Check_Button* enable_file_output;
Fl_Check_Button* persistent_trail;
Fl_Check_Button* autorange_button;
Fl_Check_Button* x_zoom_only;
Fl_Text_Buffer * tbuff;
Fl_Text_Buffer * sbuff;
Fl_Choice* menu;
AxisRangeInput* axis_boxes;

//end of declarations_____________________________


//Spinning top window class
class MyGlutWindow : public Fl_Glut_Window {

    void init()
    {
        glClearColor(0.0, 0.0, 0.0, 0.0);

        //right-click menu on spinnning top's viewport
        glutCreateMenu(createMenu);
        glutAddMenuEntry(plus_label, '+');
        glutAddMenuEntry(minus_label, '-');      
        glutAddMenuEntry(t_label, 't');
        glutAddMenuEntry(r_label, 'r');
        glutAddMenuEntry(v_label, 'v');
        glutAddMenuEntry(a_label, 'a');
        glutAddMenuEntry(f_label, 'f');
        glutAddMenuEntry(s_label, 's');
        glutAddMenuEntry(x_label, 'x');
        glutAddMenuEntry(h_label, 'h');
        glutAttachMenu(GLUT_RIGHT_BUTTON);

        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
        glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.1);
        glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.05);
        //glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.002);
        glShadeModel(GL_SMOOTH);
        glEnable(GL_LINE_SMOOTH);
        glDepthFunc(GL_LEQUAL);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_CULL_FACE);
        glEnable(GL_TEXTURE_2D);

        glPolygonOffset(-2.0, -1.0); //important: avoids depth buffer aliasing for shadows   
    }


    void FixViewport(int W, int H) 
    {
        if (W == 0 || H == 0) return;

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(45.0, (GLdouble)W / (GLdouble)H, 0.01, 100.);
        glViewport(0, 0, pixel_w(), pixel_h());  //Use the whole window for rendering

        FPS_display_width = (double)W / Fl::w() / (0.5 * 0.9), //size is fixed (does not vary if resized)
        FPS_display_height = (double)H / Fl::h() / (0.8 * 0.95);

        glMatrixMode(GL_MODELVIEW);
    }


    void draw() // DRAW METHOD
    {
        this->make_current();
        static bool first_time = true;
        if (first_time) { valid(1); init();  FixViewport(w(), h()); first_time = false; }
       
        displaySpinningTop();
    }

    void resize(int X, int Y, int W, int H) 
    {
        this->make_current();
        Fl_Gl_Window::resize(X, Y, W, H);
        FixViewport(W, H);
        redraw();
    }

public:


    // OPENGL WINDOW CONSTRUCTOR
    MyGlutWindow(int X, int Y, int W, int H, const char* L = 0, bool use_mouse_keyboard = false) : Fl_Glut_Window(X, Y, W, H, L)
    {
        mode(FL_RGB | FL_DEPTH | FL_STENCIL | FL_MULTISAMPLE);
        Fl::add_idle(idleSpinningTop, this);
        if (use_mouse_keyboard) {
            this->mouse = mouseFunction;
            this->motion = mouseMotionCallback;
            this->keyboard = keyboardFunction;
        }
        end();
    }

};


//graph window class
class MyGlutWindow2 : public Fl_Glut_Window {

    void init()
    {
        glClearColor(0.0, 0.0, 0.0, 0.0);
        glEnable(GL_DEPTH_TEST);
    }

    void FixViewport(int W, int H) {
        if (W == 0 || H == 0) return;

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        setAxisRange();
        glViewport(0, 0, pixel_w(), pixel_h());  //Use the whole window for rendering
        glMatrixMode(GL_MODELVIEW);

    }
    // DRAW METHOD
    void draw()
    {
        this->make_current();
        static bool first_time = true;
        if (first_time) { valid(1); init();  FixViewport(w(), h()); first_time = false; }

        drawGraph();
    }

    void resize(int X, int Y, int W, int H)
    {
        this->make_current();
        Fl_Gl_Window::resize(X, Y, W, H);
        FixViewport(W, H);
        redraw();
    }

public:


    // OPENGL WINDOW CONSTRUCTOR
    MyGlutWindow2(int X, int Y, int W, int H, const char* L = 0) : Fl_Glut_Window(X, Y, W, H, L)
    {
        
        mode(FL_RGB | FL_DEPTH);
        Fl::add_idle(idleGraph, this);
        this->mouse = graphmouseFunction;
        this->motion = graphmouseMotionCallback;
        end();
    }

};


void idleSpinningTop(void*)
{
    animationLoop(); //updates data
#ifdef SLEEP  
    spinning_top_window->redraw();
#else
    if (goDraw) {
        spinning_top_window->redraw();
        goDraw = false;
    }
#endif
}


void idleGraph(void*)
{    
    graph_window->redraw();
}


void mouseFunction(int button, int state, int x, int y)
{
    lastX = x;
    lastY = y;
    mousePressed = state == GLUT_DOWN;
    if (Fl::event() == FL_MOUSEWHEEL || Fl::event() == FL_ZOOM_GESTURE) mouseWheelFunction(Fl::event_dy());
}


void mouseWheelFunction(int wheel_motion)
{
    if (camera_distance + zoom_scale_factor * wheel_motion < 0.001 || camera_distance + zoom_scale_factor * wheel_motion> 50.) return;
    camera_distance += zoom_scale_factor * wheel_motion;
    cameraPos[0] = camera_distance * sin(M_PI / 180. * theta) * cos(M_PI / 180. * phi);
    cameraPos[1] = camera_distance * sin(M_PI / 180. * theta) * sin(M_PI / 180. * phi);
    cameraPos[2] = camera_distance * cos(M_PI / 180. * theta);
}


void mouseMotionCallback(int xpos, int ypos)
{
    if (mousePressed) //mousePressed is updated by mouseFunction
    {
        float dx = xpos - lastX;  //x-coordinates go from left to right
        float dy = lastY - ypos; // reversed since y-coordinates go from bottom to top
        lastX = xpos;
        lastY = ypos;

        float xSensitivity = 2. * 360. / Fl::w(); //sensitivities are normalized to window size to produce the same result 
        float ySensitivity = 2. * 180. / Fl::h(); //regardless of pixel density inside the window 
        dx *= xSensitivity;
        dy *= ySensitivity;

        theta += dy; //if mouse motion is upward, camera should go downward, thus theta is incremented     
        phi -= dx; //if mouse motion is to the right, camera should go to the left, thus phi is decremented

        // make sure that when theta is out of bounds, screen doesn't get flipped
        if (theta > 179.9f) theta = 179.9f;
        if (theta < 0.1f) theta = 0.1;

        cameraPos[0] = camera_distance * sin(M_PI / 180. * theta) * cos(M_PI / 180. * phi);
        cameraPos[1] = camera_distance * sin(M_PI / 180. * theta) * sin(M_PI / 180. * phi);
        cameraPos[2] = camera_distance * cos(M_PI / 180. * theta);
    }
}


void graphmouseFunction(int button, int state, int x, int y)
{
    graphlastX = x;
    graphlastY = y;
    mousePressed = state == GLUT_DOWN;
    if (Fl::event() == FL_MOUSEWHEEL || Fl::event() == FL_ZOOM_GESTURE) graphmouseWheelFunction(Fl::event_dy());
}


void graphmouseWheelFunction(int wheel_motion)
{
    float zoom_scale = 0.03;
    
    float delta_x = xmax_graph - xmin_graph;
    float delta_y = ymax_graph - ymin_graph;
    
    xmin_graph -= zoom_scale * delta_x * wheel_motion;
    xmax_graph += zoom_scale * delta_x * wheel_motion;

    if (!x_zoom_only->value())
    {
        ymin_graph -= zoom_scale * delta_y * wheel_motion;
        ymax_graph += zoom_scale * delta_y * wheel_motion;
    }

    setAxisRange();
}


void graphmouseMotionCallback(int xpos, int ypos)
{
    if (mousePressed) //mousePressed is updated by mouseFunction
    {
        float dx = graphlastX - xpos; 
        float dy = ypos - graphlastY;
        graphlastX = xpos;
        graphlastY = ypos;

        float motion_scale = 1.;
        float delta_x = xmax_graph - xmin_graph;
        float delta_y = ymax_graph - ymin_graph;
        float xSensitivity = motion_scale * delta_x / graph_window->pixel_w();
        float ySensitivity = motion_scale * delta_y / graph_window->pixel_h();
        dx *= xSensitivity;
        dy *= ySensitivity;

        xmin_graph += dx, xmax_graph += dx;
        ymin_graph += dy, ymax_graph += dy;
        setAxisRange();
    }
}


void graphChoiceMenuCallback(Fl_Widget* f)
{
    graphID = ((Fl_Choice*)f)->value();
    last_index = 0; //reset for autorange
}

#ifdef AXIS_BOXES
void autorangeButtonCallback(Fl_Widget* f)
{
    if(((Fl_Check_Button*)f)->value()) axis_boxes->deactivate();
    else axis_boxes->activate();
    axis_boxes->setvalues(xmin_graph, xmax_graph, ymin_graph, ymax_graph);
}
#endif


void createMenu(int key) { 
    keyboardFunction((unsigned char)key, 0, 0); 
}


void activateButtons()
{
    mass_slider->activate();
    for (int i = 0; i < 3; ++i) initial_angular_velocity[i]->activate();
    
    time_increment->activate();
    total_time->activate();
    initial_angle->activate();
    low_height->activate();
    high_height->activate();
    cone_radius->activate();
    cylinder_height->activate();
    start_button->activate();
    reset_button->activate();
    enable_file_output->activate();
}


void deactivateButtons()
{
    mass_slider->deactivate();
    for (int i = 0; i < 3; ++i) initial_angular_velocity[i]->deactivate();
    time_increment->deactivate();
    total_time->deactivate();
    initial_angle->deactivate();
    low_height->deactivate();
    high_height->deactivate();
    cone_radius->deactivate();
    cylinder_height->deactivate();
    start_button->deactivate();
    reset_button->deactivate();
    enable_file_output->deactivate();
}


void start_callback(Fl_Widget* f) {

    rotation_data.clear(), trail_points.clear(), W_vectors.clear(), L_vectors.clear(), T_vectors.clear();
    rcm_data.clear(), L_lab_data.clear(), L_top_data.clear(), w_lab_data.clear(), w_top_data.clear(), theta_phi_phidot_data.clear(), Energy_data.clear();
    last_index = 0;

    pause_animation = false, pause_button->clear();

    deactivateButtons();

    setInitialConditions();

    startAlgorithm();

}


void pause_callback(Fl_Widget* f) { 
    pause_animation= !pause_animation; 
}


void stop_callback(Fl_Widget* f) {

    run_animation = false, pause_animation = false, pause_button->clear();

    rotation_data.clear(), trail_points.clear(), W_vectors.clear(), L_vectors.clear(), T_vectors.clear();

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
    activateButtons();
}


void reset_callback(Fl_Widget* f)
{
   
    for (int i = 0; i < 3; ++i)
    {
        initial_angular_velocity[i]->value(initial_angular_velocity_def[i]);
        additional_force[i]->value(additional_force_def[i]);
    }
    mass_slider->value(mass_def),
    time_increment->value(time_increment_def),
    total_time->value(total_time_def),
    initial_angle->value(initial_angle_def),
    low_height->value(low_height_def),
    high_height->value(high_height_def),
    cone_radius->value(cone_radius_def),
    cylinder_height->value(cylinder_height_def),
    FPS_slider->value(FPS_default),
    enable_gravity->value(enable_gravity_def);

}


void main_window_cb(Fl_Widget* widget, void*)
{
#if defined SLEEP && (defined(_WIN32) || defined(WIN32)) && defined CHANGE_SYSTEM_TIMER_RESOLUTION
    timeEndPeriod(1);
#endif

    return exit(EXIT_SUCCESS);
}


int CreateMyWindow(int argc, char** argv) {
    
#if defined SLEEP && (defined(_WIN32) || defined(WIN32)) && defined CHANGE_SYSTEM_TIMER_RESOLUTION
    timeBeginPeriod(1); //allows better control over FPS in Windows. 
    //Use with caution, since it changes the minimum resolution of periodic timers of the whole system and thus increases cpu usage.
#endif

    Fl::scheme("gtk+"); //style, comment to get traditional look

    Fl::use_high_res_GL(1); //for Apple retina display
    
    //SETTING CONTAINER WINDOW_________________________________________________________________________________   
    int w_ext = 0.9 * Fl::w(); //container window sizes
    int h_ext = 0.85 * Fl::h(); 
    double x0 = w_ext * 0.82; //x position of sliders (pixel count from main window's left edge)
    double y0 = 0.04 * h_ext;
    double delta_h = h_ext * 0.054; //height step between sliders' upper-left corner 
    double slider_width = w_ext * 0.16; 
    double slider_height = h_ext * 0.027; // 0.042;
    double button_width = 0.9 * slider_width;
    double button_height = 2.5*h_ext*0.042;

    main_window = new Fl_Window((Fl::w()- w_ext)/2, (Fl::h()-h_ext)/2, w_ext, h_ext, "Lagrange spinning top");
    main_window->resizable(main_window);
    main_window->callback(main_window_cb);
    main_window->show(argc, argv);
    main_window->begin();
    

    //SETTING SUBWINDOWS AND SLIDERS___________________________________________________________________________

    //subwindows
    spinning_top_window = new MyGlutWindow(0.025 * w_ext, 0.025 * h_ext, 0.45 * w_ext, 0.95 * h_ext, "Top_window", true);
    graph_window = new MyGlutWindow2(x0 - 1.88 * slider_width, 0.025 * h_ext, 0.9 * slider_width+ x0 - 0.94 * slider_width - (x0 - 1.88 * slider_width) , 0.5 * main_window->h());


    //sliders
    mass_slider = new SliderInput(x0, y0, slider_width, slider_height, "mass (kg)", &values[0]);
    time_increment = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "dt (s)", &values[1]);
    total_time = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "total time (s)", &values[2]);
    initial_angle = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "theta_0 (rad)", &values[3]);
    initial_angular_velocity[0] = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "omega1 (rad/s)", &values[4]);
    initial_angular_velocity[1] = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "omega2 (rad/s)", &values[5]);
    initial_angular_velocity[2] = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "omega3 (rad/s)", &values[6]);
    additional_force[0] = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "Fx_add (N)", &values[7]);
    additional_force[1] = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "Fy_add (N)", &values[8]);
    additional_force[2] = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "Fz_add (N)", &values[9]);
    low_height = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "h_1 (m)", &values[10]);
    high_height = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "h_2 (m)", &values[11]);
    cone_radius = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "cone radius (m)", &values[12]);
    cylinder_height = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "cylinder height (m)", &values[13]);
    FPS_slider = new SliderInput(x0, y0 += delta_h, slider_width, slider_height, "FPS max", &values[14]);


    //buttons
    start_button = new Fl_Button(x0, y0 += 1.1*delta_h, button_width, button_height, "@#>\tStart");
    start_button->color(FL_GREEN);
    pause_button = new Fl_Toggle_Button(x0 - 0.94 * slider_width, y0, button_width, button_height, "@#||\tPause/Unpause");
    pause_button->color(FL_YELLOW);
    stop_button = new Fl_Button(x0 - 1.88 * slider_width, y0, button_width, button_height, "@#square\tStop");
    stop_button->color(FL_RED);
    reset_button = new Fl_Button(x0 - 1.88 * slider_width, y0 - 0.75* button_height, button_width, 0.6*button_height, "Reset sliders");
    reset_button->color(FL_CYAN);
    vector_legend_box = new Fl_Text_Display(x0 - 1.88 * slider_width, 0.65 * h_ext, button_width, 0.9 * button_height, "Vectors legend");
    enable_gravity = new Fl_Check_Button(x0 - 0.7*slider_width, 0.78 * h_ext, 0.4*slider_width, slider_height, "Enable gravity");
    enable_file_output = new Fl_Check_Button(x0 - 0.7*slider_width, 0.65 * h_ext, 0.5*slider_width, slider_height, "Output data on file");
    persistent_trail = new Fl_Check_Button(x0 - 0.7 * slider_width, 0.81 * h_ext, 0.4 * slider_width, slider_height, "Persistent trail");
    autorange_button = new Fl_Check_Button(x0 - 0.7 * slider_width, 0.028 * h_ext + 0.5 * main_window->h(), 0.5 * slider_width, slider_height, "Auto range");
    x_zoom_only = new Fl_Check_Button(x0 - 0.7 * slider_width, 0.028 * h_ext + 0.5 * main_window->h()+slider_height, 0.5 * slider_width, slider_height, "X axis only zoom");
    menu = new Fl_Choice(x0 - 1.88 * slider_width, 0.028 * h_ext + 0.5 * main_window->h(), 0.9*slider_width, slider_height);
#ifdef AXIS_BOXES
    axis_boxes = new AxisRangeInput(x0 - 1.88 * slider_width, 0.045 * h_ext + 0.5 * main_window->h() + slider_height, 1.1*slider_width, slider_height,
        &xmin_graph, &xmax_graph, &ymin_graph, &ymax_graph, "xmin, xmax, ymin, ymax");
    axis_boxes->deactivate();
#endif


    for (auto i : menu_labels) menu->add(i);
    menu->callback(graphChoiceMenuCallback);
#ifdef AXIS_BOXES
    autorange_button->callback(autorangeButtonCallback);
#endif

    //setting default values and bounds
    for (int i = 0; i < 3; ++i) 
    {
        initial_angular_velocity[i]->bounds(-500., 500.);
        initial_angular_velocity[i]->value(initial_angular_velocity_def[i]);

        additional_force[i]->bounds(-10., 10.);
        additional_force[i]->value(additional_force_def[i]);
    }

    mass_slider->bounds(0.001, 5.);
    mass_slider->value(mass_def);

    time_increment->bounds(1.e-4, 0.01);
    time_increment->step(0.00001);
    time_increment->value(time_increment_def);

    total_time->bounds(1., 600.);
    total_time->value(total_time_def);

    initial_angle->bounds(0, 2 * M_PI);
    initial_angle->value(initial_angle_def);

    low_height->bounds(0.003, .2);
    low_height->value(low_height_def);

    high_height->bounds(0.003, .2);
    high_height->value(high_height_def);

    cone_radius->bounds(0.001, .2);
    cone_radius->value(cone_radius_def);

    cylinder_height->bounds(0, .1);
    cylinder_height->value(cylinder_height_def);

    FPS_slider->bounds(20, 360);
    FPS_slider->step(1);
    FPS_slider->value(FPS_default);

    enable_gravity->value(enable_gravity_def);
    enable_file_output->value(enable_file_output_def);
    persistent_trail->value(persistent_trail_def);

    start_button->callback(start_callback);
    pause_button->callback(pause_callback);
    stop_button->callback(stop_callback);
    reset_button->callback(reset_callback);
    autorange_button->value(1);
    x_zoom_only->value(0);

  
    //Vectors legend***************************************************
    Fl_Text_Display::Style_Table_Entry sTable[] = {
      // FONT COLOR      FONT FACE   FONT SIZE
      // --------------- ----------- --------------
      {  FL_RED,         FL_HELVETICA_BOLD, 17 }, // A - Red
      {  0xffdc1e00,      FL_HELVETICA_BOLD, 17 }, // B - Yellow
      {  FL_DARK_GREEN,  FL_HELVETICA_BOLD, 17 }, // C - Green
      {  FL_BLUE,        FL_HELVETICA_BOLD, 17 }, // D - Blue
    };

    tbuff = new Fl_Text_Buffer;      // text buffer
    sbuff = new Fl_Text_Buffer;      // style buffer
    vector_legend_box->buffer(tbuff);
    int sTable_size = sizeof(sTable) / sizeof(sTable[0]);
    vector_legend_box->highlight_data(sbuff, sTable, sTable_size, 'A', 0, 0);

    // Text
    tbuff->text("Torque\nAngular velocity\nAngular momentum\n");
    // Style for text
    sbuff->text("BBBBBB\nCCCCCCCCCCCCCCCC\nDDDDDDDDDDDDDDDD\n");
    //******************************************************************
     
    main_window->end();
    
    main_window->show();
    spinning_top_window->show();
    graph_window->show();

    return (Fl::run());
}