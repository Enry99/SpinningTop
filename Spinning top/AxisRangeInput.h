#ifndef AXIS_RANGE_INPUT
#define AXIS_RANGE_INPUT

#include <FL/Fl.H>
#include <FL/gl.h>
#include <FL/glu.h>
#include <FL/glut.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Float_Input.H>

class myFl_Float_Input : public Fl_Float_Input
{
    float* handled_val;
public:
    
  
    myFl_Float_Input(int X, int Y, int W, int H, float* value, const char* l = 0) : Fl_Float_Input(X, Y, W, H, l), handled_val(value)
    {

    }

    float handledvalue() { return *handled_val; }
    void handledvalue(float val) { *handled_val = val; }
};

class AxisRangeInput : public Fl_Group {
    myFl_Float_Input
        * xmin_box,
        * xmax_box,
        * ymin_box,
        * ymax_box;

    void box_callback(myFl_Float_Input* box)
    {
        static int recurse = 0;
        if (recurse) return;
        else {
            recurse = 1;
            double val = 0;
#if defined(WIN32)
            if (sscanf_s(box->value(), "%lf", &val) != 1) val = box->handledvalue();
#else
            if (sscanf(box->value(), "%lf", &val) != 1) val = box->handledvalue();
#endif
           
            char s[80];
#if defined(_WIN32) || defined(WIN32)
            sprintf_s(s, "%lf", val);
#else
            sprintf(s, "%lf", val);
#endif

            box->value(s);

            box->handledvalue(val);

            recurse = 0;

            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            gluOrtho2D(xmin_box->handledvalue(), xmax_box->handledvalue(), ymin_box->handledvalue(), ymax_box->handledvalue());
            glMatrixMode(GL_MODELVIEW);
           
        }
    }

    static void Input_CB(Fl_Widget* w, void* data) {
        ((AxisRangeInput*)data)->box_callback((myFl_Float_Input*)w);
    }

public:
    // CTOR
    AxisRangeInput(double x, double y, double w, double h, float* xmin, float* xmax, float* ymin, float* ymax, const char* l = 0) :
        Fl_Group(x, y, w, h, l)
    {
        xmin_box = new myFl_Float_Input(x, y, w / 4., h, xmin);
        xmax_box = new myFl_Float_Input(x+w/4., y, w / 4., h, xmax);
        ymin_box = new myFl_Float_Input(x+2*w/4, y, w / 4., h, ymin);
        ymax_box = new myFl_Float_Input(x+3*w/4, y, w / 4., h, ymax);

        xmin_box->callback(Input_CB, (void*)this);
        xmin_box->when(FL_WHEN_ENTER_KEY);
        xmax_box->callback(Input_CB, (void*)this);
        xmax_box->when(FL_WHEN_ENTER_KEY);
        ymin_box->callback(Input_CB, (void*)this);
        ymin_box->when(FL_WHEN_ENTER_KEY);
        ymax_box->callback(Input_CB, (void*)this);
        ymax_box->when(FL_WHEN_ENTER_KEY);
        end();             
    }

    void setvalues(float xmin, float xmax, float ymin, float ymax)
    {
        char s1[80];
        char s2[80];
        char s3[80];
        char s4[80];
#if defined(_WIN32) || defined(WIN32)
        sprintf_s(s1, "%lf", xmin);
        sprintf_s(s2, "%lf", xmax);
        sprintf_s(s3, "%lf", ymin);
        sprintf_s(s4, "%lf", ymax);
#else
        sprintf(s1, "%lf", xmin);
        sprintf(s2, "%lf", xmax);
        sprintf(s3, "%lf", ymin);
        sprintf(s4, "%lf", ymax);
#endif

        xmin_box->value(s1);
        xmax_box->value(s2);
        ymin_box->value(s3);
        ymax_box->value(s4);
    }
   
};

#endif