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
    double** handled_vals;
    int index;
public:
    
    myFl_Float_Input(int X, int Y, int W, int H, double** values, int index, const char* l = 0) : Fl_Float_Input(X, Y, W, H, l), handled_vals(values), index(index) { }

    float handledvalue() { return *handled_vals[index]; }
    void handledvalue(float val) { *handled_vals[index] = val; }

    void set_text_value(double val)
    {
        char s[80];
#if defined(_WIN32) || defined(WIN32)
        sprintf_s(s, "%lf", val);
#else
        sprintf(s, "%lf", val);
#endif
        this->value(s);
    }

    void update_handledvalue_from_text()
    {
        double val = 0;
#if defined(WIN32)
        if (sscanf_s(this->value(), "%lf", &val) != 1) val = *handled_vals[index];
#else
        if (sscanf(this->value(), "%lf", &val) != 1) val = *handled_vals[index];
#endif
        switch (index) //if xmin>=xmax, ymin>=ymax leaves previous value
        {
        case 0:
            if (val >= *handled_vals[1]) val = *handled_vals[index];
            break;
        case 1:
            if (val <= *handled_vals[0]) val = *handled_vals[index];
            break;
        case 2:
            if (val >= *handled_vals[3]) val = *handled_vals[index];
            break;
        case 3:
            if (val <= *handled_vals[2]) val = *handled_vals[index];
            break;
        }

        char s[80];
#if defined(_WIN32) || defined(WIN32)
        sprintf_s(s, "%lf", val);
#else
        sprintf(s, "%lf", val);
#endif

        this->value(s);

        *handled_vals[index] = val;
    }
};

class AxisRangeInput : public Fl_Group {
    myFl_Float_Input
        * xmin_box,
        * xmax_box,
        * ymin_box,
        * ymax_box;
    double* valueaddresses[4];

    void box_callback(myFl_Float_Input* box)
    {
        box->update_handledvalue_from_text();
    }

    static void Input_CB(Fl_Widget* w, void* data) {
        ((AxisRangeInput*)data)->box_callback((myFl_Float_Input*)w);
    }

public:
    // CTOR
    AxisRangeInput(double x, double y, double w, double h, double* xmin, double* xmax, double* ymin, double* ymax, const char* l = 0) :
        Fl_Group(x, y, w, h, l)
    {
        valueaddresses[0] = xmin;
        valueaddresses[1] = xmax;
        valueaddresses[2] = ymin;
        valueaddresses[3] = ymax;

        xmin_box = new myFl_Float_Input(x, y, w / 4., h, valueaddresses, 0);
        xmax_box = new myFl_Float_Input(x+w/4., y, w / 4., h, valueaddresses, 1);
        ymin_box = new myFl_Float_Input(x+2*w/4, y, w / 4., h, valueaddresses, 2);
        ymax_box = new myFl_Float_Input(x+3*w/4, y, w / 4., h, valueaddresses, 3);

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

    void set_text_values(float xmin, float xmax, float ymin, float ymax)
    {
        xmin_box->set_text_value(xmin);
        xmax_box->set_text_value(xmax);
        ymin_box->set_text_value(ymin);
        ymax_box->set_text_value(ymax);
    }

    void setvalues(float xmin, float xmax, float ymin, float ymax)
    {
        set_text_values(xmin, xmax, ymin, ymax);
        xmin_box->handledvalue(xmin);
        xmax_box->handledvalue(xmax);
        ymin_box->handledvalue(ymin);
        ymax_box->handledvalue(ymax);

    }

    ~AxisRangeInput()
    {
        delete xmin_box;
        delete xmax_box;
        delete ymin_box;
        delete ymax_box;
    }
  
};

#endif