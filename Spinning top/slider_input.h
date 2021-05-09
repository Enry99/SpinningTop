#ifndef SLIDER_INPUT_H
#define SLIDER_INPUT_H
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Float_Input.H>
#include <FL/Fl_Slider.H>
#include <stdio.h>

// sliderinput -- simple example of tying an fltk slider and input widget together
// 1.00 erco 10/17/04
//modified by us in 2020, added option to handle an external variable, updating its value with slider's value.

class SliderInput : public Fl_Group {
    Fl_Float_Input *input;
    Fl_Slider    *slider;
    double *handled_value;

    // CALLBACK HANDLERS
    //    These 'attach' the input and slider's values together.
    //
    void Slider_CB2() 
    {
        static int recurse = 0;
        if ( recurse ) return;
        else 
        {
            recurse = 1;
            char s[80];
#if defined(WIN32)
            sprintf_s(s, "%lf", slider->value() );
#else
            sprintf(s, "%lf", slider->value());
#endif
            // fprintf(stderr, "SPRINTF(%d) -> '%s'\n", (int)(slider->value()+.5), s);
            input->value(s);          // pass slider's value to input
            if (handled_value) *handled_value = slider->value();
            recurse = 0;
        }
    }

    static void Slider_CB(Fl_Widget *w, void *data) {
        ((SliderInput*)data)->Slider_CB2();
    }

    void Input_CB2() 
    {
        static int recurse = 0;
        if (recurse) return;
        else {
            recurse = 1;
            double val = 0;
#if defined(WIN32)
            if (sscanf_s(input->value(), "%lf", &val) != 1) val = 0;
#else
            if (sscanf(input->value(), "%lf", &val) != 1) val = 0;
#endif
            // fprintf(stderr, "SCANF('%s') -> %d\n", input->value(), val);
            double min = slider->minimum(); //if value is outside range, set it to the nearest bound of the slider's range
            double max = slider->maximum();
            if (val < min) val = min;
            else if (val > max) val = max;
            slider->value(val);         // pass input's value to slider

           char s[80];
#if defined(_WIN32) || defined(WIN32)
           sprintf_s(s, "%lf", val );
#else
           sprintf(s, "%lf", val);
#endif

            input->value(s);

            if (handled_value) *handled_value = val;

            recurse = 0;
        }
    }
    static void Input_CB(Fl_Widget *w, void *data) {
        ((SliderInput*)data)->Input_CB2();
    }

public:
    // CTOR
    SliderInput(double x, double y, double w, double h, const char *l=0, double* handled_val = nullptr) : Fl_Group(x,y,w,h,l), handled_value(handled_val) {
        int in_w = 90;
        input  = new Fl_Float_Input(x, y, in_w, h);
        input->callback(Input_CB, (void*)this);
        input->when(FL_WHEN_ENTER_KEY|FL_WHEN_NOT_CHANGED);

        slider = new Fl_Slider(x+in_w, y, w-in_w, h);
        slider->type(FL_HOR_NICE_SLIDER);
        slider->callback(Slider_CB, (void*)this);

        //bounds(1, 10);     // some usable default
        //value(5);          // some usable default
        end();             // close the group
    }

    // MINIMAL ACCESSORS
    double  value() const    { return slider->value() ; }
    void value(double val)   { slider->value(val); Slider_CB2(); }
    void minumum(double val) { slider->minimum(val); }
    double  minumum() const  { return slider->minimum(); }
    void maximum(double val) { slider->maximum(val); }
    double maximum() const  { return slider->maximum(); }
    void bounds(double low, double high) { slider->bounds(low, high); }
    void step(double epsilon){ slider->step(epsilon);}
    
};

#endif
