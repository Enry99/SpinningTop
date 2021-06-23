#include <FL/gl.h>
#include <FL/glu.h>
#include <FL/glut.H>
#include <FL/x.H>
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Toggle_Button.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Choice.H>
#include <vector>
#include <array>
#include "AxisRangeInput.h"


extern Fl_Choice* menu;
extern Fl_Check_Button* autorange_button;
extern AxisRangeInput* axis_boxes;
extern std::vector<std::array<double, 3>> rcm_data;
extern std::vector<std::array<double, 3>> L_lab_data;
extern std::vector<std::array<double, 3>> L_top_data;
extern std::vector<std::array<double, 3>> w_lab_data;
extern std::vector<std::array<double, 3>> w_top_data;
extern std::vector<std::array<double, 3>> theta_phi_phidot_data;
extern std::vector<double> Energy_data;
extern double delta_t;

double xmin_graph = -1;
double xmax_graph = 1;
double ymin_graph = -1;
double ymax_graph = 1;
bool autorange = true;
double data_xMax;
double data_xMin;
double data_yMax;
double data_yMin;
int last_index = 0;
int graphID = -1;

const char* menu_labels[]
{
    "xRcm vs yRcm",
    "xRcm vs zRcm",
    "yRCM vs zRcm",
    "xRcm (t)",
    "yRcm (t)",
    "zRcm (t)",
    "theta (t)",
    "phi (t)",
    "phidot (t)",
    "phidot vs phi",
    "phidot and theta (t) (green, orange)",
    "w1_top and w2_top (t) (green, orange)",
    "w1_top vs w2_top",
    "Lx_lab (t)",
    "Ly_lab (t)",
    "Lz_lab (t)",
    "Energy (t)"
};

enum menuItems
{
    xvsyRCM,
    xvszRCM,
    yvszRCM,
    xRCM_T,
    yRCM_T,
    zRCM_T,
    THETA_T,
    PHI_T,
    PHIDOT_T,
    PHIDOTvsPHI,
    PHIDOTandTHETA_T,
    W1_topandW2_top_T,
    W1_topvsW2_top,
    LXlab_T,
    LYlab_T,
    LZlab_T,
    ENERGY_T
};


void setAxisRange()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(xmin_graph, xmax_graph, ymin_graph, ymax_graph);
    if(autorange) axis_boxes->set_text_values(xmin_graph, xmax_graph, ymin_graph, ymax_graph);
    glMatrixMode(GL_MODELVIEW);
}

template<class T>
void xvsy_plot(std::vector<T>& x_data, int x_index, std::vector<T>& y_data, int y_index)
{
    if (autorange_button->value() && y_data.size() && x_data.size())
    {
        if (last_index == 0)
        {
            data_xMax = x_data[0][x_index];
            data_xMin = x_data[0][x_index];
            data_yMax = y_data[0][y_index];
            data_yMin = y_data[0][y_index];
        }
        for (; last_index < y_data.size(); ++last_index)
        {
            data_xMax = (std::max)(data_xMax, x_data[last_index][x_index]);
            data_xMin = (std::min)(data_xMin, x_data[last_index][x_index]);
            data_yMax = (std::max)(data_yMax, y_data[last_index][y_index]);
            data_yMin = (std::min)(data_yMin, y_data[last_index][y_index]);
        }

        xmax_graph = data_xMax + abs(data_xMax) * 0.1;
        xmin_graph = data_xMin - abs(data_xMin) * 0.1;
        ymax_graph = data_yMax + abs(data_yMax) * 0.1;
        ymin_graph = data_yMin - abs(data_yMin) * 0.1;

        setAxisRange();
    }

    //draw axis
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(xmin_graph, 0);
    glVertex2f(xmax_graph, 0);
    glEnd();

    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(0, ymin_graph);
    glVertex2f(0, ymax_graph);
    glEnd();

    if (x_data.size() && y_data.size())
    {
        glBegin(GL_POINTS);
        glColor3f(0, 1, 0);

        for (int i = 0; i < rcm_data.size(); ++i)
        {
            glVertex2f(x_data[i][x_index], y_data[i][y_index]);
        }
        glEnd();

    }
}

template<class T>
void time_single_plot(std::vector<T>& y_data, int y_index)
{    
    if (autorange_button->value() && y_data.size())
    {
        if (last_index == 0)
        {
            data_xMax = 0+delta_t;
            data_xMin = 0;
            data_yMax = y_data[0][y_index];
            data_yMin = y_data[0][y_index];
        }
        for (; last_index < y_data.size(); ++last_index)
        {
            data_yMax = (std::max)(data_yMax, y_data[last_index][y_index]);
            data_yMin = (std::min)(data_yMin, y_data[last_index][y_index]);
        }

        xmax_graph = delta_t*y_data.size()+1;
        xmin_graph = 0 - 0.25;
        ymax_graph = data_yMax + abs(data_yMax) * 0.1;
        ymin_graph = data_yMin - abs(data_yMin) * 0.1;

        setAxisRange();
    }

    //draw axis
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(xmin_graph, 0);
    glVertex2f(xmax_graph, 0);
    glEnd();

    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(0, ymin_graph);
    glVertex2f(0, ymax_graph);
    glEnd();

    if (y_data.size())
    {
        glBegin(GL_POINTS);
        glColor3f(0, 1, 0);

        for (int i = 0; i < rcm_data.size(); ++i)
        {
            glVertex2f(i*delta_t, y_data[i][y_index]);
        }
        glEnd();

    }
}

template<class T>
void time_single_plot(std::vector<T>& y_data)
{
    if (autorange_button->value() && y_data.size())
    {
        if (last_index == 0)
        {
            data_xMax = 0 + delta_t;
            data_xMin = 0;
            data_yMax = y_data[0];
            data_yMin = y_data[0];
        }
        for (;last_index < y_data.size(); ++last_index)
        {
            data_yMax = (std::max)(data_yMax, y_data[last_index]);
            data_yMin = (std::min)(data_yMin, y_data[last_index]);
        }

        xmax_graph = delta_t * y_data.size() + 1;
        xmin_graph = 0 - 0.25;
        ymax_graph = data_yMax + abs(data_yMax) * 0.1;
        ymin_graph = data_yMin - abs(data_yMin) * 0.1;

        setAxisRange();
    }

    //draw axis
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(xmin_graph, 0);
    glVertex2f(xmax_graph, 0);
    glEnd();

    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(0, ymin_graph);
    glVertex2f(0, ymax_graph);
    glEnd();

    if (y_data.size())
    {
        glBegin(GL_POINTS);
        glColor3f(0, 1, 0);

        for (int i = 0; i < rcm_data.size(); ++i)
        {
            glVertex2f(i * delta_t, y_data[i]);
        }
        glEnd();

    }
}

template<class T>
void time_double_plot(std::vector<T>& y1_data, int y1_index, std::vector<T>& y2_data, int y2_index)
{
    if (autorange_button->value() && y1_data.size() && y2_data.size())
    {
        if (last_index == 0)
        {
            data_xMax = 0 + delta_t;
            data_xMin = 0;
            data_yMax = (std::max)(y1_data[0][y1_index], y2_data[0][y2_index]);
            data_yMin = (std::min)(y1_data[0][y1_index], y2_data[0][y2_index]);
        }
        for (; last_index < y1_data.size(); ++last_index)
        {
            data_yMax = (std::max)(data_yMax, (std::max)(y1_data[last_index][y1_index], y2_data[last_index][y2_index]));
            data_yMin = (std::min)(data_yMin, (std::min)(y1_data[last_index][y1_index], y2_data[last_index][y2_index]));
        }

        xmax_graph = delta_t * y1_data.size() + 1;
        xmin_graph = 0 - 0.25;
        ymax_graph = data_yMax + abs(data_yMax) * 0.1;
        ymin_graph = data_yMin - abs(data_yMin) * 0.1;

        setAxisRange();
    }

    //draw axis
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(xmin_graph, 0);
    glVertex2f(xmax_graph, 0);
    glEnd();

    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(0, ymin_graph);
    glVertex2f(0, ymax_graph);
    glEnd();

    if (y1_data.size() && y2_data.size())
    {
        glBegin(GL_POINTS);
        glColor3f(0, 1, 0);
        for (int i = 0; i < rcm_data.size(); ++i)
        {
            glVertex2f(i * delta_t, y1_data[i][y1_index]);
        }
        glEnd();

        glBegin(GL_POINTS);
        glColor3f(1, 0.55, 0);
        for (int i = 0; i < rcm_data.size(); ++i)
        {
            glVertex2f(i * delta_t, y2_data[i][y2_index]);
        }
        glEnd();

    }
}


void drawGraph() {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    switch (graphID)
    {
    case xvsyRCM:
        xvsy_plot(rcm_data, 0, rcm_data, 1);
        break;

    case xvszRCM:
        xvsy_plot(rcm_data, 0, rcm_data, 2);
        break;

    case yvszRCM:
        xvsy_plot(rcm_data, 1, rcm_data, 2);
        break;

    case xRCM_T:
        time_single_plot(rcm_data, 0);
        break;

    case yRCM_T:
        time_single_plot(rcm_data, 1);
        break;

    case zRCM_T:
        time_single_plot(rcm_data, 2);
        break;

    case THETA_T:
        time_single_plot(theta_phi_phidot_data, 0);
        break;

    case PHI_T:
        time_single_plot(theta_phi_phidot_data, 1);
        break;

    case PHIDOT_T:
        time_single_plot(theta_phi_phidot_data, 2);
        break;

    case PHIDOTvsPHI:
        xvsy_plot(theta_phi_phidot_data, 1, theta_phi_phidot_data, 2);
        break;

    case PHIDOTandTHETA_T:
        time_double_plot(theta_phi_phidot_data, 2, theta_phi_phidot_data, 0);
        break;

    case W1_topandW2_top_T:
        time_double_plot(w_top_data, 0, w_top_data, 1);
        break;

    case W1_topvsW2_top:
        xvsy_plot(w_top_data, 1, w_top_data, 0);
        break;

    case LXlab_T:
        time_single_plot(L_lab_data, 0);
        break;

    case LYlab_T:
        time_single_plot(L_lab_data, 1);
        break;

    case LZlab_T:
        time_single_plot(L_lab_data, 2);
        break;

    case ENERGY_T:
        time_single_plot(Energy_data);
        break;

    default:
        //draw axis
        glBegin(GL_LINES);
        glColor3f(1, 0, 0);
        glVertex2f(xmin_graph, 0);
        glVertex2f(xmax_graph, 0);
        glEnd();

        glBegin(GL_LINES);
        glColor3f(1, 0, 0);
        glVertex2f(0, ymin_graph);
        glVertex2f(0, ymax_graph);
        glEnd();
    }

    glFlush();

}


