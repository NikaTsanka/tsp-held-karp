#include <iostream>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <limits>
#include <fstream>

using namespace std;

vector<pair<int, int> > points;

// here are our X variables
Display *dis;
int screen;
Window win;
GC gc;

// here are our X routines declared
void init_x();
void close_x();
void redraw();

// my methods
void drawing_board(bool);
void draw_point(Colormap, char [], XColor &, int);
void draw_point(Colormap, char [], XColor &, int, int);

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-noreturn"
void drawing_board(bool mode) {
    // the XEvent declaration
    XEvent event;
    // to handle KeyPress Events
    KeySym key;
    // for KeyPress Events
    char text[255];
    // create and open the window
    init_x();
    // init colors
    XColor color;
    Colormap colormap;
    char dodger_blue[] = "#1E90FF";
    char color_black[] = "#000000";
    char firebrick[] = "#B22222";
    colormap = DefaultColormap(dis, 0);

    // vars
    int x, y;

    // start the main loop
    cout << "You can use \'q\' or \'Q\' to terminate the program at any time. :)\n";
    while(1) {
        /* get the next event and stuff it into our event variable.
           Note:  only events we set the mask for are detected!
        */
        XNextEvent(dis, &event);
        switch (event.type) {
            case Expose:
                if (event.type == Expose && event.xexpose.count == 0) {
                    // the window was exposed redraw it
                    redraw();
                    // now put the servers on the screen
                    if (mode) {
                        for (int i = 0; i < points.size(); ++i) {
                            draw_point(colormap, color_black, color, i);
                        }
                    }
                }
                break;

            case KeyPress:
                if (XLookupString(&event.xkey, text, 255, &key, 0) == 1) {
                    /* use the XLookupString routine to convert the invent
                       KeyPress data into regular text.  Weird but necessary...
                    */
                    if (text[0] == 'q' || text[0] == 'Q') {
                        close_x();
                    }
                    printf("You pressed the %c key!\n", text[0]);
                }
                break;

            case ButtonPress:
                // tell where the mouse Button was Pressed
                x = event.xbutton.x;
                y = event.xbutton.y;
                if (not mode) {
                    // catch left mouse click
                    if (event.xbutton.button == Button1) {
                        // store the point
                        points.push_back(make_pair(x, y));
                        // draw the point
                        draw_point(colormap, firebrick, color, x, y);

                    } else if (event.xbutton.button == Button3 ) {
                        // calculate here
                        /*cout << "Number of points: " << points.size() << endl;

                        for (int i = 0; i < points.size(); ++i) {
                            printf("x = %d, y = %d \n", points[i].first, points[i].second);
                        }*/
                    }
                } // mode
                break;

            default:
                break;

        }
    }
}
#pragma clang diagnostic pop
int main(int argc, char *argv[]) {
    //std::cout << "Hello, World!" << std::endl;
    if (argc == 2) {
        printf("Mode: Command-Line\n");

        int x, y;

        ifstream input_file(argv[1]);

        while(input_file >> x >> y) {
            points.push_back(make_pair(x, y));
        }
        input_file.close();

        /*cout << "Number of points: " << points.size() << endl;

        for (int i = 0; i < points.size(); ++i) {
            printf("x = %d, y = %d \n", points[i].first, points[i].second);
        }*/

        drawing_board(true);

    } else {
        printf("Mode: Interactive\n");
        drawing_board(false);

    }


}

void draw_point(Colormap colormap, char color_name[], XColor &color, int x, int y) {
    XParseColor(dis, colormap, color_name, &color);
    XAllocColor(dis, colormap, &color);
    XSetForeground(dis, gc, color.pixel);
    XFillArc(dis, win, gc, x-(15/2), y-(15/2), 8, 8, 0, 360*64);
}

void draw_point(Colormap colormap, char color_name[], XColor &color, int server) {
    XParseColor(dis, colormap, color_name, &color);
    XAllocColor(dis, colormap, &color);
    XSetForeground(dis, gc, color.pixel);
    XFillArc(dis, win, gc, points[server].first-(15/2), points[server].second-(15/2), 8, 8, 0, 360*64);
}

void init_x() {
    // get the colors black and white (see section for details)
    unsigned long black,white;

    dis=XOpenDisplay((char *)0);
    screen=DefaultScreen(dis);
    black=BlackPixel(dis,screen), white=WhitePixel(dis, screen);

    win=XCreateSimpleWindow(dis, DefaultRootWindow(dis), 0, 0,
                            700, 700, 5, black, white);

    XSetStandardProperties(dis,win,"TSP Simulation - Held-Karp Algorithm",NULL,None,NULL,0,NULL);
    XSelectInput(dis, win, ExposureMask|ButtonPressMask|KeyPressMask);

    gc=XCreateGC(dis, win, 0,0);

    XSetBackground(dis,gc,white);
    XSetForeground(dis,gc,black);

    XClearWindow(dis, win);
    XMapRaised(dis, win);
};

void close_x() {
    points.clear();
    XFreeGC(dis, gc);
    XDestroyWindow(dis,win);
    XCloseDisplay(dis);
    exit(1);
};

void redraw() {
    XClearWindow(dis, win);
}