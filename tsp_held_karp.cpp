#include <iostream>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <list>


using namespace std;

// vars
vector<pair<int, int> > points;
vector<vector<int> > adj_Matrix;
vector<vector<int> > combinations;
const int DEFAULT_VAL = 0;

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
void held_karp(unsigned long);
vector<vector<int> > gen_combinations(int, int);
int dist(int, int, int, int);
void draw_point(Colormap, char [], XColor &, int);
void draw_point(Colormap, char [], XColor &, int, int);
void draw_line(Colormap, char [], XColor &, int, int);

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
    char color_black[] = "#000000";
    char firebrick[] = "#B22222";
    colormap = DefaultColormap(dis, 0);

    // vars
    int x, y;
    unsigned long n;
    bool input = true;
    bool display_result = true;


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
                    if (event.xbutton.button == Button1 and input) {
                        // store the point
                        points.push_back(make_pair(x, y));
                        // draw the point
                        draw_point(colormap, firebrick, color, x, y);

                    }
                } // mode

                if (event.xbutton.button == Button3 and display_result) {
                    input = false;
                    // init the matrix with 0
                    n = points.size();
                    cout << "Number of points: " << n << endl;
                    adj_Matrix.resize(n, vector<int>(n, DEFAULT_VAL));
                    //printf("adj_Matrix[%d][%d]\n", (int) adj_Matrix.size(), (int) adj_Matrix[0].size());

                    // calc distances and fill the matrix
                    for (int i = 0; i < n; ++i) {
                        for (int j = i + 1; j < n; ++j) {
                            adj_Matrix[i][j] = adj_Matrix[j][i] = dist(points[i].first, points[i].second,
                                                                       points[j].first, points[j].second);
                        }
                    }

                    /*0,  2,  9, 10
                      1,  0,  6,  4
                      15, 7,  0,  8
                      6,  3, 12,  0
                      */

                    /*for (int i = 0; i < n; ++i) {
                        for (int j = 0; j < n; ++j) {
                            printf("%d,",adj_Matrix[i][j]);
                        }
                        cout << "\n";
                    }*/

                    /*for (int i = 0; i < points.size(); ++i) {
                        printf("x = %d, y = %d \n", points[i].first, points[i].second);
                    }*/




                    if (n <= 15) {
                        held_karp(n);
                        display_result = false;
                    } else {
                        // break them in sets of 15

                    }
                }
                break;

            default:
                break;

        }
    }
}
#pragma clang diagnostic pop

void held_karp(unsigned long n) {
    // init colors
    XColor color;
    Colormap colormap;
    char dodger_blue[] = "#1E90FF";
    colormap = DefaultColormap(dis, 0);

    map<pair<int, int>, pair<int, int> > C;
    // then do this
    /*# Set transition cost from initial state
        for k in range(1, n):
            C[(1 << k, k)] = (dists[0][k], 0)*/

    for (int k = 1; k < n; ++k) {
        C.insert(make_pair(make_pair(1 << k, k), make_pair(adj_Matrix[0][k], 0)));
    }

    /*for(map<pair<int, int>, pair<int, int> >::const_iterator it = C.begin(); it != C.end(); ++it) {
        std::cout << "{(" << it->first.first << "," << it->first.second << ") : ("
                  << it->second.first << "," << it->second.second << ")} ";
    }
    cout << "\n";*/

    for (int i = 2; i < n; ++i) {
        // here generate combinations. sorted
        combinations = gen_combinations((int) (n - 1), i);

        //cout << "size. = " << combinations.size() << " size[0]. = " << combinations[0].size() << endl;

        for (int j = 0; j < combinations.size(); ++j) {

            int bits = 0;
            // Set bits for all nodes in this subset
            for (int k = 0; k < combinations[j].size(); ++k) {

                bits |= 1 << combinations[j][k];

                //printf("[%d]", combinations[j][k]);
                //printf("%d\n", bits);
            }
            //cout << endl;

            /*
             *  // Find the lowest cost to get to this subset
                for k in subset:
                    prev = bits & ~(1 << k)

                    res = []
                    for m in subset:
                        if m == 0 or m == k:
                            continue
                        res.append((C[(prev, m)][0] + dists[m][k], m))
                    C[(bits, k)] = min(res)*/

            for (int k = 0; k < combinations[j].size(); ++k) {
                int prev = bits & ~(1 << combinations[j][k]);

                //cout << "prev: " << prev << endl;

                list<pair<int, int> > integer_list;

                for (int m = 0; m < combinations[j].size(); ++m) {
                    //cout << combinations[j][m] << " m and k " << combinations[j][k] << endl;

                    if (combinations[j][m] == 0 || combinations[j][m] == combinations[j][k]) {
                        //cout << "yes\n";
                        continue;
                    }
                    int list_i = 0;

                    map<pair<int, int>, pair<int, int> >::const_iterator it;
                    it = C.find(make_pair(prev, combinations[j][m]));
                    if (it != C.end()) {
                        // ... use itr to read the values ...
                        list_i = it->second.first;
                        //cout << list_i << endl;
                    }
                    //cout << prev << " , " << combinations[j][m] << endl;
                    //cout << list_i << " + " << adj_Matrix[combinations[j][m]][combinations[j][k]] << " \n";
                    //cout << "(" << list_i + adj_Matrix[combinations[j][m]][combinations[j][k]] << "," << combinations[j][m] << ")\n";
                    integer_list.push_back(make_pair(list_i + adj_Matrix[combinations[j][m]][combinations[j][k]], combinations[j][m]));
                }
                // min
                /*list<pair<int, int> >::iterator iter;
                for (iter = integer_list.begin(); iter != integer_list.end(); ++iter) {
                    cout << (*iter).first << "," << (*iter).second << endl;
                }*/

                list<pair<int, int> >::iterator it = min_element(integer_list.begin(), integer_list.end());
                //cout << "min: " << (*it).first << "," << (*it).second << endl;
                //C[(bits, k)] = min(res)*/

                C.insert(make_pair(make_pair(bits, combinations[j][k]), make_pair((*it).first, (*it).second) ));

                /*for(map<pair<int, int>, pair<int, int> >::const_iterator itr = C.begin(); itr != C.end(); ++itr) {
                    std::cout << "{(" << itr->first.first << "," << itr->first.second << ") : ("
                              << itr->second.first << "," << itr->second.second << ")} ";
                }
                cout << "\n";*/
                integer_list.clear();

            }


        }

    }
    // We're interested in all bits but the least significant (the start state)
    int bits = (int)(pow(2, (double)n) - 1) - 1;

    //cout << bits << endl;

    // calculate
    list<pair<int, int> > res;

    for (int k = 1; k < n; ++k) {
        int list_i = 0;
        map<pair<int, int>, pair<int, int> >::const_iterator it;
        it = C.find(make_pair(bits, k));
        if (it != C.end()) {
            // ... use itr to read the values ...
            list_i = it->second.first;
            //cout << list_i << endl;
        }
        res.push_back(make_pair(list_i + adj_Matrix[k][0], k));
    }
    /*list<pair<int, int> >::iterator iter;
    for (iter = res.begin(); iter != res.end(); ++iter) {
        cout << (*iter).first << "," << (*iter).second << endl;
    }*/

    list<pair<int, int> >::iterator it = min_element(res.begin(), res.end());
    //cout << "min: " << (*it).first << "," << (*it).second << endl;
    int opt = (*it).first, parent = (*it).second;

    //cout << "parent: " << parent << endl;

    //Backtrack to find full path
    list<int> path;

    for (int i = 0; i < n - 1; ++i) {
        path.push_back(parent);
        int new_bits = bits & ~(1 << parent);
        //cout << "bits: " << bits << " parent: " << parent << endl;
        int list_parent = 0;
        map<pair<int, int>, pair<int, int> >::const_iterator iter;
        iter = C.find(make_pair(bits, parent));
        if (iter != C.end()) {
            // ... use itr to read the values ...
            list_parent = iter->second.second;
            //cout << list_parent << endl;
        }
        parent = list_parent;
        bits = new_bits;

    }
    // add manual starting and ending points.
    path.push_back(0);
    path.push_front(0);

    path.reverse();

    //int *path_index = (int *) malloc(n + 1);
    int path_index[n + 1];
    int index = 0;

    //cout << "opt: " << opt << " ( ";

    for(std::list<int>::iterator list_iter = path.begin();
        list_iter != path.end(); list_iter++)
    {
        //std::cout<<*list_iter << " ";
        path_index[index] = *list_iter;
        index++;
    }
    //cout << ")\n";

    for (int l = 0; l < index - 1; ++l) {
        //cout << l << ": " << path_index[l] << endl;
        draw_line(colormap, dodger_blue, color, path_index[l], path_index[l + 1]);
    }

    cout << "Total length of the tour: " << opt << endl;

    // clear containers
    C.clear();
    path.clear();
    res.clear();
}

vector<vector<int> > gen_combinations(int n, int r) {
    //http://stackoverflow.com/questions/9430568/generating-combinations-in-c
    int count = 0;
    vector<vector<int> > container;
    vector<int> row;
    vector<bool> v((unsigned long) n);
    fill(v.begin(), v.begin() + r, true);

    do {
        for (int i = 0; i < n; ++i) {
            if (v[i]) {
                count++;
                //std::cout << (i + 1) << " ";
                row.push_back(i + 1);
                if (count == r) {
                    container.push_back(row);
                    row.clear();
                    count = 0;
                }
            }
        }
        //std::cout << "\n";
    } while (prev_permutation(v.begin(), v.end()));

    return container;
}

int dist(int x1, int y1, int x2, int y2) {
    return (int) sqrt(pow(x2 - x1, 2.0f) + pow(y2 - y1, 2.0f));
}

void draw_point(Colormap colormap, char color_name[], XColor &color, int index) {
    XParseColor(dis, colormap, color_name, &color);
    XAllocColor(dis, colormap, &color);
    XSetForeground(dis, gc, color.pixel);
    XFillArc(dis, win, gc, points[index].first-(15/2), points[index].second-(15/2), 8, 8, 0, 360*64);
}

void draw_point(Colormap colormap, char color_name[], XColor &color, int x, int y) {
    XParseColor(dis, colormap, color_name, &color);
    XAllocColor(dis, colormap, &color);
    XSetForeground(dis, gc, color.pixel);
    XFillArc(dis, win, gc, x-(15/2), y-(15/2), 8, 8, 0, 360*64);
}

void draw_line(Colormap colormap, char color_name[], XColor &color, int i, int j) {
    XParseColor(dis, colormap, color_name, &color);
    XAllocColor(dis, colormap, &color);
    XSetForeground(dis, gc, color.pixel);
    XDrawLine(dis,win,gc,points[i].first,points[i].second,
              points[j].first,points[j].second);
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
    combinations.clear();
    points.clear();
    adj_Matrix.clear();
    XFreeGC(dis, gc);
    XDestroyWindow(dis,win);
    XCloseDisplay(dis);
    exit(1);
};

void redraw() {
    XClearWindow(dis, win);
}