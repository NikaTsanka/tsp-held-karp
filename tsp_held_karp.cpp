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

// global vars
vector<pair<int, int> > coordinates;
vector<vector<int> > adj_Matrix;
vector<vector<int> > path_indices;
vector<pair<int, int> > connecting_indices;
vector<vector<pair<int, int> > > points;
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
int held_karp(unsigned long, int);
vector<vector<int> > gen_combinations(int, int);
int dist(int, int, int, int);
bool compare_pair(pair<int, int> &, pair<int, int> &);
int connect_slices();
void draw_point(Colormap, char [], XColor &, int);
void draw_point(Colormap, char [], XColor &, int, int);
void draw_line(Colormap, char [], XColor &, int, int, vector<pair<int, int> >);
void draw_line(Colormap c, char [], XColor &, int, int, int, int);
void draw_line(Colormap c, char [], XColor &, int, int);

int main(int argc, char *argv[]) {
    if (argc == 2) {
        printf("Mode: Command-Line\n");
        int x, y;
        ifstream input_file(argv[1]);
        while(input_file >> x >> y) {
            coordinates.push_back(make_pair(x, y));
        }
        input_file.close();
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
    char crimson[] = "#DC143C";
    char navy_blue[] = "#000080";
    char blue[] = "#0000FF";
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
                        for (int i = 0; i < coordinates.size(); ++i) {
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
                        coordinates.push_back(make_pair(x, y));
                        // draw the point
                        draw_point(colormap, navy_blue, color, x, y);
                    }
                } // mode

                if (event.xbutton.button == Button3 and display_result) {
                    input = false;
                    // init the matrix with 0
                    // n number of coordinates
                    n = coordinates.size();
                    cout << "Number of coordinates: " << n << endl;
                    if (n <= 15) {
                        adj_Matrix.resize(n, vector<int>(n, DEFAULT_VAL));
                        // calc distances and fill the matrix
                        for (int i = 0; i < n; ++i) {
                            for (int j = i + 1; j < n; ++j) {
                                adj_Matrix[i][j] = adj_Matrix[j][i] = dist(coordinates[i].first, coordinates[i].second,
                                                                           coordinates[j].first, coordinates[j].second);
                            }
                        }
                        // do the calculation
                        held_karp(n, 1);
                        // draw the path
                        for (int k = 0; k < path_indices.size(); ++k) {
                            for (int i = 0; i < path_indices[k].size() - 1; ++i) {
                                draw_line(colormap, crimson, color, path_indices[k][i], path_indices[k][i + 1]);
                            }
                        }
                        display_result = false;
                    } else {
                        // first sort them
                        sort(coordinates.begin(), coordinates.end(), compare_pair);
                        // break them in sets of 15
                        int k = (int) (n / 15);
                        int dif = (int) (n - (k * 15));
                        int div = k;
                        int end = 0;
                        unsigned long points_size = 0;
                        int total_tour = 0;
                        int g = (dif == 0) ? k : k + 1;

                        for (int i = 0; i < g; ++i) {
                            // slice and fill matrix with it and calculate
                            div--;
                            if (div < 0) {
                                end = 0;
                                points_size = (unsigned long) dif;
                            } else {
                                end = (int) (n - (n - (((div) * 15) + dif)));
                                points_size = 15;
                            }
                            vector<pair<int, int> > points_slice;

                            copy(coordinates.begin() + (i * 15), coordinates.end() - end, back_inserter(points_slice));

                            adj_Matrix.resize(points_size, vector<int>(points_size, DEFAULT_VAL));
                            // calc distances and fill the matrix
                            for (int indx = 0; indx < points_size; ++indx) {
                                for (int j = indx + 1; j < points_size; ++j) {
                                    adj_Matrix[indx][j] = adj_Matrix[j][indx] = dist(points_slice[indx].first, points_slice[indx].second,
                                                                                     points_slice[j].first, points_slice[j].second);
                                }
                            }
                            // do the calculation
                            total_tour += held_karp(points_size, i + 1);
                            // update
                            points.push_back(points_slice);
                            // clear
                            points_slice.clear();
                            adj_Matrix.clear();
                        }
                        // update total
                        total_tour += connect_slices();

                        cout << "Total tour: " << total_tour << " (includes new connecting edges)" << endl;
                        // draw
                        int num = (int) connecting_indices.size();
                        for (int m = 0; m < path_indices.size(); ++m) {
                            for (int i = 0; i < path_indices[m].size() - 1; ++i) {
                                if (m == 0) {
                                    if (!(i == 0 && i == connecting_indices[m].first)) {
                                        draw_line(colormap, crimson, color, path_indices[m][i], path_indices[m][i + 1],
                                                  points[m]);
                                    }
                                    if (i + 1 == connecting_indices[m].first) {
                                        for (int j = 0; j < path_indices[m + 1].size(); ++j) {
                                            if (j == connecting_indices[m].second) {
                                                draw_line(colormap, blue, color,
                                                          points[m][path_indices[m][i + 1]].first,
                                                          points[m][path_indices[m][i + 1]].second,
                                                          points[m + 1][path_indices[m + 1][j]].first,
                                                          points[m + 1][path_indices[m + 1][j]].second);

                                                draw_line(colormap, blue, color,
                                                          points[m][path_indices[m][i + 2]].first,
                                                          points[m][path_indices[m][i + 2]].second,
                                                          points[m + 1][path_indices[m + 1][j + 1]].first,
                                                          points[m + 1][path_indices[m + 1][j + 1]].second);
                                            }
                                        }
                                        // jump
                                        i+=1;
                                    }
                                } else {
                                    if (m != num) {
                                        if (!(i == 0 && i == connecting_indices[m].first ||
                                              i == 0 && i == connecting_indices[m - 1].second)) {
                                            draw_line(colormap, crimson, color, path_indices[m][i], path_indices[m][i + 1],
                                                      points[m]);
                                        }
                                        if (i + 1 == connecting_indices[m - 1].second) {
                                            i+=1;
                                        }
                                        if (i + 1 == connecting_indices[m].first) {
                                            //continue;
                                            for (int j = 0; j < path_indices[m + 1].size(); ++j) {
                                                if (j == connecting_indices[m].second) {
                                                    draw_line(colormap, blue, color,
                                                              points[m][path_indices[m][i + 1]].first,
                                                              points[m][path_indices[m][i + 1]].second,
                                                              points[m + 1][path_indices[m + 1][j]].first,
                                                              points[m + 1][path_indices[m + 1][j]].second);

                                                    draw_line(colormap, blue, color,
                                                              points[m][path_indices[m][i + 2]].first,
                                                              points[m][path_indices[m][i + 2]].second,
                                                              points[m + 1][path_indices[m + 1][j + 1]].first,
                                                              points[m + 1][path_indices[m + 1][j + 1]].second);
                                                }
                                            }
                                            i+=1;
                                        }
                                    }
                                    if (num == m) {
                                        if (!(i == 0 && i == connecting_indices[m - 1].second)) {
                                            draw_line(colormap, crimson, color, path_indices[m][i], path_indices[m][i + 1],
                                                      points[m]);
                                        }
                                        if (i + 1 == connecting_indices[m - 1].second) {
                                            i+=1;
                                        }
                                    }
                                }
                            }
                        }
                        display_result = false;
                    }
                }
                break;

            default:
                break;
        }
    }
}
#pragma clang diagnostic pop

int held_karp(unsigned long n, int group) {
    map<pair<int, int>, pair<int, int> > map_container;
    // then do this
    for (int k = 1; k < n; ++k) {
        map_container.insert(make_pair(make_pair(1 << k, k), make_pair(adj_Matrix[0][k], 0)));
    }

    vector<vector<int> > combinations;
    for (int i = 2; i < n; ++i) {
        // here generate combinations. sorted
        combinations = gen_combinations((int) (n - 1), i);
        for (int j = 0; j < combinations.size(); ++j) {
            int bits = 0;
            // set bits for all nodes in this subset
            for (int k = 0; k < combinations[j].size(); ++k) {
                bits |= 1 << combinations[j][k];
            }
            // find the lowest cost to get to this subset
            for (int k = 0; k < combinations[j].size(); ++k) {
                int prev = bits & ~(1 << combinations[j][k]);

                list<pair<int, int> > integer_list;
                for (int m = 0; m < combinations[j].size(); ++m) {
                    if (combinations[j][m] == 0 || combinations[j][m] == combinations[j][k]) {
                        continue;
                    }
                    int list_i = 0;

                    map<pair<int, int>, pair<int, int> >::const_iterator it;
                    it = map_container.find(make_pair(prev, combinations[j][m]));
                    if (it != map_container.end()) {
                        // use itr to read the values
                        list_i = it->second.first;
                    }
                    integer_list.push_back(make_pair(list_i + adj_Matrix[combinations[j][m]][combinations[j][k]], combinations[j][m]));
                }
                // min
                list<pair<int, int> >::iterator it = min_element(integer_list.begin(), integer_list.end());
                map_container.insert(make_pair(make_pair(bits, combinations[j][k]), make_pair((*it).first, (*it).second) ));
                integer_list.clear();
            }
        }
    }
    // we're interested in all bits but the least significant (the start state)
    int bits = (int)(pow(2, (double)n) - 1) - 1;

    // calculate
    list<pair<int, int> > list_of_pairs_result;
    for (int k = 1; k < n; ++k) {
        int list_i = 0;
        map<pair<int, int>, pair<int, int> >::const_iterator it;
        it = map_container.find(make_pair(bits, k));
        if (it != map_container.end()) {
            // use itr to read the values
            list_i = it->second.first;
        }
        list_of_pairs_result.push_back(make_pair(list_i + adj_Matrix[k][0], k));
    }
    list<pair<int, int> >::iterator it = min_element(list_of_pairs_result.begin(), list_of_pairs_result.end());
    int opt = (*it).first, parent = (*it).second;

    // backtrack to find full path
    list<int> path;
    for (int i = 0; i < n - 1; ++i) {
        path.push_back(parent);
        int new_bits = bits & ~(1 << parent);
        int list_parent = 0;
        map<pair<int, int>, pair<int, int> >::const_iterator iter;
        iter = map_container.find(make_pair(bits, parent));
        if (iter != map_container.end()) {
            // use itr to read the values
            list_parent = iter->second.second;
        }
        parent = list_parent;
        bits = new_bits;

    }
    // add manual starting and ending coordinates.
    path.push_back(0);
    path.push_front(0);

    path.reverse(); // doesn't really matter.

    vector<int> path_slice;
    for(std::list<int>::iterator list_iter = path.begin();
        list_iter != path.end(); list_iter++) {
        path_slice.push_back(*list_iter);
    }

    cout << "Total length of the tour #" << group << ": "  << opt << endl;
    path_indices.push_back(path_slice);

    // clear containers
    combinations.clear();
    map_container.clear();
    path.clear();
    list_of_pairs_result.clear();
    return opt;
}

vector<vector<int> > gen_combinations(int n, int r) {
    int count = 0;
    vector<vector<int> > container;
    vector<int> row;
    vector<bool> v((unsigned long) n);
    fill(v.begin(), v.begin() + r, true);

    do {
        for (int i = 0; i < n; ++i) {
            if (v[i]) {
                count++;
                row.push_back(i + 1);
                if (count == r) {
                    container.push_back(row);
                    row.clear();
                    count = 0;
                }
            }
        }
    } while (prev_permutation(v.begin(), v.end()));
    return container;
}

int dist(int x1, int y1, int x2, int y2) {
    return (int) sqrt(pow(x2 - x1, 2.0f) + pow(y2 - y1, 2.0f));
}

bool compare_pair(pair<int, int>& i, pair<int, int>& j) {
    return i.first < j.first;
}

int connect_slices() {
    // finds optimal connections between the groups
    int total = 0;
    for (int g = 0; g < path_indices.size() - 1; ++g) { // number of groups
        // pick an edge from group i
        int pi = 0;
        int qj = 0;

        int overall_min_distance = 99999;
        int min_distance = 99999;
        for (int i = 0; i < path_indices[g].size() - 1; ++i) {
            int pix =  points[g][path_indices[g][i]].first,     piy = points[g][path_indices[g][i]].second;
            int pi1x = points[g][path_indices[g][i + 1]].first, pi1y = points[g][path_indices[g][i + 1]].second;
            for (int j = 0; j < path_indices[g + 1].size() - 1; ++j) {
                int qjx =  points[g + 1][path_indices[g + 1][j]].first,     qjy = points[g + 1][path_indices[g + 1][j]].second;
                int qj1x = points[g + 1][path_indices[g + 1][j + 1]].first, qj1y = points[g + 1][path_indices[g + 1][j + 1]].second;
                int dist_pi_qj = dist(pix, piy, qjx, qjy);
                int dist_pi1_qj1 = dist(pi1x, pi1y, qj1x, qj1y);
                int distance =  dist_pi_qj + dist_pi1_qj1;

                if (distance < min_distance) {
                    min_distance = distance;
                    pi = i;
                    qj = j;
                }
            }
            if (min_distance < overall_min_distance) {
                overall_min_distance = min_distance;
            }
        }
        connecting_indices.push_back(make_pair(pi, qj));
        total += overall_min_distance;
    }
    return total;
}

void draw_point(Colormap colormap, char color_name[], XColor &color, int index) {
    XParseColor(dis, colormap, color_name, &color);
    XAllocColor(dis, colormap, &color);
    XSetForeground(dis, gc, color.pixel);
    XFillArc(dis, win, gc, coordinates[index].first-(15/2), coordinates[index].second-(15/2), 8, 8, 0, 360*64);
}

void draw_point(Colormap colormap, char color_name[], XColor &color, int x, int y) {
    XParseColor(dis, colormap, color_name, &color);
    XAllocColor(dis, colormap, &color);
    XSetForeground(dis, gc, color.pixel);
    XFillArc(dis, win, gc, x-(15/2), y-(15/2), 8, 8, 0, 360*64);
}

void draw_line(Colormap colormap, char color_name[], XColor &color, int i, int j, vector<pair<int, int> > current_points) {
    XParseColor(dis, colormap, color_name, &color);
    XAllocColor(dis, colormap, &color);
    XSetForeground(dis, gc, color.pixel);
    XDrawLine(dis,win,gc,current_points[i].first,current_points[i].second,
              current_points[j].first,current_points[j].second);
}

void draw_line(Colormap colormap, char color_name[], XColor &color, int ix, int iy, int jx, int jy) {
    XParseColor(dis, colormap, color_name, &color);
    XAllocColor(dis, colormap, &color);
    XSetForeground(dis, gc, color.pixel);
    XDrawLine(dis,win,gc,ix,iy,jx,jy);
}

void draw_line(Colormap colormap, char color_name[], XColor &color, int i, int j) {
    XParseColor(dis, colormap, color_name, &color);
    XAllocColor(dis, colormap, &color);
    XSetForeground(dis, gc, color.pixel);
    XDrawLine(dis,win,gc,coordinates[i].first,coordinates[i].second,coordinates[j].first,coordinates[j].second);
}

void init_x() {
    // get the colors black and white (see section for details)
    unsigned long black,white;

    dis=XOpenDisplay((char *)0);
    screen=DefaultScreen(dis);
    black=BlackPixel(dis,screen), white=WhitePixel(dis, screen);

    win=XCreateSimpleWindow(dis, DefaultRootWindow(dis), 0, 0,
                            500, 500, 5, black, white);

    XSetStandardProperties(dis,win,"TSP Simulation - Held-Karp Algorithm",NULL,None,NULL,0,NULL);
    XSelectInput(dis, win, ExposureMask|ButtonPressMask|KeyPressMask);

    gc=XCreateGC(dis, win, 0,0);

    XSetBackground(dis,gc,white);
    XSetForeground(dis,gc,black);

    XClearWindow(dis, win);
    XMapRaised(dis, win);
};

void close_x() {
    // clear();
    coordinates.clear();
    adj_Matrix.clear();
    path_indices.clear();
    connecting_indices.clear();
    points.clear();

    XFreeGC(dis, gc);
    XDestroyWindow(dis,win);
    XCloseDisplay(dis);
    cout << "Program was terminated successfully";
    exit(1);
};

void redraw() {
    XClearWindow(dis, win);
}