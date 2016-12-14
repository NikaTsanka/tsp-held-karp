/*
 * CSC I0600 - Fundamental Algorithms - Fall Semester 2016
 * Homework Project 3
 * Nika Tsankashvili
 * */

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
struct Points {
    int x, y;
};
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

// new methods
void drawing_board(bool);
int held_karp(unsigned long);
vector<vector<int> > gen_combinations(int, int);
int dist(int, int, int, int);
bool compare_pair(const pair<int, int> &, const pair<int, int> &);
bool compare_dist(const Points &, const Points &);
bool check_intersection(int, int, int, int, int, int, int, int);
double det(int, int, int, int, int, int);
int connect_slices(int);
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

// new methods
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
    colormap = DefaultColormap(dis, 0);

    // vars
    int x, y;
    unsigned long n;
    bool input = true;
    bool display_result = true;

    // start the main loop
    cout << "You can use \'q\' or \'Q\' to terminate the program at any time. :)\n";
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-noreturn"
    while(1) {
        /* get the next event and stuff it into our event variable.
           Note:  only events we set the mask for are detected!
        */
        XNextEvent(dis, &event);
        switch (event.type) {
            case Expose:
                if (event.type == Expose && event.xexpose.count == 0) {
                    // the window was exposed redraw it
                    //redraw();
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

                    int total_tour = 0;

                    // init the matrix with 0
                    // n number of coordinates
                    n = coordinates.size();
                    cout << "Number of coordinates: " << n << endl;

                    if (n == 0) {
                        cout << "No coordinates given.\n";
                        close_x();
                    }

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
                        total_tour = held_karp(n);
                        // draw the path
                        for (int k = 0; k < path_indices.size(); ++k) {
                            for (int i = 0; i < path_indices[k].size() - 1; ++i) {
                                draw_line(colormap, crimson, color, path_indices[k][i], path_indices[k][i + 1]);
                            }
                        }

                        total_tour += dist(coordinates[path_indices[0][0]].first, coordinates[path_indices[0][0]].second,
                                           coordinates[path_indices[0][1]].first, coordinates[path_indices[0][1]].second);
                        total_tour += dist(coordinates[path_indices[0].size()].first, coordinates[path_indices[0].size()].second,
                                           coordinates[path_indices[0].size() - 1].first, coordinates[path_indices[0].size() - 1].second);

                        cout << "Total tour: " << total_tour << endl;

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
                            total_tour += held_karp(points_size);

                            // update
                            points.push_back(points_slice);
                            // clear
                            points_slice.clear();
                            adj_Matrix.clear();
                        }
                        // update total

                        total_tour = connect_slices(total_tour);

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
#pragma clang diagnostic pop
}
int held_karp(unsigned long n) {
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

                list<Points> integer_list;
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
                    Points p;
                    p.x = list_i + adj_Matrix[combinations[j][m]][combinations[j][k]], p.y = combinations[j][m];
                    integer_list.push_back(p);
                }
                // min
                list<Points>::iterator it = min_element(integer_list.begin(), integer_list.end(), compare_dist);
                map_container.insert(make_pair(make_pair(bits, combinations[j][k]), make_pair((*it).x, (*it).y)));
                integer_list.clear();
            }
        }
    }
    // we're interested in all bits but the least significant (the start state)
    int bits = (int)(pow(2, (double)n) - 1) - 1;

    // calculate
    list<Points> list_of_pairs_result;
    for (int k = 1; k < n; ++k) {
        int list_i = 0;
        map<pair<int, int>, pair<int, int> >::const_iterator it;
        it = map_container.find(make_pair(bits, k));
        if (it != map_container.end()) {
            // use itr to read the values
            list_i = it->second.first;
        }
        Points p;
        p.x = list_i + adj_Matrix[k][0], p.y = k;
        list_of_pairs_result.push_back(p);
    }
    list<Points>::iterator it = min_element(list_of_pairs_result.begin(), list_of_pairs_result.end(), compare_dist);
    int opt = (*it).x, parent = (*it).y;

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

    vector<int> path_slice;
    for(std::list<int>::iterator list_iter = path.begin();
        list_iter != path.end(); list_iter++) {
        path_slice.push_back(*list_iter);
    }

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
bool compare_pair(const pair<int, int>& i, const pair<int, int>& j) {
    return i.first < j.first;
}
bool compare_dist(const Points &lhs, const Points &rhs) {
    return lhs.x < rhs.x;
}
double det(int px, int py, int qx, int qy, int rx, int ry) {
    return (px*qy)+(py*rx)+(qx*ry)-(qy*rx)-(px*ry)-(py*qx);
};
bool check_intersection(int p1x, int p1y, int q1x, int q1y, int p2x, int p2y, int q2x, int q2y) {
    /*Determinant formula
      If (((det (p, q , r) * det(p q s) ) < 0)
      And ((Det (s r p) * det ( s r q))< 0)
      Then the edge pq intersects sr
     */
    double pqr = 0, pqs = 0, srp = 0, srq = 0;
    // p1 = p, q1 = q, r = p2, s = q2

    pqr = det(p1x, p1y, // p
              q1x, q1y, // q
              p2x, p2y); // r

    pqs = det(p1x, p1y, // p
              q1x, q1y, // q
              q2x, q2y); // s

    srp = det(q2x, q2y, //s
              p2x, p2y, //r
              p1x, p1y); // p

    srq = det(q2x, q2y, //s
              p2x, p2y, //r
              q1x, q1y); //q

    if (((pqr * pqs) < 0) && ((srp * srq) < 0)) {
        return true;
    } else {
        return false;
    }
}
int connect_slices(int total) {
    int temp_total;
    // init colors
    XColor color;
    Colormap colormap;
    char blue[] = "#0000FF";
    colormap = DefaultColormap(dis, 0);
    int skip_index = 99999;
    int shorter_dist = false;
    // finds optimal connections between the groups
    //int total = 0;
    for (int g = 0; g < path_indices.size() - 1; ++g) { // number of groups
        // pick an edge from group i
        int pi = 0;
        int min_pi1 = 0;
        int qj = 0;
        int min_qj1 = 0;

        int overall_min_distance = 99999;//, overall_min_distance_j = 99999;
        int min_distance = 99999;//, min_distance_j = 99999;
        for (int i = 0; i < path_indices[g].size() - 1; ++i) {
            if (i == skip_index) {
                continue;
            }
            // edge
            int pix =  points[g][path_indices[g][i]].first,     piy = points[g][path_indices[g][i]].second;
            int pi1x = points[g][path_indices[g][i + 1]].first, pi1y = points[g][path_indices[g][i + 1]].second;
            for (int j = 0; j < path_indices[g + 1].size() - 1; ++j) {
                temp_total = total;
                // edge
                int qjx =  points[g + 1][path_indices[g + 1][j]].first,     qjy = points[g + 1][path_indices[g + 1][j]].second;
                int qj1x = points[g + 1][path_indices[g + 1][j + 1]].first, qj1y = points[g + 1][path_indices[g + 1][j + 1]].second;

                int dist_pi_qj = dist(pix, piy, qjx, qjy);
                int dist_pi1_qj1 = dist(pi1x, pi1y, qj1x, qj1y);

                int dist_pi_qj1 = dist(pix, piy, qj1x, qj1y);
                int dist_pi1_qj = dist(pi1x, pi1y, qjx, qjy);

                // minus cut edges
                temp_total -= dist(pix, piy, pi1x, pi1y);
                temp_total -= dist(qjx, qjy, qj1x, qj1y);

                // add the new edges that don't cross each other.
                if (check_intersection(pix, piy, qjx, qjy, pi1x, pi1y, qj1x, qj1y)) {
                    // if they cross try other
                    int distance_j = dist_pi_qj1 + dist_pi1_qj;
                    temp_total += distance_j;
                    if (temp_total < min_distance) {
                        min_distance = temp_total;
                        pi = i;
                        qj = j;
                    }
                } else {
                    // check if the other way and if it also doesn't cross
                    // pick the shortest
                    int distance_j_temp = 99999;
                    if (!check_intersection(pix, piy, qj1x, qj1y, pi1x, pi1y, qjx, qjy)) {
                        distance_j_temp = dist_pi_qj1 + dist_pi1_qj;
                    }
                    int distance_i = dist_pi_qj + dist_pi1_qj1;

                    if (distance_j_temp < distance_i) {
                        temp_total += distance_j_temp;
                        if (temp_total < min_distance) {
                            min_distance = temp_total;
                            pi = i;
                            qj = j;
                            shorter_dist = true;
                        }
                    } else {
                        temp_total += distance_i;
                        if (temp_total < min_distance) {
                            min_distance = temp_total;
                            pi = i;
                            qj = j;
                            shorter_dist = false;
                        }
                    }
                }
            }
            if (min_distance < overall_min_distance) {
                overall_min_distance = min_distance;
                // remember indexes of the smallest
                min_pi1 = pi;
                min_qj1 = qj;
            }
        }
        if (check_intersection(points[g][path_indices[g][min_pi1]].first, points[g][path_indices[g][min_pi1]].second,
                               points[g + 1][path_indices[g + 1][min_qj1]].first, points[g + 1][path_indices[g + 1][min_qj1]].second,
                               points[g][path_indices[g][min_pi1 + 1]].first, points[g][path_indices[g][min_pi1 + 1]].second,
                               points[g + 1][path_indices[g + 1][min_qj1 + 1]].first, points[g + 1][path_indices[g + 1][min_qj1 + 1]].second)) {

            draw_line(colormap, blue, color, points[g][path_indices[g][min_pi1 + 1]].first, points[g][path_indices[g][min_pi1 + 1]].second,
                       points[g + 1][path_indices[g + 1][min_qj1]].first, points[g + 1][path_indices[g + 1][min_qj1]].second);
            draw_line(colormap, blue, color, points[g][path_indices[g][min_pi1]].first, points[g][path_indices[g][min_pi1]].second,
                      points[g + 1][path_indices[g + 1][min_qj1 + 1]].first, points[g + 1][path_indices[g + 1][min_qj1 + 1]].second);
            XFlush(dis);
        } else {
            if (shorter_dist) {
                draw_line(colormap, blue, color, points[g][path_indices[g][min_pi1 + 1]].first, points[g][path_indices[g][min_pi1 + 1]].second,
                          points[g + 1][path_indices[g + 1][min_qj1]].first, points[g + 1][path_indices[g + 1][min_qj1]].second);
                draw_line(colormap, blue, color, points[g][path_indices[g][min_pi1]].first, points[g][path_indices[g][min_pi1]].second,
                          points[g + 1][path_indices[g + 1][min_qj1 + 1]].first, points[g + 1][path_indices[g + 1][min_qj1 + 1]].second);
                XFlush(dis);
            } else {
                draw_line(colormap, blue, color, points[g][path_indices[g][min_pi1]].first, points[g][path_indices[g][min_pi1]].second,
                          points[g + 1][path_indices[g + 1][min_qj1]].first, points[g + 1][path_indices[g + 1][min_qj1]].second);
                draw_line(colormap, blue, color, points[g][path_indices[g][min_pi1 + 1]].first, points[g][path_indices[g][min_pi1 + 1]].second,
                          points[g + 1][path_indices[g + 1][min_qj1 + 1]].first, points[g + 1][path_indices[g + 1][min_qj1 + 1]].second);
                XFlush(dis);
            }
        }
        skip_index = min_qj1;
        connecting_indices.push_back(make_pair(min_pi1, min_qj1));
        total = overall_min_distance;
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

// X routines
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
}
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
    cout << "Program was terminated successfully\n";
    exit(1);
}