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
/*
 * https://www.youtube.com/watch?v=-JjA4BLQyqE
 * https://www.codeproject.com/articles/762581/held-karp-algorithm-implementation-in-csharp
 * http://gebweb.net/blogpost/2011/06/24/the-dynamic-programming-algorithm-for-the-travelling-salesman-problem/
 * */

// vars
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
void held_karp(unsigned long, int);
vector<vector<int> > gen_combinations(int, int);
int dist(int, int, int, int);
bool compare_pair(pair<int, int> &, pair<int, int> &);
void connect_slices();
void draw_point(Colormap, char [], XColor &, int);
void draw_point(Colormap, char [], XColor &, int, int);
void draw_line(Colormap, char [], XColor &, int, int, vector<pair<int, int> >);
void draw_line(Colormap c, char [], XColor &, int, int, int, int);
void draw_line(Colormap c, char [], XColor &, int, int);



int main(int argc, char *argv[]) {
    //std::cout << "Hello, World!" << std::endl;
    if (argc == 2) {
        printf("Mode: Command-Line\n");
        int x, y;
        ifstream input_file(argv[1]);
        while(input_file >> x >> y) {
            coordinates.push_back(make_pair(x, y));
        }
        input_file.close();

        /*cout << "Number of coordinates: " << coordinates.size() << endl;

        for (int i = 0; i < coordinates.size(); ++i) {
            printf("x = %d, y = %d \n", coordinates[i].first, coordinates[i].second);
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
    char navy_blue[] = "#000080";
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
                    /*for (int i = 0; i < coordinates.size(); ++i) {
                        printf("x = %d, y = %d \n", coordinates[i].first, coordinates[i].second);
                    }*/

                    if (n <= 15) {
                        adj_Matrix.resize(n, vector<int>(n, DEFAULT_VAL));
                        //printf("adj_Matrix[%d][%d]\n", (int) adj_Matrix.size(), (int) adj_Matrix[0].size());

                        // calc distances and fill the matrix
                        for (int i = 0; i < n; ++i) {
                            for (int j = i + 1; j < n; ++j) {
                                adj_Matrix[i][j] = adj_Matrix[j][i] = dist(coordinates[i].first, coordinates[i].second,
                                                                           coordinates[j].first, coordinates[j].second);
                            }
                        }

                        /*for (int i = 0; i < n; ++i) {
                            for (int j = 0; j < n; ++j) {
                                printf("%d,",adj_Matrix[i][j]);
                            }
                            cout << "\n";
                        }*/
                        held_karp(n, 1);
                        //cout << "path_indices.size(): " << path_indices.size() << endl;
                        for (int k = 0; k < path_indices.size(); ++k) {
                            for (int i = 0; i < path_indices[k].size() - 1; ++i) {
                                draw_line(colormap, firebrick, color, path_indices[k][i], path_indices[k][i + 1]);

                            }
                        }
                        display_result = false;
                    } else {
                        // first sort them
                        sort(coordinates.begin(), coordinates.end(), compare_pair);
                        /*for (int i = 0; i < coordinates.size(); ++i) {
                            printf("x = %d, y = %d \n", coordinates[i].first, coordinates[i].second);
                        }*/
                        // break them in sets of 15

                        int k = (int) (n / 15);
                        int dif = (int) (n - (k * 15));
                        int div = k;
                        int end = 0;
                        unsigned long points_size = 0;
                        //unsigned long path_size = 0;
                        //cout << k << endl;
                        //cout << dif << endl;

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
                            //cout << "range " << i * 15 << " - " << end << " size: " << points_size <<"\n";
                            vector<pair<int, int> > points_slice;//(coordinates.begin() + (i * 15), coordinates.end() - end);

                            copy(coordinates.begin() + (i * 15), coordinates.end() - end, back_inserter(points_slice));

                            //printf("points_slice size: %d\n", (int) points_slice.size());


                            /*for (int l = 0; l < points_slice.size(); ++l) {
                                printf("x = %d, y = %d \n", points_slice[l].first, points_slice[l].second);
                            }*/

                            adj_Matrix.resize(points_size, vector<int>(points_size, DEFAULT_VAL));
                            //printf("adj_Matrix[%d][%d]\n", (int) adj_Matrix.size(), (int) adj_Matrix[0].size());
                            //path_size += points_size;
                            //cout << "path_size: " << path_size << endl;
                            //path_indices.resize(path_size);


                            /*for (int w1 = 0; w1 < points_size; ++w1) {
                                for (int j = 0; j < points_size; ++j) {
                                    printf("%d,",adj_Matrix[w1][j]);
                                }
                                cout << "\n";
                            }*/

                            // calc distances and fill the matrix
                            for (int indx = 0; indx < points_size; ++indx) {
                                for (int j = indx + 1; j < points_size; ++j) {
                                    adj_Matrix[indx][j] = adj_Matrix[j][indx] = dist(points_slice[indx].first, points_slice[indx].second,
                                                                                     points_slice[j].first, points_slice[j].second);
                                }
                            }

                            held_karp(points_size, i + 1);

                            points.push_back(points_slice);

                            /*for (int w1 = 0; w1 < points_size; ++w1) {
                                for (int j = 0; j < points_size; ++j) {
                                    printf("%d,",adj_Matrix[w1][j]);
                                }
                                cout << "\n";
                            }*/
                            //path_indices.clear();
                            points_slice.clear();
                            adj_Matrix.clear();


                        }

                        connect_slices();

                        /*cout << "points.size(): " << path_indices.size() << endl;
                        for (int p = 0; p < points.size(); ++p) {
                            for (int i = 0; i < points[p].size(); ++i) {
                                cout << i << ": x " << points[p][i].first << " , y " << points[p][i].second  << endl;
                            }
                        }*/

                        /*for (int l = 0; l < connecting_indices.size(); ++l) {
                            printf("p = %d, q = %d \n", connecting_indices[l].first, connecting_indices[l].second);
                        }*/


                        int num = (int) connecting_indices.size();
                        //bool check = true;
                        //cout << "path_indices.size(): " << path_indices.size() - 1 << endl;
                        for (int m = 0; m < path_indices.size(); ++m) {
                            //cout << m << ": " << path_indices[m] << endl;
                            for (int i = 0; i < path_indices[m].size() - 1; ++i) {
                                /*cout << i << ": x " << points[m][path_indices[m][i]].first << ", y "
                                                    << points[m][path_indices[m][i + 1]].second << endl;*/

                                if (m == 0) {
                                    if (!(i == 0 && i == connecting_indices[m].first)) {
                                        draw_line(colormap, firebrick, color, path_indices[m][i], path_indices[m][i + 1], points[m]);
                                    }
                                    if (i + 1 == connecting_indices[m].first) {
                                        //continue;
                                        for (int j = 0; j < path_indices[m + 1].size(); ++j) {
                                            if (j == connecting_indices[m].second) {
                                                draw_line(colormap, firebrick, color,
                                                          points[m][path_indices[m][i + 1]].first,
                                                          points[m][path_indices[m][i + 1]].second,
                                                          points[m + 1][path_indices[m + 1][j]].first,
                                                          points[m + 1][path_indices[m + 1][j]].second);

                                                draw_line(colormap, firebrick, color,
                                                          points[m][path_indices[m][i + 2]].first,
                                                          points[m][path_indices[m][i + 2]].second,
                                                          points[m + 1][path_indices[m + 1][j + 1]].first,
                                                          points[m + 1][path_indices[m + 1][j + 1]].second);
                                            }
                                        }
                                        //cout << "add\n";
                                        i+=1;
                                    }
                                    //cout << i << endl;
                                } else {
                                    if (m != num) {
                                        if (!(i == 0 && i == connecting_indices[m].first ||
                                              i == 0 && i == connecting_indices[m - 1].second)) {
                                            draw_line(colormap, firebrick, color, path_indices[m][i], path_indices[m][i + 1], points[m]);
                                        }
                                        if (i + 1 == connecting_indices[m - 1].second) {
                                            i+=1;
                                            //cout << "add\n";
                                        }
                                        if (i + 1 == connecting_indices[m].first) {
                                            //continue;
                                            for (int j = 0; j < path_indices[m + 1].size(); ++j) {
                                                if (j == connecting_indices[m].second) {
                                                    draw_line(colormap, firebrick, color,
                                                              points[m][path_indices[m][i + 1]].first,
                                                              points[m][path_indices[m][i + 1]].second,
                                                              points[m + 1][path_indices[m + 1][j]].first,
                                                              points[m + 1][path_indices[m + 1][j]].second);

                                                    draw_line(colormap, firebrick, color,
                                                              points[m][path_indices[m][i + 2]].first,
                                                              points[m][path_indices[m][i + 2]].second,
                                                              points[m + 1][path_indices[m + 1][j + 1]].first,
                                                              points[m + 1][path_indices[m + 1][j + 1]].second);
                                                }
                                            }
                                            i+=1;
                                            //cout << "add\n";
                                        }
                                        //cout << i << endl;
                                    }
                                    if (num == m) {
                                        if (!(i == 0 && i == connecting_indices[m - 1].second)) {
                                            draw_line(colormap, firebrick, color, path_indices[m][i], path_indices[m][i + 1], points[m]);
                                        }
                                        if (i + 1 == connecting_indices[m - 1].second) {
                                            i+=1;
                                            //cout << "add\n";
                                        }
                                        //cout << i << endl;
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

void held_karp(unsigned long n, int group) {
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

    vector<vector<int> > combinations;

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


            // Find the lowest cost to get to this subset

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

    // Backtrack to find full path
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
    // add manual starting and ending coordinates.
    path.push_back(0);
    path.push_front(0);

    path.reverse();

    //int *path_index = (int *) malloc(n + 1);
    //int path_index[n + 1];
    //int index = 0;

    //cout << "opt: " << opt << " ( ";

    vector<int> path_slice;

    for(std::list<int>::iterator list_iter = path.begin();
        list_iter != path.end(); list_iter++)
    {
        //std::cout<<*list_iter << " ";
        //path_index[index] = *list_iter;
        path_slice.push_back(*list_iter);
        //index++;
    }
    //cout << "\n";


    /*for (int l = 0; l < n + 1; ++l) {
        cout << l << ": " << path_slice[l] << endl;
    }*/



    cout << "Total length of the tour #" << group << ": "  << opt << endl;
    //printf("Total length of the tour #%d: %d", group, opt );
    path_indices.push_back(path_slice);

    // clear containers
    //path_slice.clear();
    combinations.clear();
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

bool compare_pair(pair<int, int>& i, pair<int, int>& j) {
    return i.first < j.first;
}

void connect_slices() {
    // finds optimal connections between the groups
    for (int g = 0; g < path_indices.size() - 1; ++g) { // number of groups
        // pick an edge from group i
        int pi = 0;
        //int pi1 = 0;
        int qj = 0;
        //int qj1 = 0;

        int overall_min_distance = 99999;
        int min_distance = 99999;
        for (int i = 0; i < path_indices[g].size() - 1; ++i) {
            int pix =  points[g][path_indices[g][i]].first,     piy = points[g][path_indices[g][i]].second;
            int pi1x = points[g][path_indices[g][i + 1]].first, pi1y = points[g][path_indices[g][i + 1]].second;
            //printf("pix=%d, piy=%d, pi1x=%d, pi1y=%d\n", pix, pix, pi1x, pi1y);
            for (int j = 0; j < path_indices[g + 1].size() - 1; ++j) {
                int qjx =  points[g + 1][path_indices[g + 1][j]].first,     qjy = points[g + 1][path_indices[g + 1][j]].second;
                int qj1x = points[g + 1][path_indices[g + 1][j + 1]].first, qj1y = points[g + 1][path_indices[g + 1][j + 1]].second;

                //int dist_pi_qj1 = dist(pix, piy, qj1x, qj1y);
                //int dist_pi1_qj = dist(pi1x, pi1y, qjx, qjy);
                int dist_pi_qj = dist(pix, piy, qjx, qjy);
                int dist_pi1_qj1 = dist(pi1x, pi1y, qj1x, qj1y);

                int distance =  dist_pi_qj + dist_pi1_qj1; //dist_pi_qj1 + dist_pi1_qj +

                if (distance < min_distance) {
                    min_distance = distance;
                    pi = i;
                    //pi1 = i + 1;
                    qj = j;
                    //qj1 = j + 1;
                    //cout << "MIN Distance: " << distance << endl;
                }
                //cout << "Distance: " << distance << endl;
            }
            //cout << "Sub MIN Distance: " << min_distance << endl;
            if (min_distance < overall_min_distance) {
                overall_min_distance = min_distance;
            }
            //cout << endl;
        }
        //cout << "----------------------OVERALL MIN Distance: " << overall_min_distance << " ----------------------" << endl;
        //printf("pi=%d, qj=%d\n", pi, qj);
        connecting_indices.push_back(make_pair(pi, qj));
    }
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
    exit(1);
};

void redraw() {
    XClearWindow(dis, win);
}