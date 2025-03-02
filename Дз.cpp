#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

using namespace std;

struct Point {
    double x, y;
};

double distance(const Point& p1, const Point& p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

double calculate_path_length(const vector<int>& path, const vector<Point>& points) {
    double length = 0.0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        length += distance(points[path[i]], points[path[i + 1]]);
    }
    return length;
}

double calculate_cycle_length(const vector<int>& path, const vector<Point>& points) {
    double length = calculate_path_length(path, points);
    length += distance(points[path.back()], points[path.front()]);
    return length;
}


void two_opt_swap(vector<int>& path, int i, int k) {
    reverse(path.begin() + i + 1, path.begin() + k + 1);
}

pair<double, vector<int>> two_opt(vector<int> path, const vector<Point>& points, bool is_cycle) {
    double best_length = is_cycle ? calculate_cycle_length(path, points) : calculate_path_length(path, points);
    bool improved = true;

    while (improved) {
        improved = false;
        for (size_t i = 0; i < path.size() - 1; ++i) {
            for (size_t k = i + 2; k < path.size(); ++k) {
                vector<int> new_path = path;
                two_opt_swap(new_path, i, k);
                double new_length = is_cycle ? calculate_cycle_length(new_path, points) : calculate_path_length(new_path, points);

                if (new_length < best_length) {
                    best_length = new_length;
                    path = new_path;
                    improved = true;
                }
            }
        }
    }
    return { best_length, path };
}

void three_opt_swap(vector<int>& path, int i, int j, int k, int permutation) {
    vector<int> new_path;

    for (int l = 0; l <= i; ++l) {
        new_path.push_back(path[l]);
    }

    if (permutation == 0) {
        for (int l = i + 1; l <= j; ++l) {
            new_path.push_back(path[l]);
        }
        for (int l = j + 1; l <= k; ++l) {
            new_path.push_back(path[l]);
        }
    }
    else if (permutation == 1) {
        for (int l = i + 1; l <= j; ++l) {
            new_path.push_back(path[l]);
        }
        for (int l = k; l >= j + 1; --l) {
            new_path.push_back(path[l]);
        }
    }
    else if (permutation == 2) {
        for (int l = j; l >= i + 1; --l) {
            new_path.push_back(path[l]);
        }
        for (int l = j + 1; l <= k; ++l) {
            new_path.push_back(path[l]);
        }
    }
    else if (permutation == 3) {
        for (int l = j; l >= i + 1; --l) {
            new_path.push_back(path[l]);
        }
        for (int l = k; l >= j + 1; --l) {
            new_path.push_back(path[l]);
        }
    }

    for (size_t l = k + 1; l < path.size(); ++l)
    {
        new_path.push_back(path[l]);
    }
    path = new_path;
}

pair<double, vector<int>> three_opt(vector<int> path, const vector<Point>& points, bool is_cycle) {
    double best_length = is_cycle ? calculate_cycle_length(path, points) : calculate_path_length(path, points);
    bool improved = true;

    while (improved) {
        improved = false;
        for (size_t i = 0; i < path.size() - 2; ++i) {
            for (size_t j = i + 1; j < path.size() - 1; ++j) {
                for (size_t k = j + 1; k < path.size(); ++k) {
                    for (int permutation = 0; permutation < 4; ++permutation) {
                        vector<int> new_path = path;
                        three_opt_swap(new_path, i, j, k, permutation);

                        double new_length = is_cycle ? calculate_cycle_length(new_path, points) : calculate_path_length(new_path, points);

                        if (new_length < best_length) {
                            best_length = new_length;
                            path = new_path;
                            improved = true;
                        }
                    }
                }
            }
        }
    }
    return { best_length, path };
}


int main() {
    int n;
    cin >> n;

    vector<Point> points(n);
    for (int i = 0; i < n; ++i) {
        cin >> points[i].x >> points[i].y;
    }

    int is_cycle;
    cout << "Enter 0 for path, 1 for cycle: ";
    cin >> is_cycle;

    vector<int> initial_path(n);
    for (int i = 0; i < n; ++i) {
        initial_path[i] = i;
    }

    pair<double, vector<int>> two_opt_result = two_opt(initial_path, points, is_cycle == 1);

    cout << two_opt_result.first << endl;
    cout << is_cycle << endl;

    for (int i = 0; i < n; ++i) {
        cout << two_opt_result.second[i] << (i == n - 1 ? "" : " ");
    }
    cout << endl;

    pair<double, vector<int>> three_opt_result = three_opt(initial_path, points, is_cycle == 1);

    cout << three_opt_result.first << endl;
    cout << is_cycle << endl;

    for (int i = 0; i < n; ++i) {
        cout << three_opt_result.second[i] << (i == n - 1 ? "" : " ");
    }
    cout << endl;

    return 0;
}