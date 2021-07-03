#include<pybind11/pybind11.h>
#include<pybind11/stl.h>
#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<string>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using vec = std::vector<double>;


struct SplineSet{
    // def: S(u) = a + b(u-x) + c(u-x)^2) + d(u-x)^3
    // y2 = a + b(x2-x1)
    // y1 = a + b(x1-x1)
    double a;
    double b;
    double c;
    double d;
    double x;
};


class SplineInterpolate
{
private:
    std::vector<SplineSet> spline_sets;
    SplineSet left_extra_set{0, 0, 0, 0, 0};
    SplineSet right_extra_set{0, 0, 0, 0, 0};
    std::string extrap_method;
    vec knobs;
    
    std::vector<SplineSet> linear_spline(vec &x, vec &y);
    std::vector<SplineSet> cubic_spline(vec &x, vec &y);
    double calc_y_within_spline(SplineSet &splineset, double x) {
        double a = splineset.a;
        double b = splineset.b;
        double c = splineset.c;
        double d = splineset.d;
        double u = x - splineset.x;

        return a + b * u + c * pow(u, 2) + d * pow(u, 3);
    }
    
public:
    SplineInterpolate(std::string &method, std::string &extraploate_method, vec &x, vec &y){
        int n = x.size() - 1;

        if (method == "linear") {
            spline_sets = linear_spline(x, y);
        }
        else if (method == "cubic") {
            spline_sets = cubic_spline(x, y);
        }

        if (extraploate_method == "linear") {
            left_extra_set.a = y[0];
            left_extra_set.b = (y[1] - y[0]) / (x[1] - x[0]);
            left_extra_set.x = x[0];

            right_extra_set.a = y[n-1];
            right_extra_set.b = (y[n] - y[n-1]) / (x[n] - x[n-1]);
            right_extra_set.x = x[n-1];
        }
        else if (extraploate_method == "flat") {
            left_extra_set.a = y[0];
            right_extra_set.a = y[n];
        }
        
        for (int i = 0; i < x.size(); i++) {
            knobs.push_back(x[i]);
        }
    }

    double F(double x){
        int pos = std::lower_bound(knobs.begin(), knobs.end(), x) - knobs.begin();
        int n = knobs.size() - 1;
        
        if (pos == 0)
            return calc_y_within_spline(left_extra_set, x);
        else if (pos == n+1)
            return calc_y_within_spline(right_extra_set, x);
        else {
            return calc_y_within_spline(spline_sets[pos-1], x);
        }
    }
};

std::vector<SplineSet> SplineInterpolate::linear_spline(vec &x, vec &y)
{
    int n = x.size()-1;
    vec a, b;
    a.insert(a.begin(), y.begin(), y.end());

    for (int i = 0; i < n; i++) {
        double b_i = (y[i+1] - y[i]) / (x[i+1] - x[i]);
        b.push_back(b_i);
    }
    
    std::vector<SplineSet> output_set(n);
    
    for(int i = 0; i < n; ++i)
    {
        output_set[i].a = a[i];
        output_set[i].b = b[i];
        output_set[i].c = 0;
        output_set[i].d = 0;
        output_set[i].x = x[i];
    }
    return output_set;
}

std::vector<SplineSet> SplineInterpolate::cubic_spline(vec &x, vec &y)
{
    // reference: https://stackoverflow.com/a/19216702
    int n = x.size()-1;
    vec a;
    a.insert(a.begin(), y.begin(), y.end());
    vec b(n);
    vec d(n);
    vec h;

    for(int i = 0; i < n; ++i)
        h.push_back(x[i+1]-x[i]);

    vec alpha;
    alpha.push_back(0);
    for(int i = 1; i < n; ++i)
        alpha.push_back( 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1]  );

    vec c(n+1);
    vec l(n+1);
    vec mu(n+1);
    vec z(n+1);
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for(int i = 1; i < n; ++i)
    {
        l[i] = 2 *(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    for(int j = n-1; j >= 0; --j)
    {
        c[j] = z [j] - mu[j] * c[j+1];
        b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
        d[j] = (c[j+1]-c[j])/3/h[j];
    }

    std::vector<SplineSet> output_set(n);
    
    for(int i = 0; i < n; ++i)
    {
        output_set[i].a = a[i];
        output_set[i].b = b[i];
        output_set[i].c = c[i];
        output_set[i].d = d[i];
        output_set[i].x = x[i];
    }
    return output_set;
}

namespace py = pybind11;

PYBIND11_MODULE(cpp_spline_interp, m) {
    py::class_<SplineInterpolate>(m, "SplineInterpolate")
        .def(py::init<std::string &, std::string &, vec &, vec &>())
        .def("F", &SplineInterpolate::F);


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
