#include <string>
#include <vector>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <cmath>
#include <iostream>
#include <tuple>
#include <unistd.h>

using namespace std;

/*
==================================================================================
* function name: fitness
* user defined fitness function
*    
* @para x: start x
* @para v: search direction
* @para lam: the lambda (ratio, scaler)
* return values: 
*  fintness function value
==================================================================================
*/
double fitness(vector<double> vx){
    double x, y, z;
    x = vx[0];
    y = vx[1];
    z = vx[2];
    double f = 1.0 / (1.0 + (x-y)*(x-y)) + sin(0.5 * y * z) + 
    exp(-((x+z)/y-2)*((x+z)/y-2));
    usleep(5000);
    return -f;
}

/*
==================================================================================
* function name: fitness_direction
* an interface to user defined fitness function, for one direction search
*    
* @para x: start x
* @para v: search direction
* @para lam: the lambda (ratio, scaler)
* return values: 
*  fintness function value
==================================================================================
*/
double fitness_direction(vector<double>& x, double lam, vector<double>& v){
    int n = x.size();
    vector<double> x_new(n, 0.0);
    for (int i = 0; i < n; ++i){
        x_new[i] = x[i] + lam * v[i];
    }
    return fitness(x_new);
}
/*
==================================================================================
* function name: bracket
* bracket the interval of minimum point
*    
* @para x: start x
* @para v: search direction
* @para x1, start point
* @para h: initial search increment used in 'bracket'
* return values: 
*  the interval, upper and lower 
==================================================================================
*/

tuple<double, double> bracket(vector<double>& x,
                              vector<double>& v,
                              double x1,
                              double h){
    double c = 1.618033989;
    double x2 = h + x1;
    double f1, f2;
    f1 = fitness_direction(x, x1, v);
    f2 = fitness_direction(x, x2, v);

    if(f2 > f1){
        h = -h;
        x2 = x1 + h;
        f2 = fitness_direction(x, x2, v);
        if(f2 > f1){
            return tuple<double, double>(x2, x1-h);
        }
    }

    for (int i = 0; i < 100; ++i){
        h = c * h;
        double x3 = x2 + h;
        double f3 = fitness_direction(x, x3, v);
        if(f3 > f2){
            return tuple<double, double>(x1, x3);
        }
        x1 = x2;
        x2 = x3;
        f1 = f2;
        f2 = f3;
    }
    cout << "Bracket did not find a mimimum" << endl;
    exit(0);
    return tuple<double, double>(0.0, 0.0);
    //return tuple<int, unsigned, string, int>(1, 2, "3", 4);
}
/*
==================================================================================
* function name: golden_section_search
* The golden section search is a technique for finding the extremum 
* (minimum or maximum) of a strictly unimodal function by successively narrowing 
* the range of values inside which the extremum is known to exist. 
* i.e., find lambda to min f(x+lam*v)
*    
* @para x: start x
* @para v: search direction
* @para a, b: search interval [a,b]
* @para tol: tolerente of variables, i.e., the final interval width
* return values: 
*   lambda and minimum value found 
==================================================================================
*/
tuple<double, double> golden_section_search(vector<double>& x,
                              vector<double>& v,
                              double a,
                              double b,
                              double tol = 1.0e-9){
    // compute the number of telescoping operations required to 
    // reduce h from |b − a| to an error tolerance
    int nIter = ceil(-2.078087*log(tol/abs(b-a)));
    double R = 0.618033989;   // golden ratio
    double C = 1.0 - R;
    // First telescoping
    double x1 = R * a + C * b;
    double x2 = C * a + R * b;
    double f1, f2;
    f1 = fitness_direction(x, x1, v);
    f2 = fitness_direction(x, x2, v);
    for (int i = 0; i < nIter; ++i){
        if (f1 > f2){
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = C * a + R * b;
            f2 = fitness_direction(x, x2, v);
        }
        else{
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = R * a + C * b;
            f1 = fitness_direction(x, x1, v);
        }
    }
    if(f1 < f2){
        return tuple<double, double>(x1, f1);
    }
    else{
        return tuple<double, double>(x2, f2);
    }
}

/*
==================================================================================
* function name: mse
* calculating Mean squared error
*    
* @para v1: the first vector
* @para v2: the second vector
* return values: 
*   Mean squared error 
==================================================================================
*/
double mse(vector<double> v1, vector<double> v2){
    double len = 0.0;
    for (int i = 0; i < v1.size(); ++i){
        len += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    return sqrt(len / v1.size());
}
/*
==================================================================================
* function name: min_powell
* Powell's method of minimizing user-supplied function
* without calculating its derivatives
*    
* @para x: starting point 
* @para h: initial search increment used in 'bracket'
* @para tolerate:
* @para maxit: maximum iterations
* return values: 
*   a set of oarameters which will carry out the minimum (local)    
==================================================================================
*/
vector< double > min_powell(vector< double >& x, // initial value
                  double h,
                  double tolerate,
                  int maxit) 
{
    int n = x.size();                 // Number of design variables
    vector<double> df(n, 0);          // Decreases of fitness stored here
    // direction vectors v stored here by rows
    vector<vector<double> > u(n, vector<double>(n, 0.0) );
    // set direction vectors
    for (int i = 0; i < n; ++i){
        u[i][i] = 1.0;
    }

    // main iteration loop
    for (int j = 0; j < maxit; j++) { 
        vector< double > x_old = x;
        double fitness_old = fitness(x_old);
        vector<double> fitness_dir_min(n+1, 0.0);
        fitness_dir_min[0] = fitness_old;
        for (int i = 0; i < n; ++i){
            vector<double> v = u[i];
            double a, b, s;
            tie(a, b) = bracket(x, v, 0.0, h);
            tie(s, fitness_dir_min[i+1]) = golden_section_search(x, v, a, b);
            for (int i = 0; i < n; ++i){
                x[i] = x[i] + s * v[i];
            }
        }
        for (int i = 0; i < n; ++i){
            df[i] = fitness_dir_min[i] - fitness_dir_min[i+1];
        }
        // Last line gloden section search in the cycle    
        vector<double> v(n);
        for (int i = 0; i < n; ++i){
            v[i] = x[i] - x_old[i];
        }
        double a, b, s, dummy;
        tie(a, b) = bracket(x, v, 0.0, h);
        tie(s, dummy) = golden_section_search(x, v, a, b);
        // dependence among search directions
        for (int i = 0; i < n; ++i){
            x[i] = x[i] + s * v[i];
        }
        // Check for convergence
        if(mse(x, x_old) < tolerate){
            cout << "found minimize value at " << j+1 << " step" << " with value: " << fitness(x) << endl;
            return x;
        }
        // Identify biggest decrease & update search directions
        int i_max = 0;
        for (int i = 1; i < n; ++i){
            if(df[i] > df[i_max]){
                i_max = i;
            }
        }
        for (int i = i_max; i < n-1; ++i){
            u[i] = u[i+1];
        }
        u[n-1] = v;
    }
    cout << "Powell did not converge" << endl;
    return vector< double > (n, 0.0);
}

int main(int argc, char const *argv[])
{
    vector<double> vx_min;
    vector<double> vx_init{0.11, 0.18, 0.1};
    vx_min = min_powell(vx_init, 0.1, 1.0e-6, 30);
    cout << vx_min[0] << ", " << vx_min[1] << ", " << vx_min[2] << endl;
    return 0;
}