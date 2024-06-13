/**
 * @authors @UN-L, @Agentcastor and @TheKingHydra
 * Report made in collaboration: @UN-L - 24R51513, @Agentcastor - 24R51516 and @TheKingHydra - 24R51515
*/
#include <cstdio>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
using namespace std;

int main (int argc, char** argv) {
    int nx = 41;
    int ny = 41;
    int nt = 500;
    int nit = 50;
    double dx = 2.0 / (nx - 1);
    double dy = 2.0 / (ny - 1);
    float dt = .01;
    int rho = 1;
    float nu = .02;
    
    std::vector<std::vector<double>> u(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> v(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> p(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> b(ny, std::vector<double>(nx, 0.0));
    ofstream ufile("u.dat");
    ofstream vfile("v.dat");
    ofstream pfile("p.dat");
    for (int n = 0; n < nt; ++n) 
    {
        for (int j = 1; j < ny - 1; ++j) 
        {
            for (int i = 1; i < nx - 1; ++i) 
            {
                b[j][i] = rho * (1 / dt *
                (((u[j][i + 1] - u[j][i - 1]) / (2 * dx)) + ((v[j + 1][i] - v[j - 1][i]) / (2 * dy))) -
                std::pow(((u[j][i + 1] - u[j][i - 1]) / (2 * dx)), 2) - 
                2 * (((u[j + 1][i] - u[j - 1][i]) / (2 * dy)) * ((v[j][i + 1] - v[j][i - 1]) / (2 * dx))) - 
                std::pow(((v[j + 1][i] - v[j - 1][i]) / (2 * dy)), 2));
            }
        }
        for (int it = 0; it < nit; ++it) {
            std::vector<std::vector<double>> pn = p;
            for (int j = 1; j < ny - 1; ++j) {
                for (int i = 1; i < nx - 1; ++i) {
                    p[j][i] = (dy * dy * (pn[j][i + 1] + pn[j][i - 1]) +
                               dx * dx * (pn[j + 1][i] + pn[j - 1][i]) -
                               b[j][i] * dx * dx * dy * dy) /
                              (2 * (dx * dx + dy * dy));
                }
            }
            
            //p[:, -1] = p[:, -2]
            for (int j = 0; j < ny; ++j) p[j][nx - 1] = p[j][nx - 2];
            //p[0, :] = p[1, :]
            for (int i = 0; i < nx; ++i) p[0][i] = p[1][i];
            //p[:, 0] = p[:, 1]
            for (int j = 0; j < ny; ++j) p[j][0] = p[j][1];
            //p[-1, :] = 0
            for (int i = 0; i < nx; ++i) p[ny - 1][i] = 0.0;
        }

        std::vector<std::vector<double>> un = u;
        std::vector<std::vector<double>> vn = v;
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                u[j][i] = un[j][i] - un[j][i] * dt / dx * (un[j][i] - un[j][i - 1]) -
                        un[j][i] * dt / dy * (un[j][i] - un[j - 1][i]) -
                        dt / (2 * rho * dx) * (p[j][i + 1] - p[j][i - 1]) +
                        nu * dt / (dx * dx) * (un[j][i + 1] - 2 * un[j][i] + un[j][i - 1]) +
                        nu * dt / (dy * dy) * (un[j + 1][i] - 2 * un[j][i] + un[j - 1][i]);
                v[j][i] = vn[j][i] - vn[j][i] * dt / dx * (vn[j][i] - vn[j][i - 1]) -
                        vn[j][i] * dt / dy * (vn[j][i] - vn[j - 1][i]) -
                        dt / (2 * rho * dx) * (p[j + 1][i] - p[j - 1][i]) +
                        nu * dt / (dx * dx) * (vn[j][i + 1] - 2 * vn[j][i] + vn[j][i - 1]) +
                        nu * dt / (dy * dy) * (vn[j + 1][i] - 2 * vn[j][i] + vn[j - 1][i]);
            }
        }
        //u[0, :]  = 0
        //u[-1, :] = 1
        //v[0, :]  = 0 
        //v[-1, :] = 0
        for (int i = 0; i < nx; ++i)
        {
            u[0][i] = 0;
            u[nx-1][i] = 1;
            v[0][i] = 0;
            v[nx-1][i] = 0;
        }
        //u[:, 0]  = 0
        //u[:, -1] = 0
        //v[:, 0]  = 0
        //v[:, -1] = 0
        for (int j = 0; j < ny; ++j)
        {
            u[j][0] = 0;
            u[j][ny-1] = 0;
            v[j][0] = 0;
            v[j][ny-1] = 0;
        }
        if (n % 10 == 0) {
        for (int j=0; j<ny; j++)
            for (int i=0; i<nx; i++)
            ufile << u[j][i] << " ";
        ufile << "\n";
        for (int j=0; j<ny; j++)
            for (int i=0; i<nx; i++)
            vfile << v[j][i] << " ";
        vfile << "\n";
        for (int j=0; j<ny; j++)
            for (int i=0; i<nx; i++)
            pfile << p[j][i] << " ";
        pfile << "\n";
        }
    }
    ufile.close();
    vfile.close();
    pfile.close();
    return 0;
}
