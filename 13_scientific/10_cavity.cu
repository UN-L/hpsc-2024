/**
 * @authors @UN-L, @Agentcastor and @TheKingHydra
 * Report made in collaboration: @UN-L - 24R51513, @Agentcastor - 24R51516 and @TheKingHydra - 24R51515
*/
#include <iostream>
#include <fstream>
using namespace std;



__global__ void zeros(double *array, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int index = i * ny + j;
    if (i < nx && j < ny) {
        array[index] = 0;
    }
}

__global__ void compute_b(double *b, double *u, double *v, int nx, int ny, int rho, double dt, double dx, double dy) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= 1 && i < nx - 1 && j >= 1 && j < ny - 1) {
        int index = j * nx + i;
        int idx_i_plus_one = j * nx + (i+1);
        int idx_i_minus_one = j * nx + (i-1);
        int idx_j_plus_one = (j+1) * nx + i;
        int idx_j_minus_one = (j-1) * nx + i;

        b[index] = rho * (1 / dt *
            (((u[idx_i_plus_one] - u[idx_i_minus_one]) / (2 * dx)) + ((v[idx_j_plus_one] - v[idx_j_minus_one]) / (2 * dy))) -
            pow(((u[idx_i_plus_one] - u[idx_i_minus_one]) / (2 * dx)), 2) - 
            2 * (((u[idx_j_plus_one] - u[idx_j_minus_one]) / (2 * dy)) * ((v[idx_i_plus_one] - v[idx_i_minus_one]) / (2 * dx))) - 
            pow(((v[idx_j_plus_one] - v[idx_j_minus_one]) / (2 * dy)), 2));
    }
}

__global__ void compute_p(double *p, double *b, double *pn, int nx, int ny, double dx, double dy) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= 1 && i < nx - 1 && j >= 1 && j < ny - 1) {
        int index = j * nx + i;
        int idx_i_plus_one = j * nx + (i+1);
        int idx_i_minus_one = j * nx + (i-1);
        int idx_j_plus_one = (j+1) * nx + i;
        int idx_j_minus_one = (j-1) * nx + i;
        p[index] = (dy * dy * (pn[idx_i_plus_one] + pn[idx_i_minus_one]) +
                    dx * dx * (pn[idx_j_plus_one] + pn[idx_j_minus_one]) -
                    b[index] * dx * dx * dy * dy) /
                    (2 * (dx * dx + dy * dy));
    }
}

__global__ void compute_uv(double *u, double *v, double *p, double *un, double *vn, int nx, int ny, int rho, double dt, double nu, double dx, double dy) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= 1 && i < nx - 1 && j >= 1 && j < ny - 1) {
        int index = j * nx + i;
        int idx_i_plus_one = j * nx + (i+1);
        int idx_i_minus_one = j * nx + (i-1);
        int idx_j_plus_one = (j+1) * nx + i;
        int idx_j_minus_one = (j-1) * nx + i;
        u[index] = un[index] - un[index] * dt / dx * (un[index] - un[idx_i_minus_one]) -
                un[index] * dt / dy * (un[index] - un[idx_j_minus_one]) -
                dt / (2 * rho * dx) * (p[idx_i_plus_one] - p[idx_i_minus_one]) +
                nu * dt / (dx * dx) * (un[idx_i_plus_one] - 2 * un[index] + un[idx_i_minus_one]) +
                nu * dt / (dy * dy) * (un[idx_j_plus_one] - 2 * un[index] + un[idx_j_minus_one]);
        v[index] = vn[index] - vn[index] * dt / dx * (vn[index] - vn[idx_i_minus_one]) -
                vn[index] * dt / dy * (vn[index] - vn[idx_j_minus_one]) -
                dt / (2 * rho * dx) * (p[idx_j_plus_one] - p[idx_j_minus_one]) +
                nu * dt / (dx * dx) * (vn[idx_i_plus_one] - 2 * vn[index] + vn[idx_i_minus_one]) +
                nu * dt / (dy * dy) * (vn[idx_j_plus_one] - 2 * vn[index] + vn[idx_j_minus_one]);
    }
}

__global__ void compute_pj(double *p, int nx, int ny) {
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (j < ny) {
        int idx_minus_1 = j * nx + (nx - 1);
        int idx_minus_2 = j * nx + (nx - 2);
        int idx_0 = j * nx + 0;
        int idx_1 = j * nx + 1;
        p[idx_minus_1] = p[idx_minus_2];
        p[idx_0] = p[idx_1];
    }
}

__global__ void compute_pi(double *p, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nx) {
        int idx_minus_1 = (ny - 1) * nx + i;
        int idx_0 = 0 * nx + i;
        int idx_1 = 1 * nx + i;
        p[idx_0] = p[idx_1];
        p[idx_minus_1] = 0.0;
    }
}

__global__ void compute_uvi(double *u, double *v, int nx, int ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nx) {
        int idx_minus_1 = (ny - 1) * nx + i;
        int idx_0 = 0 * nx + i;
        u[idx_0] = 0;
        u[idx_minus_1] = 1;
        v[idx_0] = 0;
        v[idx_minus_1] = 0;
    }
}

__global__ void compute_uvj(double *u, double *v, int nx, int ny) {
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (j < ny) {
        int idx_minus_1 = j * nx + (nx - 1);
        int idx_0 = j * nx + 0;
        u[idx_0] = 0;
        u[idx_minus_1] = 0;
        v[idx_0] = 0;
        v[idx_minus_1] = 0;
    }
}

int main(int argc, char** argv) {
    int nx = 41;
    int ny = 41;
    int nt = 500;
    int nit = 50;
    double dx = 2.0 / (nx - 1);
    double dy = 2.0 / (ny - 1);
    float dt = .01;
    int rho = 1;
    float nu = .02;

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((nx + threadsPerBlock.x - 1) / threadsPerBlock.x, (ny + threadsPerBlock.y - 1) / threadsPerBlock.y);

    ofstream ufile("u.dat");
    ofstream vfile("v.dat");
    ofstream pfile("p.dat");

    double *u;
    cudaMallocManaged(&u, ny*nx*sizeof(double));
    zeros<<<numBlocks, threadsPerBlock>>>(u, nx, ny);
    cudaDeviceSynchronize();

    double *v;
    cudaMallocManaged(&v, ny*nx*sizeof(double));
    zeros<<<numBlocks, threadsPerBlock>>>(v, nx, ny);
    cudaDeviceSynchronize();

    double *p;
    cudaMallocManaged(&p, ny*nx*sizeof(double));
    zeros<<<numBlocks, threadsPerBlock>>>(p, nx, ny);
    cudaDeviceSynchronize();

    double *b;
    cudaMallocManaged(&b, ny*nx*sizeof(double));
    zeros<<<numBlocks, threadsPerBlock>>>(b, nx, ny);
    cudaDeviceSynchronize();

    double *pn;
    cudaMallocManaged(&pn, ny*nx*sizeof(double));

    double *un;
    cudaMallocManaged(&un, ny*nx*sizeof(double));

    double *vn;
    cudaMallocManaged(&vn, ny*nx*sizeof(double));

    for (int n = 0; n < nt; ++n) {
        compute_b<<<numBlocks, threadsPerBlock>>>(b, u, v, nx, ny, rho, dt, dx, dy);
        cudaDeviceSynchronize();

        for (int it = 0; it < nit; ++it) {
            cudaMemcpy(pn, p, nx * ny * sizeof(double), cudaMemcpyDeviceToDevice);
            compute_p<<<numBlocks, threadsPerBlock>>>(p, b, pn, nx, ny, dx, dy);
            cudaDeviceSynchronize();
            compute_pj<<<numBlocks, threadsPerBlock>>>(p, nx, ny);
            compute_pi<<<numBlocks, threadsPerBlock>>>(p, nx, ny);
            cudaDeviceSynchronize();
        }

        cudaMemcpy(un, u, nx * ny * sizeof(double), cudaMemcpyDeviceToDevice);
        cudaMemcpy(vn, v, nx * ny * sizeof(double), cudaMemcpyDeviceToDevice);
        compute_uv<<<numBlocks, threadsPerBlock>>>(u, v, p, un, vn, nx, ny, rho, dt, nu, dx, dy);
        cudaDeviceSynchronize();

        compute_uvi<<<numBlocks, threadsPerBlock>>>(u, v, nx, ny);
        cudaDeviceSynchronize();
        compute_uvj<<<numBlocks, threadsPerBlock>>>(u, v, nx, ny);
        cudaDeviceSynchronize();
        if (n % 10 == 0) {
            for (int j=0; j<ny; j++)
                for (int i=0; i<nx; i++)
                ufile << u[j * ny + i] << " ";
            ufile << "\n";
            for (int j=0; j<ny; j++)
                for (int i=0; i<nx; i++)
                vfile << v[j * ny + i] << " ";
            vfile << "\n";
            for (int j=0; j<ny; j++)
                for (int i=0; i<nx; i++)
                pfile << p[j * ny + i] << " ";
            pfile << "\n";
        }
    }
   
    cudaFree(u);
    cudaFree(v);
    cudaFree(p);
    cudaFree(b);
    cudaFree(pn);
    cudaFree(un);
    cudaFree(vn);
    ufile.close();
    vfile.close();
    pfile.close();
    return 0;
}