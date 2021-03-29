#include "core/Simulation.h"


namespace boltzmann {
    namespace core {
        __global__
        void collide(uint32_t xdim, uint32_t ydim,
                     bool **barrier,
                     double **n0,
                     double **nN,
                     double **nS,
                     double **nE,
                     double **nW,
                     double **nNW,
                     double **nNE,
                     double **nSW,
                     double **nSE,
                     double **density,
                     double **xvel,
                     double **yvel,
                     double **speed2,
                     double **n0_temp,
                     double **nN_temp,
                     double **nS_temp,
                     double **nE_temp,
                     double **nW_temp,
                     double **nNW_temp,
                     double **nNE_temp,
                     double **nSW_temp,
                     double **nSE_temp,
                     double **density_temp,
                     double **xvel_temp,
                     double **yvel_temp,
                     double **speed2_temp,
                     double omega) {
            uint32_t y = blockIdx.x;
            uint32_t x = threadIdx.x;

            const double four9ths = 4.0 / 9;
            const double one9th = 1.0 / 9;
            const double one36th = 1.0 / 36;

            if (y < ydim && x < xdim) {
                double n, one9thn, one36thn, vx, vy, vx2, vy2, vx3, vy3, vxvy2, v2, v215;
                if (!barrier[y][x]) {
                    n = n0[y][x] + nN[y][x] + nS[y][x] + nE[y][x] + nW[y][x] + nNW[y][x] + nNE[y][x] + nSW[y][x] +
                        nSE[y][x];
                    one9thn = 1.0 / 9.0 * n;
                    one36thn = 1.0 / 36.0 * n;
                    if (n > 0) {
                        vx = (nE[y][x] + nNE[y][x] + nSE[y][x] - nW[y][x] - nNW[y][x] - nSW[y][x]) / n;
                        vy = (nN[y][x] + nNE[y][x] + nNW[y][x] - nS[y][x] - nSE[y][x] - nSW[y][x]) / n;
                    } else {
                        vx = 0;
                        vy = 0;
                    }
                    vx3 = 3 * vx;
                    vy3 = 3 * vy;
                    vx2 = vx * vx;
                    vy2 = vy * vy;
                    vxvy2 = 2 * vx * vy;
                    v2 = vx2 + vy2;
                    v215 = 1.5 * v2;

                    density[y][x] = n;
                    xvel[y][x] = vx;
                    yvel[y][x] = vy;
                    speed2[y][x] = v2;
                    n0[y][x] += omega * (four9ths * n * (1 - v215) - n0[y][x]);
                    nE[y][x] += omega * (one9thn * (1 + vx3 + 4.5 * vx2 - v215) - nE[y][x]);
                    nW[y][x] += omega * (one9thn * (1 - vx3 + 4.5 * vx2 - v215) - nW[y][x]);
                    nN[y][x] += omega * (one9thn * (1 + vy3 + 4.5 * vy2 - v215) - nN[y][x]);
                    nS[y][x] += omega * (one9thn * (1 - vy3 + 4.5 * vy2 - v215) - nS[y][x]);
                    nNE[y][x] += omega * (one36thn * (1 + vx3 + vy3 + 4.5 * (v2 + vxvy2) - v215) - nNE[y][x]);
                    nNW[y][x] += omega * (one36thn * (1 - vx3 + vy3 + 4.5 * (v2 - vxvy2) - v215) - nNW[y][x]);
                    nSE[y][x] += omega * (one36thn * (1 + vx3 - vy3 + 4.5 * (v2 - vxvy2) - v215) - nSE[y][x]);
                    nSW[y][x] += omega * (one36thn * (1 - vx3 - vy3 + 4.5 * (v2 + vxvy2) - v215) - nSW[y][x]);
                }
            }
        }


        __global__
        void compute_curl(uint32_t xdim, uint32_t ydim, double **curl, double **yvel, double **xvel) {
            uint32_t y = blockIdx.x;
            uint32_t x = threadIdx.x;
            if (y >= 1 && y < ydim - 1 && x >= 1 && x < xdim - 1) {
                curl[y][x] = (yvel[y][x + 1] - yvel[y][x - 1]) - (xvel[y + 1][x] - xvel[y - 1][x]);
            }
        }

        __global__
        void stream(uint32_t xdim,
                    uint32_t ydim,
                    bool **barrier,
                    double **n0,
                    double **nN,
                    double **nS,
                    double **nE,
                    double **nW,
                    double **nNW,
                    double **nNE,
                    double **nSW,
                    double **nSE,
                    double **density,
                    double **xvel,
                    double **yvel,
                    double **speed2,
                    double **n0_temp,
                    double **nN_temp,
                    double **nS_temp,
                    double **nE_temp,
                    double **nW_temp,
                    double **nNW_temp,
                    double **nNE_temp,
                    double **nSW_temp,
                    double **nSE_temp,
                    double **density_temp,
                    double **xvel_temp,
                    double **yvel_temp,
                    double **speed2_temp,
                    double omega,
                    double v) {


            const double four9ths = 4.0 / 9;
            const double one9th = 1.0 / 9;
            const double one36th = 1.0 / 36;

            uint32_t y = blockIdx.x;
            uint32_t x = threadIdx.x;

            if (y > 0 && y <= ydim - 1 && x >= 0 && x < xdim - 1) {
                nN[y][x] = nN_temp[y - 1][x];
                nNW[y][x] = nNW_temp[y - 1][x + 1];
            }
            if (y > 0 && y <= ydim - 1 && x > 0 && x <= xdim - 1) {
                nE[y][x] = nE_temp[y][x - 1];
                nNE[y][x] = nNE_temp[y - 1][x - 1];
            }
            if (y >= 0 && y < ydim - 1 && x > 0 && x <= xdim - 1) {
                nS[y][x] = nS_temp[y + 1][x];
                nSE[y][x] = nSE_temp[y + 1][x - 1];
            }
            if (y >= 0 && y < ydim - 1 && x >= 0 && x < xdim - 1) {
                nW[y][x] = nW_temp[y][x + 1];
                nSW[y][x] = nSW_temp[y + 1][x + 1];
            }
            if (y >= 0 && y < ydim - 1 && x == 0) {

                nS[y][0] = nS_temp[y + 1][0];
            }
            if (y <= ydim - 1 && y > 0 && x == xdim - 1) {
                nN[y][xdim - 1] = nN_temp[y - 1][xdim - 1];
            }
            if (y >= 0 && y < ydim && x == 0 && !barrier[y][x]) {
                nE[y][0] = one9th * (1 + 3 * v + 3 * v * v);
                nNE[y][0] = one36th * (1 + 3 * v + 3 * v * v);
                nSE[y][0] = one36th * (1 + 3 * v + 3 * v * v);
            }
            if (y == 0 && x >= 0 && x < xdim) {
                n0[0][x] = four9ths * (1 - 1.5 * v * v);
                nE[0][x] = one9th * (1 + 3 * v + 3 * v * v);
                nW[0][x] = one9th * (1 - 3 * v + 3 * v * v);
                nN[0][x] = one9th * (1 - 1.5 * v * v);
                nS[0][x] = one9th * (1 - 1.5 * v * v);
                nNE[0][x] = one36th * (1 + 3 * v + 3 * v * v);
                nSE[0][x] = one36th * (1 + 3 * v + 3 * v * v);
                nNW[0][x] = one36th * (1 - 3 * v + 3 * v * v);
                nSW[0][x] = one36th * (1 - 3 * v + 3 * v * v);

                n0[ydim - 1][x] = four9ths * (1 - 1.5 * v * v);
                nE[ydim - 1][x] = one9th * (1 + 3 * v + 3 * v * v);
                nW[ydim - 1][x] = one9th * (1 - 3 * v + 3 * v * v);
                nN[ydim - 1][x] = one9th * (1 - 1.5 * v * v);
                nS[ydim - 1][x] = one9th * (1 - 1.5 * v * v);
                nNE[ydim - 1][x] = one36th * (1 + 3 * v + 3 * v * v);
                nSE[ydim - 1][x] = one36th * (1 + 3 * v + 3 * v * v);
                nNW[ydim - 1][x] = one36th * (1 - 3 * v + 3 * v * v);
                nSW[ydim - 1][x] = one36th * (1 - 3 * v + 3 * v * v);
            }
        }

        __global__
        void bounce(uint32_t xdim,
                    uint32_t ydim,
                    bool **barrier,
                    double **n0,
                    double **nN,
                    double **nS,
                    double **nE,
                    double **nW,
                    double **nNW,
                    double **nNE,
                    double **nSW,
                    double **nSE,
                    double **density,
                    double **xvel,
                    double **yvel,
                    double **speed2,
                    double **n0_temp,
                    double **nN_temp,
                    double **nS_temp,
                    double **nE_temp,
                    double **nW_temp,
                    double **nNW_temp,
                    double **nNE_temp,
                    double **nSW_temp,
                    double **nSE_temp,
                    double **density_temp,
                    double **xvel_temp,
                    double **yvel_temp,
                    double **speed2_temp,
                    double omega,
                    double v) {

            uint32_t y = blockIdx.x;
            uint32_t x = threadIdx.x;
            if(y < ydim && x < xdim && barrier[y][x]) {
                if (nN[y][x] > 0) {
                    nS[y - 1][x] += nN_temp[y][x];
                    nN[y][x] = 0;
                }
                if (nS[y][x] > 0) {
                    nN[y + 1][x] += nS_temp[y][x];
                    nS[y][x] = 0;
                }
                if (nE[y][x] > 0) {
                    nW[y][x - 1] += nE_temp[y][x];
                    nE[y][x] = 0;
                }
                if (nW[y][x] > 0) {
                    nE[y][x + 1] += nW_temp[y][x];
                    nW[y][x] = 0;
                }
                if (nNW[y][x] > 0) {
                    nSE[y - 1][x + 1] += nNW_temp[y][x];
                    nNW[y][x] = 0;
                }
                if (nNE[y][x] > 0) {
                    nSW[y - 1][x - 1] += nNE_temp[y][x];
                    nNE[y][x] = 0;
                }
                if (nSW[y][x] > 0) {
                    nNE[y + 1][x + 1] += nSW_temp[y][x];
                    nSW[y][x] = 0;
                }
                if (nSE[y][x] > 0) {
                    nNW[y + 1][x - 1] += nSE_temp[y][x];
                    nSE[y][x] = 0;
                }
            }
        }

        __global__
        void synchronize(uint32_t xdim,
                         uint32_t ydim,
                         bool **barrier,
                         double **n0,
                         double **nN,
                         double **nS,
                         double **nE,
                         double **nW,
                         double **nNW,
                         double **nNE,
                         double **nSW,
                         double **nSE,
                         double **density,
                         double **xvel,
                         double **yvel,
                         double **speed2,
                         double **n0_temp,
                         double **nN_temp,
                         double **nS_temp,
                         double **nE_temp,
                         double **nW_temp,
                         double **nNW_temp,
                         double **nNE_temp,
                         double **nSW_temp,
                         double **nSE_temp,
                         double **density_temp,
                         double **xvel_temp,
                         double **yvel_temp,
                         double **speed2_temp,
                         double omega,
                         double v) {

            uint32_t y = blockIdx.x;
            uint32_t x = threadIdx.x;

            density_temp[y][x] = density[y][x];
            xvel_temp[y][x] = xvel[y][x];
            yvel_temp[y][x] = yvel[y][x];
            speed2_temp[y][x] = speed2[y][x];
            n0_temp[y][x] = n0[y][x];
            nE_temp[y][x] = nE[y][x];
            nW_temp[y][x] = nW[y][x];
            nN_temp[y][x] = nN[y][x];
            nS_temp[y][x] = nS[y][x];
            nNE_temp[y][x] = nNE[y][x];
            nNW_temp[y][x] = nNW[y][x];
            nSE_temp[y][x] = nSE[y][x];
            nSW_temp[y][x] = nSW[y][x];
        }


        Simulation::Simulation(int width_, int height_) : xdim(width_), ydim(height_) {
            this->xdim = width_;
            this->ydim = height_;

            cudaMallocManaged(&this->n0, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nN, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nS, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nE, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nW, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nNW, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nNE, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nSW, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nSE, sizeof(double *) * ydim);
            cudaMallocManaged(&this->density, sizeof(double *) * ydim);
            cudaMallocManaged(&this->xvel, sizeof(double *) * ydim);
            cudaMallocManaged(&this->yvel, sizeof(double *) * ydim);
            cudaMallocManaged(&this->speed2, sizeof(double *) * ydim);
            cudaMallocManaged(&this->curl, sizeof(double *) * ydim);
            cudaMallocManaged(&this->barrier, sizeof(bool *) * ydim);

            cudaMallocManaged(&this->n0_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nN_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nS_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nE_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nW_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nNW_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nNE_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nSW_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->nSE_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->density_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->xvel_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->yvel_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->speed2_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->curl_temp, sizeof(double *) * ydim);
            cudaMallocManaged(&this->barrier_temp, sizeof(bool *) * ydim);

            for (int i = 0; i < this->ydim; i++) {
                cudaMallocManaged(&this->n0[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nN[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nS[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nE[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nW[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nNW[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nNE[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nSW[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nSE[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->density[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->xvel[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->yvel[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->speed2[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->curl[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->barrier[i], sizeof(barrier) * xdim);

                cudaMallocManaged(&this->n0_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nN_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nS_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nE_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nW_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nNW_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nNE_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nSW_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->nSE_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->density_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->xvel_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->yvel_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->speed2_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->curl_temp[i], sizeof(double) * xdim);
                cudaMallocManaged(&this->barrier_temp[i], sizeof(barrier) * xdim);

                cudaMemset(this->n0[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nN[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nS[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nW[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nE[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nNW[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nNE[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nSW[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nSE[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->density[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->xvel[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->yvel[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->speed2[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->curl[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->barrier[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(bool));

                cudaMemset(this->n0_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nN_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nS_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nW_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nE_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nNW_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nNE_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nSW_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->nSE_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->density_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->xvel_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->yvel_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->speed2_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->curl_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                cudaMemset(this->barrier_temp[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(bool));
            }

            cudaDeviceSynchronize();

            this->draw_circle(50, 150, 20);
            this->draw_circle(150, 200, 20);
            this->draw_circle(50, 250, 20);
            this->draw_circle(150, 300, 20);
            this->draw_circle(50, 350, 20);
            cudaDeviceSynchronize();

            this->init_fluid();
            cudaDeviceSynchronize();
        }

        Simulation::~Simulation() {
            for (int i = 0; i < this->ydim; i++) {
                cudaFree(this->n0[i]);
                cudaFree(this->nN[i]);
                cudaFree(this->nS[i]);
                cudaFree(this->nE[i]);
                cudaFree(this->nW[i]);
                cudaFree(this->nNW[i]);
                cudaFree(this->nNE[i]);
                cudaFree(this->nSW[i]);
                cudaFree(this->nSE[i]);
                cudaFree(this->density[i]);
                cudaFree(this->xvel[i]);
                cudaFree(this->yvel[i]);
                cudaFree(this->speed2[i]);
                cudaFree(this->curl[i]);
                cudaFree(this->barrier[i]);

                cudaFree(this->n0_temp[i]);
                cudaFree(this->nN_temp[i]);
                cudaFree(this->nS_temp[i]);
                cudaFree(this->nE_temp[i]);
                cudaFree(this->nW_temp[i]);
                cudaFree(this->nNW_temp[i]);
                cudaFree(this->nNE_temp[i]);
                cudaFree(this->nSW_temp[i]);
                cudaFree(this->nSE_temp[i]);
                cudaFree(this->density_temp[i]);
                cudaFree(this->xvel_temp[i]);
                cudaFree(this->yvel_temp[i]);
                cudaFree(this->speed2_temp[i]);
                cudaFree(this->curl_temp[i]);
                cudaFree(this->barrier_temp[i]);
            }
            cudaFree(this->n0);
            cudaFree(this->nN);
            cudaFree(this->nS);
            cudaFree(this->nE);
            cudaFree(this->nW);
            cudaFree(this->nNW);
            cudaFree(this->nNE);
            cudaFree(this->nSW);
            cudaFree(this->nSE);
            cudaFree(this->density);
            cudaFree(this->xvel);
            cudaFree(this->yvel);
            cudaFree(this->speed2);
            cudaFree(this->curl);
            cudaFree(this->barrier);

            cudaFree(this->n0_temp);
            cudaFree(this->nN_temp);
            cudaFree(this->nS_temp);
            cudaFree(this->nE_temp);
            cudaFree(this->nW_temp);
            cudaFree(this->nNW_temp);
            cudaFree(this->nNE_temp);
            cudaFree(this->nSW_temp);
            cudaFree(this->nSE_temp);
            cudaFree(this->density_temp);
            cudaFree(this->xvel_temp);
            cudaFree(this->yvel_temp);
            cudaFree(this->speed2_temp);
            cudaFree(this->curl_temp);
            cudaFree(this->barrier_temp);
        }

        void Simulation::init_fluid() const {
            for (int y = 0; y < ydim; y++) {
                for (int x = 0; x < xdim; x++) {
                    if (barrier[y][x]) {
                        zeroSite(x, y);
                    } else {
                        n0[y][x] = four9ths * (1 - 1.5 * v * v);
                        nE[y][x] = one9th * (1 + 3 * v + 3 * v * v);
                        nW[y][x] = one9th * (1 - 3 * v + 3 * v * v);
                        nN[y][x] = one9th * (1 - 1.5 * v * v);
                        nS[y][x] = one9th * (1 - 1.5 * v * v);
                        nNE[y][x] = one36th * (1 + 3 * v + 3 * v * v);
                        nSE[y][x] = one36th * (1 + 3 * v + 3 * v * v);
                        nNW[y][x] = one36th * (1 - 3 * v + 3 * v * v);
                        nSW[y][x] = one36th * (1 - 3 * v + 3 * v * v);
                        density[y][x] = 1;
                        xvel[y][x] = v;
                        yvel[y][x] = 0;
                        speed2[y][x] = v * v;
                    }
                }
            }
        }

        void Simulation::zeroSite(int x, int y) const {
            n0[y][x] = 0;
            nE[y][x] = 0;
            nW[y][x] = 0;
            nN[y][x] = 0;
            nS[y][x] = 0;
            nNE[y][x] = 0;
            nNW[y][x] = 0;
            nSE[y][x] = 0;
            nSW[y][x] = 0;
            xvel[y][x] = 0;
            yvel[y][x] = 0;
            speed2[y][x] = 0;


            n0_temp[y][x] = 0;
            nE_temp[y][x] = 0;
            nW_temp[y][x] = 0;
            nN_temp[y][x] = 0;
            nS_temp[y][x] = 0;
            nNE_temp[y][x] = 0;
            nNW_temp[y][x] = 0;
            nSE_temp[y][x] = 0;
            nSW_temp[y][x] = 0;
            xvel_temp[y][x] = 0;
            yvel_temp[y][x] = 0;
            speed2_temp[y][x] = 0;
        }

        void Simulation::collide() const {
            boltzmann::core::collide<<<this->ydim, this->xdim>>>(
                    xdim,
                            ydim,
                            barrier,
                            n0,
                            nN,
                            nS,
                            nE,
                            nW,
                            nNW,
                            nNE,
                            nSW,
                            nSE,
                            density,
                            xvel,
                            yvel,
                            speed2,
                            n0_temp,
                            nN_temp,
                            nS_temp,
                            nE_temp,
                            nW_temp,
                            nNW_temp,
                            nNE_temp,
                            nSW_temp,
                            nSE_temp,
                            density_temp,
                            xvel_temp,
                            yvel_temp,
                            speed2_temp,
                            omega);
        }

        void Simulation::stream() const {
            boltzmann::core::stream<<<this->ydim, this->xdim>>>(
                    xdim,
                            ydim,
                            barrier,
                            n0,
                            nN,
                            nS,
                            nE,
                            nW,
                            nNW,
                            nNE,
                            nSW,
                            nSE,
                            density,
                            xvel,
                            yvel,
                            speed2,
                            n0_temp,
                            nN_temp,
                            nS_temp,
                            nE_temp,
                            nW_temp,
                            nNW_temp,
                            nNE_temp,
                            nSW_temp,
                            nSE_temp,
                            density_temp,
                            xvel_temp,
                            yvel_temp,
                            speed2_temp,
                            omega,
                            v);
        }

        void Simulation::bounce() const {
            boltzmann::core::bounce<<<this->ydim, this->xdim>>>(
                    xdim,
                            ydim,
                            barrier,
                            n0,
                            nN,
                            nS,
                            nE,
                            nW,
                            nNW,
                            nNE,
                            nSW,
                            nSE,
                            density,
                            xvel,
                            yvel,
                            speed2,
                            n0_temp,
                            nN_temp,
                            nS_temp,
                            nE_temp,
                            nW_temp,
                            nNW_temp,
                            nNE_temp,
                            nSW_temp,
                            nSE_temp,
                            density_temp,
                            xvel_temp,
                            yvel_temp,
                            speed2_temp,
                            omega,
                            v);
        }

        void Simulation::compute_curl() const {
            boltzmann::core::compute_curl<<<this->ydim, this->xdim>>>(xdim, ydim, curl, yvel, xvel);
        }

        void Simulation::draw_barrier(int x, int y) const {
            barrier[y][x] = true;
            zeroSite(x, y);
        }

        void Simulation::draw_circle(int x_center, int y_center, int radius) const {
            for (int y = 0; y < ydim; y++) {
                for (int x = 0; x < xdim; x++) {
                    if (std::sqrt((x - x_center) * (x - x_center) + (y - y_center) * (y - y_center)) < radius) {
                        draw_barrier(x, y);
                    }
                }
            }
        }

        void Simulation::synchronize() const {
            boltzmann::core::synchronize<<<this->ydim, this->xdim>>>(
                    xdim,
                            ydim,
                            barrier,
                            n0,
                            nN,
                            nS,
                            nE,
                            nW,
                            nNW,
                            nNE,
                            nSW,
                            nSE,
                            density,
                            xvel,
                            yvel,
                            speed2,
                            n0_temp,
                            nN_temp,
                            nS_temp,
                            nE_temp,
                            nW_temp,
                            nNW_temp,
                            nNE_temp,
                            nSW_temp,
                            nSE_temp,
                            density_temp,
                            xvel_temp,
                            yvel_temp,
                            speed2_temp,
                            omega,
                            v);
        }
    }
}
