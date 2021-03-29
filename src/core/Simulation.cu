#include "core/Simulation.h"


namespace boltzmann {
    namespace core {
        __global__
        void collide(uint32_t xdim, uint32_t ydim, bool **barrier, double **n0, double **nN, double **nS, double **nE,
                     double **nW,
                     double **nNW, double **nNE, double **nSW, double **nSE, double **density, double **xvel,
                     double **yvel,
                     double **speed2, double omega) {
            uint32_t y = blockIdx.x;
            uint32_t x = threadIdx.x;

            if(y < ydim && x < xdim) {
                double n, one9thn, one36thn, vx, vy, vx2, vy2, vx3, vy3, vxvy2, v2, v215;
                if (!barrier[y][x]) {
                    n = n0[y][x] + nN[y][x] + nS[y][x] + nE[y][x] + nW[y][x] + nNW[y][x] + nNE[y][x] + nSW[y][x] +
                        nSE[y][x];
                    density[y][x] = n;        // macroscopic density may be needed for plotting
                    one9thn = 1.0 / 9.0 * n;
                    one36thn = 1.0 / 36.0 * n;
                    if (n > 0) {
                        vx = (nE[y][x] + nNE[y][x] + nSE[y][x] - nW[y][x] - nNW[y][x] - nSW[y][x]) / n;
                        vy = (nN[y][x] + nNE[y][x] + nNW[y][x] - nS[y][x] - nSE[y][x] - nSW[y][x]) / n;
                    } else {
                        vx = 0;
                        vy = 0;
                    }
                    xvel[y][x] = vx;        // may be needed for plotting
                    yvel[y][x] = vy;        // may be needed for plotting
                    vx3 = 3 * vx;
                    vy3 = 3 * vy;
                    vx2 = vx * vx;
                    vy2 = vy * vy;
                    vxvy2 = 2 * vx * vy;
                    v2 = vx2 + vy2;
                    speed2[y][x] = v2;        // may be needed for plotting
                    v215 = 1.5 * v2;
                    n0[y][x] += omega * (4.0 / 9.0 * n * (1 - v215) - n0[y][x]);
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
        void compute_curl(uint32_t xdim, uint32_t ydim, double**curl, double**yvel, double**xvel) {
            uint32_t y = blockIdx.x;
            uint32_t x = threadIdx.x;
            if(y >= 1 && y < ydim - 1 && y < ydim && x < xdim) {
                curl[y][x] = (yvel[y][x + 1] - yvel[y][x - 1]) - (xvel[y + 1][x] - xvel[y - 1][x]);
            }

            if(y >= 1 && y < ydim - 1 && x == 0) {
                curl[y][0] = 2 * (yvel[y][1] - yvel[y][0]) - (xvel[y + 1][0] - xvel[y - 1][0]);
            }

            if(y >= 1 && y < ydim - 1 && x == xdim - 1) {
                curl[y][xdim - 1] = 2 * (yvel[y][xdim - 1] - yvel[y][xdim - 2]) - (xvel[y + 1][xdim - 1] - xvel[y - 1][xdim - 1]);
            }
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
        }

        void Simulation::collide() const {
            boltzmann::core::collide<<<this->ydim, this->xdim>>>(xdim, ydim, barrier, n0, nN, nS, nE, nW,
                    nNW, nNE, nSW, nSE, density, xvel, yvel,
                    speed2, omega);
        }

        void Simulation::stream() const {
            for (int y = ydim - 1; y > 0; y--) {
                for (int x = 0; x < xdim - 1; x++) {        // first start in NW corner...
                    nN[y][x] = nN[y - 1][x];        // move the north-moving particles
                    nNW[y][x] = nNW[y - 1][x + 1];    // and the northwest-moving particles
                }
            }
            for (int y = ydim - 1; y > 0; y--) {
                for (int x = xdim - 1; x > 0; x--) {        // now start in NE corner...
                    nE[y][x] = nE[y][x - 1];        // move the east-moving particles
                    nNE[y][x] = nNE[y - 1][x - 1];    // and the northeast-moving particles
                }
            }
            for (int y = 0; y < ydim - 1; y++) {
                for (int x = xdim - 1; x > 0; x--) {        // now start in SE corner...
                    nS[y][x] = nS[y + 1][x];        // move the south-moving particles
                    nSE[y][x] = nSE[y + 1][x - 1];    // and the southeast-moving particles
                }
            }
            for (int y = 0; y < ydim - 1; y++) {
                for (int x = 0; x < xdim - 1; x++) {        // now start in the SW corner...
                    nW[y][x] = nW[y][x + 1];        // move the west-moving particles
                    nSW[y][x] = nSW[y + 1][x + 1];    // and the southwest-moving particles
                }
            }

            // We missed a few at the left and right edges:
            for (int y = 0; y < ydim - 1; y++) {
                nS[y][0] = nS[y + 1][0];
            }
            for (int y = ydim - 1; y > 0; y--) {
                nN[y][xdim - 1] = nN[y - 1][xdim - 1];
            }

            // Stream particles in from the non-existent space to the left
            for (int y = 0; y < ydim; y++) {
                if (!barrier[y][0]) {
                    nE[y][0] = one9th * (1 + 3 * v + 3 * v * v);
                    nNE[y][0] = one36th * (1 + 3 * v + 3 * v * v);
                    nSE[y][0] = one36th * (1 + 3 * v + 3 * v * v);
                }
            }

            for (int y = 0; y < ydim; y++) {
                if (!barrier[y][0]) {
                    nW[y][xdim - 1] = one9th * (1 - 3 * v + 3 * v * v);
                    nNW[y][xdim - 1] = one36th * (1 - 3 * v + 3 * v * v);
                    nSW[y][xdim - 1] = one36th * (1 - 3 * v + 3 * v * v);
                }
            }

            // Now handle top and bottom edges:
            for (int x = 0; x < xdim; x++) {
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

        void Simulation::bounce() const {
            for (int y = 0; y < ydim; y++) {
                for (int x = 0; x < xdim; x++) {
                    if (barrier[y][x]) {
                        if (nN[y][x] > 0) {
                            nS[y - 1][x] += nN[y][x];
                            nN[y][x] = 0;
                        }
                        if (nS[y][x] > 0) {
                            nN[y + 1][x] += nS[y][x];
                            nS[y][x] = 0;
                        }
                        if (nE[y][x] > 0) {
                            nW[y][x - 1] += nE[y][x];
                            nE[y][x] = 0;
                        }
                        if (nW[y][x] > 0) {
                            nE[y][x + 1] += nW[y][x];
                            nW[y][x] = 0;
                        }
                        if (nNW[y][x] > 0) {
                            nSE[y - 1][x + 1] += nNW[y][x];
                            nNW[y][x] = 0;
                        }
                        if (nNE[y][x] > 0) {
                            nSW[y - 1][x - 1] += nNE[y][x];
                            nNE[y][x] = 0;
                        }
                        if (nSW[y][x] > 0) {
                            nNE[y + 1][x + 1] += nSW[y][x];
                            nSW[y][x] = 0;
                        }
                        if (nSE[y][x] > 0) {
                            nNW[y + 1][x - 1] += nSE[y][x];
                            nSE[y][x] = 0;
                        }
                    }
                }
            }
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

        void Simulation::compute_curl() const {
            boltzmann::core::compute_curl<<<this->ydim, this->xdim>>>(xdim, ydim, curl, yvel, xvel);
        }
    }
}
