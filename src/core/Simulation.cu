#include "core/Simulation.h"
#include "core/Kernels.h"

namespace boltzmann {
    namespace core {
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


            this->draw_square(300, 0,310, 300);
            this->draw_square(350, 200,360, 499);
            this->synchronize();
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

        void Simulation::draw_square(int x1, int y1, int x2, int y2) const {

        }
    }
}
