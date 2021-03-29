#include "core/Simulation.h"

namespace boltzmann {
    namespace core {
        Simulation::Simulation(int width_, int height_) : xdim(width_), ydim(height_) {
            this->xdim = width_;
            this->ydim = height_;
            this->n0 = new double *[this->ydim];
            this->nN = new double *[this->ydim];
            this->nS = new double *[this->ydim];
            this->nE = new double *[this->ydim];
            this->nW = new double *[this->ydim];
            this->nNW = new double *[this->ydim];
            this->nNE = new double *[this->ydim];
            this->nSW = new double *[this->ydim];
            this->nSE = new double *[this->ydim];
            this->density = new double *[this->ydim];
            this->xvel = new double *[this->ydim];
            this->yvel = new double *[this->ydim];
            this->speed2 = new double *[this->ydim];
            this->curl = new double *[this->ydim];
            this->barrier = new bool *[this->ydim];

            for (int i = 0; i < this->ydim; i++) {
                this->n0[i] = new double[this->xdim];
                this->nN[i] = new double[this->xdim];
                this->nS[i] = new double[this->xdim];
                this->nE[i] = new double[this->xdim];
                this->nW[i] = new double[this->xdim];
                this->nNW[i] = new double[this->xdim];
                this->nNE[i] = new double[this->xdim];
                this->nSW[i] = new double[this->xdim];
                this->nSE[i] = new double[this->xdim];
                this->density[i] = new double[this->xdim];
                this->xvel[i] = new double[this->xdim];
                this->yvel[i] = new double[this->xdim];
                this->speed2[i] = new double[this->xdim];
                this->curl[i] = new double[this->xdim];
                this->barrier[i] = new bool[this->xdim];

                std::memset(this->n0[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->nN[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->nS[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->nW[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->nE[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->nNW[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->nNE[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->nSW[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->nSE[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->density[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->xvel[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->yvel[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->speed2[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->curl[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(double));
                std::memset(this->barrier[i], 0, static_cast<unsigned long>(this->xdim) * sizeof(bool));
            }

            this->draw_circle(50, 150, 20);
            this->draw_circle(150, 200, 20);
            this->draw_circle(50, 250, 20);
            this->draw_circle(150, 300, 20);
            this->draw_circle(50, 350, 20);

            //this->draw_square(300, 0,310, 300);
            //this->draw_square(350, 200,360, 499);

            this->init_fluid();
        }

        Simulation::~Simulation() {
            for (int i = 0; i < this->ydim; i++) {
                delete this->n0[i];
                delete this->nN[i];
                delete this->nS[i];
                delete this->nE[i];
                delete this->nW[i];
                delete this->nNW[i];
                delete this->nNE[i];
                delete this->nSW[i];
                delete this->nSE[i];
                delete this->density[i];
                delete this->xvel[i];
                delete this->yvel[i];
                delete this->speed2[i];
                delete this->curl[i];
                delete this->barrier[i];
            }
            delete this->n0;
            delete this->nN;
            delete this->nS;
            delete this->nE;
            delete this->nW;
            delete this->nNW;
            delete this->nNE;
            delete this->nSW;
            delete this->nSE;
            delete this->density;
            delete this->xvel;
            delete this->yvel;
            delete this->speed2;
            delete this->curl;
            delete this->barrier;
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

        void Simulation::collide() {
#pragma omp parallel for
            for (int y = 0; y < ydim; y++) {
                for (int x = 0; x < xdim; x++) {
                    double n, one9thn, one36thn, vx, vy, vx2, vy2, vx3, vy3, vxvy2, v2, v215;
                    if (!barrier[y][x]) {
                        n = n0[y][x] + nN[y][x] + nS[y][x] + nE[y][x] + nW[y][x] + nNW[y][x] + nNE[y][x] + nSW[y][x] +
                            nSE[y][x];
                        density[y][x] = n;        // macroscopic density may be needed for plotting
                        one9thn = one9th * n;
                        one36thn = one36th * n;
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
            for (int y = 1; y < ydim - 1; y++) {
                for (int x = 1; x < xdim - 1; x++) {
                    curl[y][x] = (yvel[y][x + 1] - yvel[y][x - 1]) - (xvel[y + 1][x] - xvel[y - 1][x]);
                }
            }
            for (int y = 1; y < ydim - 1; y++) {
                curl[y][0] = 2 * (yvel[y][1] - yvel[y][0]) - (xvel[y + 1][0] - xvel[y - 1][0]);
                curl[y][xdim - 1] =
                        2 * (yvel[y][xdim - 1] - yvel[y][xdim - 2]) - (xvel[y + 1][xdim - 1] - xvel[y - 1][xdim - 1]);
            }
        }

        void Simulation::draw_square(int x1, int y1, int x2, int y2) const {
            for (int y = 0; y < ydim; y++) {
                for (int x = 0; x < xdim; x++) {
                    if (x > x1 && y > y1 && x < x2 && y < y2) {
                        draw_barrier(x, y);
                    }
                }
            }
        }
    }
}
