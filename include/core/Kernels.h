#ifndef LATTICE_BOLTZMANN_KERNELS_H
#define LATTICE_BOLTZMANN_KERNELS_H

namespace boltzmann {
    namespace core {
        const double four9ths = 4.0 / 9;
        const double one9th = 1.0 / 9;
        const double one36th = 1.0 / 36;

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

            if (y < ydim && x < xdim) {
                double n, one9thn, one36thn, vx, vy, vx2, vy2, vx3, vy3, vxvy2, v2, v215;
                if (!barrier[y][x]) {
                    n = n0[y][x] + nN[y][x] + nS[y][x] + nE[y][x] + nW[y][x] + nNW[y][x] + nNE[y][x] + nSW[y][x] +
                        nSE[y][x];
                    one9thn = one9th * n;
                    one36thn = one36th * n;
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
            if (y >= 1 && y < ydim - 1 && x == 0) {
                curl[y][0] = 2 * (yvel[y][1] - yvel[y][0]) - (xvel[y + 1][0] - xvel[y - 1][0]);
            }
            if (y >= 1 && y < ydim - 1 && x == xdim - 1) {
                curl[y][xdim - 1] =
                        2 * (yvel[y][xdim - 1] - yvel[y][xdim - 2]) - (xvel[y + 1][xdim - 1] - xvel[y - 1][xdim - 1]);
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

            if (y >= 0 && y < ydim && !barrier[y][0] && x == xdim - 1) {
                nW[y][xdim - 1] = one9th * (1 - 3 * v + 3 * v * v);
                nNW[y][xdim - 1] = one36th * (1 - 3 * v + 3 * v * v);
                nSW[y][xdim - 1] = one36th * (1 - 3 * v + 3 * v * v);
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
            if (y < ydim && x < xdim && barrier[y][x]) {
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

    }
}

#endif //LATTICE_BOLTZMANN_KERNELS_H
