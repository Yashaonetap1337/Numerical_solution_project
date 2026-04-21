#include "MaderSolver.h"
#include "euler_utils.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

using Field2D = std::vector<std::vector<double>>;

static double s_gamm;
static double s_Rsp;
static double s_P_min;
static double s_Q_chem;
static double s_VISC;
static double s_GASW;
static double s_MINGRHO;
static double s_P_thresh;
static int    s_Nx, s_Ny, s_fict;

static Field2D s_mass_fraction;
static std::vector<std::vector<bool>> s_reacted;
static bool s_state_initialized = false;

static constexpr double k_eps    = 1e-14;
static constexpr double WALL_RHO = 100.0;

static inline bool IsWall(double rho) { return rho > WALL_RHO; }


static void FillGhost(Field2D& F, int NX, int NY) {
    const int il = s_fict, ir = s_Nx + s_fict - 1;
    const int jb = s_fict, jt = s_Ny + s_fict - 1;
    for (int j = 0; j < NY; ++j)
        for (int g = 0; g < s_fict; ++g) {
            F[g][j]          = F[il][j];
            F[NX - 1 - g][j] = F[ir][j];
        }
    for (int i = 0; i < NX; ++i)
        for (int g = 0; g < s_fict; ++g) {
            F[i][g]          = F[i][jb];
            F[i][NY - 1 - g] = F[i][jt];
        }
}

// CHECK: MADER_EOS
static void EOS(const std::vector<std::vector<State>>& W,
                Field2D& rho, Field2D& U, Field2D& V,
                Field2D& P, Field2D& T, Field2D& I,
                int NX, int NY)
{
    const double gm1 = s_gamm - 1.0;
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const double r = std::max(W[i][j].rho, 1e-14);
            const double p = W[i][j].p;
            rho[i][j] = r;
            U[i][j]   = W[i][j].u;
            V[i][j]   = W[i][j].v;
            I[i][j]   = p / (gm1 * r);
            P[i][j]   = std::max(gm1 * r * I[i][j], s_P_min);
            T[i][j]   = P[i][j] / (r * s_Rsp);
        }
    }
}

// CHECK: MADER_ARRHENIUS
static void InstantReaction(const Field2D& P, const Field2D& rho, Field2D& I)
{
    for (int i = s_fict; i < s_Nx + s_fict; ++i) {
        for (int j = s_fict; j < s_Ny + s_fict; ++j) {
            if (IsWall(rho[i][j])) continue;
            if (s_reacted[i][j])   continue;

            const double mf = s_mass_fraction[i][j];
            if (mf <= s_GASW) { s_reacted[i][j] = true; continue; }

            if (P[i][j] > s_P_thresh) {
                I[i][j]              += s_Q_chem * mf;
                s_mass_fraction[i][j] = 0.0;
                s_reacted[i][j]       = true;
            }
        }
    }
}

// CHECK: MADER_VISC
static void ComputeViscosity(const Field2D& U, const Field2D& V,
                              const Field2D& rho,
                              Field2D& Q1, Field2D& Q2,
                              Field2D& Q3, Field2D& Q4)
{
    for (int i = s_fict; i < s_Nx + s_fict; ++i)
        for (int j = s_fict; j < s_Ny + s_fict; ++j) {
            if (IsWall(rho[i][j])) { Q1[i][j] = Q2[i][j] = 0.0; continue; }

            Q1[i][j] = (!IsWall(rho[i][j + 1]) && V[i][j] >= V[i][j + 1])
                ? s_VISC * rho[i][j] * (V[i][j] - V[i][j + 1]) : 0.0;

            Q2[i][j] = (!IsWall(rho[i + 1][j]) && U[i][j] >= U[i + 1][j])
                ? s_VISC * rho[i][j] * (U[i][j] - U[i + 1][j]) : 0.0;
        }
    for (int i = s_fict; i < s_Nx + s_fict; ++i)
        for (int j = s_fict; j < s_Ny + s_fict; ++j) {
            Q3[i][j] = (j > s_fict) ? Q1[i][j - 1] : 0.0;
            Q4[i][j] = (i > s_fict) ? Q2[i - 1][j] : 0.0;
        }
}

// CHECK: MADER_VELOCITY
static void VelocityTilde(const Field2D& U, const Field2D& V,
                           const Field2D& P, const Field2D& rho,
                           const Field2D& Q1, const Field2D& Q2,
                           const Field2D& Q3, const Field2D& Q4,
                           Field2D& U_tilde, Field2D& V_tilde,
                           const std::vector<double>& x,
                           const std::vector<double>& y,
                           double dt)
{
    for (int i = s_fict; i < s_Nx + s_fict; ++i)
        for (int j = s_fict; j < s_Ny + s_fict; ++j) {
            if (IsWall(rho[i][j])) {
                U_tilde[i][j] = 0.0; V_tilde[i][j] = 0.0; continue;
            }
            const double dR = x[i] - x[i - 1];
            const double dZ = y[j] - y[j - 1];
            const double r  = std::max(rho[i][j], 1e-14);

            V_tilde[i][j] = V[i][j]
                - dt / (r * 2.0 * dZ) * ((P[i][j + 1] - P[i][j - 1]) +
                                         (Q3[i][j]    - Q1[i][j]));
            U_tilde[i][j] = U[i][j]
                - dt / (r * 2.0 * dR) * ((P[i + 1][j] - P[i - 1][j]) +
                                         (Q4[i][j]    - Q2[i][j]));
        }
}

// CHECK: MADER_ZIP
static void ZIPEnergy(const Field2D& I, Field2D& I_tilde,
                       const Field2D& P, const Field2D& rho,
                       const Field2D& U, const Field2D& V,
                       const Field2D& U_tilde, const Field2D& V_tilde,
                       const Field2D& Q1, const Field2D& Q2,
                       const Field2D& Q3, const Field2D& Q4,
                       const std::vector<double>& x,
                       const std::vector<double>& y,
                       double dt)
{
    for (int i = s_fict; i < s_Nx + s_fict; ++i)
        for (int j = s_fict; j < s_Ny + s_fict; ++j) {
            if (IsWall(rho[i][j])) { I_tilde[i][j] = I[i][j]; continue; }
            const double dR = x[i] - x[i - 1];
            const double dZ = y[j] - y[j - 1];
            const double r  = std::max(rho[i][j], 1e-14);

            const double U1 = U[i - 1][j] + U_tilde[i - 1][j];
            const double U2 = U[i + 1][j] + U_tilde[i + 1][j];
            const double V1 = V[i][j - 1] + V_tilde[i][j - 1];
            const double V2 = V[i][j + 1] + V_tilde[i][j + 1];
            const double T3 = U[i][j]     + U_tilde[i][j];
            const double T1 = V[i][j]     + V_tilde[i][j];

            const double zip = (P[i][j]  / dR) * (U2 - U1)
                             + (Q4[i][j] / dR) * (U2 - T3)
                             + (Q2[i][j] / dR) * (T3 - U1)
                             + (P[i][j]  / dZ) * (V2 - V1)
                             + (Q3[i][j] / dZ) * (V2 - T1)
                             + (Q1[i][j] / dZ) * (T1 - V1);

            I_tilde[i][j] = std::max(I[i][j] - dt / (4.0 * r) * zip, 0.0);
        }
}

static void TransportAndRepartition(
    Field2D& rho, Field2D& U, Field2D& V, Field2D& I,
    const Field2D& U_tilde, const Field2D& V_tilde,
    const Field2D& I_tilde,
    const std::vector<double>& x, const std::vector<double>& y,
    double dt, int NX, int NY)
{
    Field2D DM (NX, std::vector<double>(NY, 0.0));
    Field2D DE (NX, std::vector<double>(NY, 0.0));
    Field2D DW (NX, std::vector<double>(NY, 0.0));
    Field2D DPU(NX, std::vector<double>(NY, 0.0));
    Field2D DPV(NX, std::vector<double>(NY, 0.0));

    // CHECK: MADER_DONOR
    for (int i = s_fict; i < s_Nx + s_fict; ++i) {
        for (int j = s_fict; j < s_Ny + s_fict; ++j) {
            if (IsWall(rho[i][j])) continue;
            const double dR = x[i] - x[i - 1];
            const double dZ = y[j] - y[j - 1];

            bool skip_r = (i > s_fict && IsWall(rho[i - 1][j]));
            if (!skip_r) {
                double denom = 1.0 + (U_tilde[i - 1][j] - U_tilde[i][j]) * dt / dR;
                if (std::abs(denom) < 1e-14) denom = (denom >= 0 ? 1e-14 : -1e-14);
                const double alpha = 0.5 * (U_tilde[i - 1][j] + U_tilde[i][j])
                                     * dt / dR / denom;

                const int di = (alpha >= 0.0) ? i - 1 : i;
                const int ai = (alpha >= 0.0) ? i     : i - 1;

                if (!IsWall(rho[di][j])) {
                    const double rho_d = std::max(rho[di][j], 1e-14);
                    const double DMASS = rho_d * std::abs(alpha);
                    const double E_d   = I_tilde[di][j]
                                       + 0.5 * (U_tilde[di][j] * U_tilde[di][j]
                                              + V_tilde[di][j] * V_tilde[di][j]);
                    const double mf_d  = s_mass_fraction[di][j];

                    DM [ai][j] += DMASS;      DM [di][j] -= DMASS;
                    DE [ai][j] += E_d * DMASS; DE [di][j] -= E_d * DMASS;
                    DPU[ai][j] += U_tilde[di][j] * DMASS;
                    DPU[di][j] -= U_tilde[di][j] * DMASS;
                    DPV[ai][j] += V_tilde[di][j] * DMASS;
                    DPV[di][j] -= V_tilde[di][j] * DMASS;
                    DW [ai][j] += mf_d * DMASS;
                    DW [di][j] -= mf_d * DMASS;
                }
            }

            bool skip_z = (j > s_fict && IsWall(rho[i][j - 1]));
            if (!skip_z) {
                double denom = 1.0 + (V_tilde[i][j - 1] - V_tilde[i][j]) * dt / dZ;
                if (std::abs(denom) < 1e-14) denom = (denom >= 0 ? 1e-14 : -1e-14);
                const double beta = 0.5 * (V_tilde[i][j - 1] + V_tilde[i][j])
                                    * dt / dZ / denom;

                const int dj = (beta >= 0.0) ? j - 1 : j;
                const int aj = (beta >= 0.0) ? j     : j - 1;

                if (!IsWall(rho[i][dj])) {
                    const double rho_d = std::max(rho[i][dj], 1e-14);
                    const double DMASS = rho_d * std::abs(beta);
                    const double E_d   = I_tilde[i][dj]
                                       + 0.5 * (U_tilde[i][dj] * U_tilde[i][dj]
                                              + V_tilde[i][dj] * V_tilde[i][dj]);
                    const double mf_d  = s_mass_fraction[i][dj];

                    DM [i][aj] += DMASS;      DM [i][dj] -= DMASS;
                    DE [i][aj] += E_d * DMASS; DE [i][dj] -= E_d * DMASS;
                    DPU[i][aj] += U_tilde[i][dj] * DMASS;
                    DPU[i][dj] -= U_tilde[i][dj] * DMASS;
                    DPV[i][aj] += V_tilde[i][dj] * DMASS;
                    DPV[i][dj] -= V_tilde[i][dj] * DMASS;
                    DW [i][aj] += mf_d * DMASS;
                    DW [i][dj] -= mf_d * DMASS;
                }
            }
        }
    }

    // CHECK: MADER_REPARTITION
    for (int i = s_fict; i < s_Nx + s_fict; ++i) {
        for (int j = s_fict; j < s_Ny + s_fict; ++j) {
            if (IsWall(rho[i][j])) continue;

            const double rho_new = rho[i][j] + DM[i][j];
            if (rho_new <= s_MINGRHO) {
                rho[i][j] = s_MINGRHO;
                U[i][j]   = 0.0; V[i][j] = 0.0; I[i][j] = 0.0;
                s_mass_fraction[i][j] = 0.0;
                continue;
            }

            const double U_new = (rho[i][j] * U_tilde[i][j] + DPU[i][j]) / rho_new;
            const double V_new = (rho[i][j] * V_tilde[i][j] + DPV[i][j]) / rho_new;
            const double E_old = I_tilde[i][j]
                               + 0.5 * (U_tilde[i][j] * U_tilde[i][j]
                                      + V_tilde[i][j] * V_tilde[i][j]);
            double I_new = (rho[i][j] * E_old + DE[i][j]) / rho_new
                         - 0.5 * (U_new * U_new + V_new * V_new);
            I_new = std::max(I_new, 0.0);

            const double mf_new = (rho[i][j] * s_mass_fraction[i][j] + DW[i][j])
                                  / rho_new;

            rho[i][j]             = rho_new;
            U[i][j]               = U_new;
            V[i][j]               = V_new;
            I[i][j]               = I_new;
            s_mass_fraction[i][j] = s_reacted[i][j]
                                    ? 0.0
                                    : std::max(mf_new, 0.0);
        }
    }
}

OmegaField make_omega_field(const Grid& grid) {
    return OmegaField(grid.Nx + 2 * grid.num_fict,
                      std::vector<double>(grid.Ny + 2 * grid.num_fict, 1.0));
}

void MaderTimeStep(const Grid& grid, const Config& cfg,
                   const std::vector<std::vector<State>>& W,
                         std::vector<std::vector<State>>& W_new,
                   OmegaField& omega, double dt, double /*t*/)
{
    s_gamm    = cfg.phys.gamma;
    s_VISC    = cfg.phys.VISC;
    s_GASW    = cfg.phys.GASW;
    s_MINGRHO = cfg.phys.MINGRHO;
    s_P_min   = 1e-10;
    s_Nx      = grid.Nx;
    s_Ny      = grid.Ny;
    s_fict    = grid.num_fict;

    s_Rsp = cfg.phys.Rgas * 1e-5 / cfg.phys.M_molar;

    s_Q_chem = cfg.phys.left_bottom.p
             / (cfg.phys.right_bottom.rho * 2.0 * (s_gamm - 1.0));

    s_P_thresh = 0.15 * cfg.phys.left_bottom.p;

    const int NX = s_Nx + 2 * s_fict;
    const int NY = s_Ny + 2 * s_fict;
    const std::vector<double>& x = grid.x_centers;
    const std::vector<double>& y = grid.y_centers;

    if (!s_state_initialized) {
        s_mass_fraction.assign(NX, std::vector<double>(NY, 1.0));
        s_reacted.assign(NX, std::vector<bool>(NY, false));
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < NY; ++j) {
                s_mass_fraction[i][j] = omega[i][j];
                s_reacted[i][j]       = (omega[i][j] <= s_GASW);
            }
        s_state_initialized = true;

        static bool init_printed = false;
        if (!init_printed) {
            init_printed = true;
            const double D_CJ = std::sqrt(cfg.phys.left_bottom.p *
                                          (s_gamm + 1.0) / cfg.phys.right_bottom.rho);
            std::cout << "[MADER INIT]"
                      << " Q_chem=" << s_Q_chem
                      << " D_CJ=" << D_CJ
                      << " P_CJ=" << cfg.phys.left_bottom.p
                      << " rho0=" << cfg.phys.right_bottom.rho
                      << " gamma=" << s_gamm
                      << " Rsp=" << s_Rsp
                      << " P_thresh=" << s_P_thresh << "\n";
        }
    }

    Field2D rho(NX, std::vector<double>(NY));
    Field2D U  (NX, std::vector<double>(NY));
    Field2D V  (NX, std::vector<double>(NY));
    Field2D P  (NX, std::vector<double>(NY));
    Field2D T  (NX, std::vector<double>(NY));
    Field2D I  (NX, std::vector<double>(NY));

    EOS(W, rho, U, V, P, T, I, NX, NY);
    FillGhost(P, NX, NY);
    FillGhost(T, NX, NY);
    FillGhost(I, NX, NY);

    InstantReaction(P, rho, I);
    FillGhost(I, NX, NY);

    Field2D Q1(NX, std::vector<double>(NY, 0.0));
    Field2D Q2(NX, std::vector<double>(NY, 0.0));
    Field2D Q3(NX, std::vector<double>(NY, 0.0));
    Field2D Q4(NX, std::vector<double>(NY, 0.0));
    ComputeViscosity(U, V, rho, Q1, Q2, Q3, Q4);
    FillGhost(Q1, NX, NY); FillGhost(Q2, NX, NY);
    FillGhost(Q3, NX, NY); FillGhost(Q4, NX, NY);

    Field2D U_tilde = U;
    Field2D V_tilde = V;
    VelocityTilde(U, V, P, rho, Q1, Q2, Q3, Q4,
                  U_tilde, V_tilde, x, y, dt);
    FillGhost(U_tilde, NX, NY);
    FillGhost(V_tilde, NX, NY);

    Field2D I_tilde = I;
    ZIPEnergy(I, I_tilde, P, rho, U, V, U_tilde, V_tilde,
              Q1, Q2, Q3, Q4, x, y, dt);
    FillGhost(I_tilde, NX, NY);

    TransportAndRepartition(rho, U, V, I,
                             U_tilde, V_tilde, I_tilde,
                             x, y, dt, NX, NY);

    W_new = W;
    const double gm1 = s_gamm - 1.0;
    for (int i = s_fict; i < s_Nx + s_fict; ++i) {
        for (int j = s_fict; j < s_Ny + s_fict; ++j) {
            if (IsWall(rho[i][j])) { W_new[i][j] = W[i][j]; continue; }

            const double r_n = std::max(rho[i][j], s_MINGRHO);
            const double u_n = U[i][j];
            const double v_n = V[i][j];
            const double I_n = std::max(I[i][j], 0.0);
            const double P_n = std::max(gm1 * r_n * I_n, s_P_min);

            W_new[i][j].rho = r_n;
            W_new[i][j].u   = u_n;
            W_new[i][j].v   = v_n;
            W_new[i][j].p   = P_n;
        }
    }

    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NY; ++j)
            omega[i][j] = s_mass_fraction[i][j];
}
