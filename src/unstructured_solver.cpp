// unstructured_solver.cpp
#include "unstructured_solver.h"
#include "euler_utils.h"
#include "choice_of_riemann_solvers.h"   // äëÿ compute_interface_flux
#include <iostream>
#include <fstream>
//CHECK: GHOST_STATE
// ------------------------------------------------------------
//  ghost_state – âû÷èñëåíèå ïðèçðà÷íîãî ñîñòîÿíèÿ äëÿ ãðàíèöû
// ------------------------------------------------------------
static State ghost_state(const State& WL, double nx, double ny,
    BoundaryType bc, const Config& cfg) {
    double un = WL.u * nx + WL.v * ny;
    switch (bc) {
    case BoundaryType::WALL: {
        State Wg = WL;
        Wg.u -= 2.0 * un * nx;
        Wg.v -= 2.0 * un * ny;
        return Wg;
    }
    case BoundaryType::INFLOW: {
        State Wg = cfg.phys.inflow_state;
        Wg.p = WL.p;                  
        return Wg;
    }
    case BoundaryType::OUTFLOW: {
        // Ôèêñèðóåì ñòàòè÷åñêîå äàâëåíèå, îñòàëüíîå ýêñòðàïîëèðóåì èç âíóòðåííåé ÿ÷åéêè
        State Wg = WL;
        Wg.p = cfg.phys.outflow_pressure;  // íóæíî äîáàâèòü â Config
        // Ïëîòíîñòü è ñêîðîñòü íå ìåíÿåì (èëè ìîæíî ýêñòðàïîëèðîâàòü)
        return Wg;
    }
    case BoundaryType::FREE:
    default:
        return WL;
    }
}
//CHECK: UNSTRUCT_SCHEMES
// ------------------------------------------------------------
//  compute_residuals_unstructured – îáùàÿ äëÿ âñåõ ìåòîäîâ,
//  çàâèñèò îò cfg.method, cfg.riemann_solver_type
// ------------------------------------------------------------
void compute_residuals_unstructured(Mesh& mesh, const Config& cfg) {
    mesh.zero_res();
    const double gamma = cfg.phys.gamma;

    for (const Face& f : mesh.faces) {
        // Ëåâàÿ ÿ÷åéêà âñåãäà ñóùåñòâóåò
        const Cell& left_cell = mesh.cells[f.left];
        State W_left = consToPhys(left_cell.U, gamma);

        State W_right;
        if (f.is_boundary()) {
            W_right = ghost_state(W_left, f.nx, f.ny, f.bc, cfg);
        }
        else {
            const Cell& right_cell = mesh.cells[f.right];
            W_right = consToPhys(right_cell.U, gamma);
        }

        // Âû÷èñëÿåì ïîòîê ÷åðåç ãðàíü ñ ó÷¸òîì íîðìàëè
        Flux flux = compute_interface_flux(W_left, W_right, cfg, f.nx, f.ny);

        // Äîáàâëÿåì âêëàä â íåâÿçêè
        mesh.cells[f.left].res.rho -= flux.rho_f * f.length;
        mesh.cells[f.left].res.rhou -= flux.rhou_f * f.length;
        mesh.cells[f.left].res.rhov -= flux.rhov_f * f.length;
        mesh.cells[f.left].res.E -= flux.E_f * f.length;

        if (!f.is_boundary()) {
            mesh.cells[f.right].res.rho += flux.rho_f * f.length;
            mesh.cells[f.right].res.rhou += flux.rhou_f * f.length;
            mesh.cells[f.right].res.rhov += flux.rhov_f * f.length;
            mesh.cells[f.right].res.E += flux.E_f * f.length;
        }
    }
}

// ------------------------------------------------------------
//  compute_dt_unstructured – óñëîâèå Êóðàíòà
// ------------------------------------------------------------
double compute_dt_unstructured(const Mesh& mesh, const Config& cfg) {
    double CFL = cfg.grid.CFL;
    double dt = 1e20;
    const double gamma = cfg.phys.gamma;

    for (const Face& f : mesh.faces) {
        // Ïðîâåðÿåì îáå ÿ÷åéêè, ïðèëåãàþùèå ê ãðàíè
        for (int side = 0; side < 2; ++side) {
            int ci = (side == 0) ? f.left : f.right;
            if (ci < 0) continue;
            const Cell& cell = mesh.cells[ci];
            double rho = cell.U.rho;
            double u = cell.U.rhou / rho;
            double v = cell.U.rhov / rho;
            double p = (gamma - 1.0) * (cell.U.E - 0.5 * rho * (u * u + v * v));
            double a = std::sqrt(gamma * p / rho);
            double un = u * f.nx + v * f.ny;
            double lambda = std::abs(un) + a;
            double h = cell.vol / f.length;   // õàðàêòåðíûé ðàçìåð
            dt = std::min(dt, CFL * h / (lambda + 1e-12));
        }
    }
    return dt;
}

// ------------------------------------------------------------
//  update_solution_unstructured – ïðîñòîé Ýéëåð (äëÿ íà÷àëà)
// ------------------------------------------------------------
void update_solution_unstructured(Mesh& mesh, double dt, const Config& cfg) {
    const double gamma = cfg.phys.gamma;
    for (auto& cell : mesh.cells) {
        cell.U.rho += dt / cell.vol * cell.res.rho;
        cell.U.rhou += dt / cell.vol * cell.res.rhou;
        cell.U.rhov += dt / cell.vol * cell.res.rhov;
        cell.U.E += dt / cell.vol * cell.res.E;

        // Ïðåäîõðàíèòåëè
        if (cell.U.rho < 1e-12) cell.U.rho = 1e-12;
        double rho = cell.U.rho;
        double u = cell.U.rhou / rho;
        double v = cell.U.rhov / rho;
        double p = (gamma - 1.0) * (cell.U.E - 0.5 * rho * (u * u + v * v));
        if (p < 1e-12) p = 1e-12;
        // Ïðè íåîáõîäèìîñòè ìîæíî ïåðåñ÷èòàòü U.E èç p (äëÿ óñòîé÷èâîñòè)
    }
}

// ------------------------------------------------------------
//  time_step_unstructured – îäèí øàã ïî âðåìåíè
// ------------------------------------------------------------
void time_step_unstructured(Mesh& mesh, double dt, const Config& cfg) {
    switch (cfg.time_integrator) {
    case TimeIntegrator::EULER:
        compute_residuals_unstructured(mesh, cfg);
        update_solution_unstructured(mesh, dt, cfg);
        break;
    case TimeIntegrator::TVD_RK3: {
        // Ñîõðàíÿåì íà÷àëüíîå ñîñòîÿíèå
        std::vector<Conserved> U0(mesh.cells.size());
        for (size_t i = 0; i < mesh.cells.size(); ++i)
            U0[i] = mesh.cells[i].U;

        // Ñòàäèÿ 1
        compute_residuals_unstructured(mesh, cfg);
        for (size_t i = 0; i < mesh.cells.size(); ++i) {
            mesh.cells[i].U = U0[i] + (dt / mesh.cells[i].vol) * mesh.cells[i].res;
        }

        // Ñòàäèÿ 2
        compute_residuals_unstructured(mesh, cfg);
        for (size_t i = 0; i < mesh.cells.size(); ++i) {
            Conserved U1 = U0[i] + (dt / mesh.cells[i].vol) * mesh.cells[i].res;
            mesh.cells[i].U = 0.75 * U0[i] + 0.25 * U1;
        }

        // Ñòàäèÿ 3
        compute_residuals_unstructured(mesh, cfg);
        for (size_t i = 0; i < mesh.cells.size(); ++i) {
            Conserved U2 = mesh.cells[i].U;
            mesh.cells[i].U = (1.0 / 3.0) * U0[i] + (2.0 / 3.0) * U2 +
                (2.0 / 3.0) * (dt / mesh.cells[i].vol) * mesh.cells[i].res;
        }
        break;
    }
    default:
        throw std::runtime_error("Unknown time integrator for unstructured mesh");
    }
}

// ------------------------------------------------------------
//  init_mesh – íà÷àëüíûå óñëîâèÿ ïî ÷åòûð¸ì êâàäðàíòàì
// ------------------------------------------------------------
void init_mesh(Mesh& mesh, const Config& cfg) {
    const double x_dia = cfg.grid.x_diaphragm;
    const double y_dia = cfg.grid.y_diaphragm;
    const double gamma = cfg.phys.gamma;

    for (auto& cell : mesh.cells) {
        State init_state;
        if (cell.cx < x_dia && cell.cy < y_dia)
            init_state = cfg.phys.left_bottom;
        else if (cell.cx < x_dia && cell.cy >= y_dia)
            init_state = cfg.phys.left_top;
        else if (cell.cx >= x_dia && cell.cy < y_dia)
            init_state = cfg.phys.right_bottom;
        else
            init_state = cfg.phys.right_top;

        cell.U = physToCons(init_state, gamma);
    }
}
// ------------------------------------------------------------
//  interpolate_to_nodes – âçâåøåííîå óñðåäíåíèå ñ ÿ÷ååê íà óçëû
// ------------------------------------------------------------
static std::vector<State> interpolate_to_nodes(const Mesh& mesh, double gamma) {
    const int n_nodes = mesh.nodes.size();
    std::vector<State> node_state(n_nodes, { 0.0, 0.0, 0.0, 0.0 });
    std::vector<double> weight(n_nodes, 0.0);

    for (const auto& cell : mesh.cells) {
        State W = consToPhys(cell.U, gamma);
        for (int nid : cell.node_ids) {
            node_state[nid].rho += W.rho * cell.vol;
            node_state[nid].u += W.u * cell.vol;
            node_state[nid].v += W.v * cell.vol;
            node_state[nid].p += W.p * cell.vol;
            weight[nid] += cell.vol;
        }
    }

    for (int i = 0; i < n_nodes; ++i) {
        if (weight[i] > 1e-30) {
            node_state[i].rho /= weight[i];
            node_state[i].u /= weight[i];
            node_state[i].v /= weight[i];
            node_state[i].p /= weight[i];
        }
    }
    return node_state;
}

// ------------------------------------------------------------
//  save_snapshot_unstructured – òåïåðü ïèøåò çíà÷åíèÿ â óçëàõ
// ------------------------------------------------------------
void save_snapshot_unstructured(const Mesh& mesh, const Config& cfg,
    int step, double t, const std::string& filename) {
    std::string full = cfg.output.snapshots_directory + "/" + filename;
    std::ofstream f(full);
    if (!f.is_open()) {
        std::cerr << "Cannot open " << full << std::endl;
        return;
    }

    const double gamma = cfg.phys.gamma;
    auto node_state = interpolate_to_nodes(mesh, gamma);

    // Çàãîëîâîê — äîáàâëÿåì n0/n1/n2 äëÿ ïåðåäà÷è òîïîëîãèè â Python
    f << "x,y,rho,u,v,p,E,step,time\n";
    for (int i = 0; i < (int)mesh.nodes.size(); ++i) {
        const auto& nd = mesh.nodes[i];
        const State& W = node_state[i];
        // E ïåðåñ÷èòûâàåì èç óçëîâûõ p/rho/u/v
        double E = W.p / (gamma - 1.0) + 0.5 * W.rho * (W.u * W.u + W.v * W.v);
        f << nd.x << "," << nd.y << ","
            << W.rho << "," << W.u << "," << W.v << "," << W.p << ","
            << E << "," << step << "," << t << "\n";
    }
    f.close();
}
// Âàðèàíò: äâà ôàéëà — nodes + triangles (îäèí ðàç, íå â êàæäîì ñíèìêå)
void save_mesh_topology(const Mesh& mesh, const Config& cfg) {
    std::string full = cfg.output.snapshots_directory + "/mesh_topology.csv";
    std::ofstream f(full);
    f << "n0,n1,n2\n";
    for (const auto& cell : mesh.cells)
        f << cell.node_ids[0] << "," << cell.node_ids[1] << "," << cell.node_ids[2] << "\n";
    f.close();
}
// ------------------------------------------------------------
//  run_unstructured – ãëàâíûé öèêë
// ------------------------------------------------------------
void run_unstructured(Mesh& mesh, const Config& cfg, const std::string& output_filename) {
    init_mesh(mesh, cfg);
    save_mesh_topology(mesh, cfg);

    double t = 0.0;
    int step = 0;

    // Íàñòðîéêè âûâîäà ñíèìêîâ
    int next_snapshot_step = 0;
    double next_snapshot_time = 0.0;
    if (cfg.output.snapshot_output != SnapshotOutputType::NONE) {
        // Ñîõðàíÿåì íà÷àëüíîå ñîñòîÿíèå
        save_snapshot_unstructured(mesh, cfg, 0, 0.0, "initial_state.csv");
        if (cfg.output.snapshot_output == SnapshotOutputType::BY_STEPS)
            next_snapshot_step = cfg.output.snapshot_interval_steps;
        else if (cfg.output.snapshot_output == SnapshotOutputType::BY_TIME)
            next_snapshot_time = cfg.output.snapshot_interval_time;
    }

    while (t < cfg.grid.t_final) {
        double dt = compute_dt_unstructured(mesh, cfg);
        if (t + dt > cfg.grid.t_final) dt = cfg.grid.t_final - t;

        time_step_unstructured(mesh, dt, cfg);
        t += dt;
        ++step;

        // Ñîõðàíåíèå ñíèìêîâ
        bool save = false;
        std::string snap_name;
        if (cfg.output.snapshot_output == SnapshotOutputType::BY_STEPS && step >= next_snapshot_step) {
            save = true;
            snap_name = "snapshot_step_" + std::to_string(step) + ".csv";
            next_snapshot_step += cfg.output.snapshot_interval_steps;
        }
        else if (cfg.output.snapshot_output == SnapshotOutputType::BY_TIME && t >= next_snapshot_time) {
            save = true;
            std::string ts = std::to_string(t);
            std::replace(ts.begin(), ts.end(), '.', '_');
            snap_name = "snapshot_time_" + ts + ".csv";
            next_snapshot_time += cfg.output.snapshot_interval_time;
        }
        if (save) save_snapshot_unstructured(mesh, cfg, step, t, snap_name);

        if (step % 100 == 0) {
            std::cout << "Unstructured step " << step << ", t = " << t << std::endl;
        }
    }

    // Ñîõðàíÿåì ôèíàëüíîå ñîñòîÿíèå
    if (cfg.output.snapshot_output != SnapshotOutputType::NONE)
        save_snapshot_unstructured(mesh, cfg, step, t, "final_state.csv");

    // Âûâîä ðåçóëüòàòà â CSV (ìîæíî ñîõðàíèòü êàê snapshot, íî îñòàâèì äëÿ ñîâìåñòèìîñòè)
    std::ofstream fout(output_filename);
    if (!fout.is_open()) throw std::runtime_error("Cannot open " + output_filename);
    fout << "x,y,rho,u,v,p,E\n";
    for (const auto& cell : mesh.cells) {
        State W = consToPhys(cell.U, cfg.phys.gamma);
        fout << cell.cx << "," << cell.cy << ","
            << W.rho << "," << W.u << "," << W.v << "," << W.p << ","
            << cell.U.E << "\n";
    }
    fout.close();
    std::cout << "Results saved to " << output_filename << std::endl;
}
