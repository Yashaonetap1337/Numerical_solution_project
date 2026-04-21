#include "mpi_utils.h"
#include <algorithm>
#include <iostream>


// CHECK: DECOMPOSITION
// Рекурсивная координатная бисекция (RCB)
static void rcb(std::vector<int>& ids, const Mesh& mesh, int K, int offset, std::vector<int>& part) {
    if (K <= 1 || ids.empty()) {
        for (int id : ids) part[id] = offset;
        return;
    }
    double x_min = 1e9, x_max = -1e9, y_min = 1e9, y_max = -1e9;
    for (int id : ids) {
        x_min = std::min(x_min, mesh.cells[id].cx);
        x_max = std::max(x_max, mesh.cells[id].cx);
        y_min = std::min(y_min, mesh.cells[id].cy);
        y_max = std::max(y_max, mesh.cells[id].cy);
    }

    bool split_x = (x_max - x_min) > (y_max - y_min);
    std::sort(ids.begin(), ids.end(), [&](int a, int b) {
        return split_x ? (mesh.cells[a].cx < mesh.cells[b].cx) : (mesh.cells[a].cy < mesh.cells[b].cy);
        });

    int K_left = K / 2;
    int K_right = K - K_left;
    int split_idx = (ids.size() * K_left) / K;

    std::vector<int> left_ids(ids.begin(), ids.begin() + split_idx);
    std::vector<int> right_ids(ids.begin() + split_idx, ids.end());

    rcb(left_ids, mesh, K_left, offset, part);
    rcb(right_ids, mesh, K_right, offset + K_left, part);
}


void decompose_mesh(Mesh& global_mesh, Mesh& local_mesh, HaloInfo& halo, int rank, int size) {
    int global_n = global_mesh.cells.size();
    std::vector<int> part(global_n, 0);
    std::vector<int> all_ids(global_n);
    for (int i = 0; i < global_n; i++) all_ids[i] = i;

    // 1. Разбиваем сетку (RCB)
    rcb(all_ids, global_mesh, size, 0, part);

    // 2. Определяем, какие глобальные ячейки нам нужны (свои + соседи/ghost)
    std::map<int, int> global_to_local;
    local_mesh.nodes = global_mesh.nodes; // Для простоты узлы копируем все

    // Сначала добавляем СВОИ ячейки
    for (int i = 0; i < global_n; i++) {
        if (part[i] == rank) {
            global_to_local[i] = local_mesh.cells.size();
            local_mesh.cells.push_back(global_mesh.cells[i]);
            local_mesh.cells.back().face_ids.clear(); // очистим, пересоберем ниже
        }
    }
    halo.n_local = local_mesh.cells.size();

    // Затем ищем GHOST-ячейки
    for (const Face& f : global_mesh.faces) {
        if (f.right >= 0) {
            int pL = part[f.left];
            int pR = part[f.right];

            if (pL == rank && pR != rank) {
                if (global_to_local.find(f.right) == global_to_local.end()) {
                    global_to_local[f.right] = local_mesh.cells.size();
                    local_mesh.cells.push_back(global_mesh.cells[f.right]);
                    halo.recv_ids[pR].push_back(global_to_local[f.right]);
                }
                halo.send_ids[pR].push_back(global_to_local[f.left]);
            }
            else if (pR == rank && pL != rank) {
                if (global_to_local.find(f.left) == global_to_local.end()) {
                    global_to_local[f.left] = local_mesh.cells.size();
                    local_mesh.cells.push_back(global_mesh.cells[f.left]);
                    halo.recv_ids[pL].push_back(global_to_local[f.left]);
                }
                halo.send_ids[pL].push_back(global_to_local[f.right]);
            }
        }
    }
    halo.n_ghost = local_mesh.cells.size() - halo.n_local;

    // Убираем дубликаты из send_ids (так как ячейка могла иметь несколько общих граней с одним соседом)
    for (auto& pair : halo.send_ids) {
        std::sort(pair.second.begin(), pair.second.end());
        pair.second.erase(std::unique(pair.second.begin(), pair.second.end()), pair.second.end());
    }

    // 3. Перестраиваем грани для локальной сетки
    for (const Face& f : global_mesh.faces) {
        bool L_in = global_to_local.count(f.left);
        bool R_in = (f.right >= 0) && global_to_local.count(f.right);

        // Нам нужна грань, если хотя бы одна из ее ячеек - НАША (локальная)
        // Грани между двумя ghost-ячейками нам не нужны
        bool L_is_mine = L_in && (global_to_local[f.left] < halo.n_local);
        bool R_is_mine = R_in && (global_to_local[f.right] < halo.n_local);

        if (L_is_mine || R_is_mine) {
            Face local_f = f;
            local_f.left = global_to_local[f.left];
            if (f.right >= 0) {
                local_f.right = global_to_local[f.right];
            }

            int fid = local_mesh.faces.size();
            local_mesh.faces.push_back(local_f);

            if (L_is_mine) local_mesh.cells[local_f.left].face_ids.push_back(fid);
            if (R_is_mine) local_mesh.cells[local_f.right].face_ids.push_back(fid);
        }
    }
}

// CHECK: HALO_EXCHANGE
void exchange_halo(Mesh& local_mesh, const HaloInfo& halo, MPI_Comm comm) {
    const int FIELDS = 4; // rho, rhou, rhov, E

    std::vector<MPI_Request> reqs;
    std::map<int, std::vector<double>> send_bufs;
    std::map<int, std::vector<double>> recv_bufs;

    // 1. Инициируем получение (Irecv)
    for (const auto& pair : halo.recv_ids) {
        int neighbor = pair.first;
        int count = pair.second.size() * FIELDS;
        recv_bufs[neighbor].resize(count);

        MPI_Request req;
        MPI_Irecv(recv_bufs[neighbor].data(), count, MPI_DOUBLE, neighbor, 0, comm, &req);
        reqs.push_back(req);
    }

    // 2. Упаковываем и отправляем (Isend)
    for (const auto& pair : halo.send_ids) {
        int neighbor = pair.first;
        const auto& ids = pair.second;
        send_bufs[neighbor].resize(ids.size() * FIELDS);

        for (size_t i = 0; i < ids.size(); ++i) {
            const Conserved& U = local_mesh.cells[ids[i]].U;
            send_bufs[neighbor][i * FIELDS + 0] = U.rho;
            send_bufs[neighbor][i * FIELDS + 1] = U.rhou;
            send_bufs[neighbor][i * FIELDS + 2] = U.rhov;
            send_bufs[neighbor][i * FIELDS + 3] = U.E;
        }

        MPI_Request req;
        MPI_Isend(send_bufs[neighbor].data(), send_bufs[neighbor].size(), MPI_DOUBLE, neighbor, 0, comm, &req);
        reqs.push_back(req);
    }

    // 3. Ждем завершения всех пересылок
    if (!reqs.empty()) {
        std::vector<MPI_Status> stats(reqs.size());
        MPI_Waitall(reqs.size(), reqs.data(), stats.data());
    }

    // 4. Распаковываем полученные данные в ghost-ячейки
    for (const auto& pair : halo.recv_ids) {
        int neighbor = pair.first;
        const auto& ids = pair.second;
        for (size_t i = 0; i < ids.size(); ++i) {
            Conserved& U = local_mesh.cells[ids[i]].U;
            U.rho = recv_bufs[neighbor][i * FIELDS + 0];
            U.rhou = recv_bufs[neighbor][i * FIELDS + 1];
            U.rhov = recv_bufs[neighbor][i * FIELDS + 2];
            U.E = recv_bufs[neighbor][i * FIELDS + 3];
        }
    }
}