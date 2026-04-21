// mesh_io.cpp
#include "mesh_io.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <iostream>

Mesh load_msh(const std::string& filename,
    const std::map<int, BoundaryType>& tag_to_bc) {
    std::ifstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("Cannot open " + filename);

    Mesh mesh;
    std::string line;

    // Óáèðàåì \r äëÿ ôàéëîâ â Windows-ôîðìàòå (CRLF)
    auto getline_clean = [&](std::string& s) -> bool {
        if (!std::getline(f, s)) return false;
        if (!s.empty() && s.back() == '\r') s.pop_back();
        return true;
        };
    auto seek_section = [&](const std::string& sec) {
        while (getline_clean(line) && line != sec) {}
        };

    // -------- $Nodes --------
    seek_section("$Nodes");
    int nnodes; f >> nnodes;
    mesh.nodes.resize(nnodes);
    for (int i = 0; i < nnodes; ++i) {
        int id; double x, y, z;
        f >> id >> x >> y >> z;
        mesh.nodes[id - 1] = { x, y };
    }

    // -------- $Elements --------
    seek_section("$Elements");
    int nelem; f >> nelem;

    std::map<std::pair<int, int>, int> edge_phys;

    // Òàáëèöà ÷èñëà óçëîâ äëÿ òèïîâ ýëåìåíòîâ MSH v2
    static const std::map<int, int> type_nodes = {
        {1,2},{2,3},{3,4},{4,4},{5,8},{6,6},{7,5},
        {8,3},{9,6},{10,9},{11,10},{15,1},{16,8},{17,20}
    };

    for (int i = 0; i < nelem; ++i) {
        int id, type, ntags;
        f >> id >> type >> ntags;
        int phys_tag = 0;
        for (int t = 0; t < ntags; ++t) {
            int tag; f >> tag;
            if (t == 0) phys_tag = tag;
        }
        if (type == 1) {
            int n0, n1; f >> n0 >> n1;
            edge_phys[{std::min(n0, n1) - 1, std::max(n0, n1) - 1}] = phys_tag;
        }
        else if (type == 2) {
            int n0, n1, n2; f >> n0 >> n1 >> n2;
            mesh.cells.push_back({ {n0 - 1, n1 - 1, n2 - 1} });
        }
        else {
            auto it = type_nodes.find(type);
            int nn = (it != type_nodes.end()) ? it->second : 0;
            for (int k = 0; k < nn; ++k) { int dummy; f >> dummy; }
        }
    }

    cell_geometry(mesh);
    build_faces(mesh);

    for (Face& face : mesh.faces) {
        if (face.is_boundary()) {
            auto it = edge_phys.find({ std::min(face.n0,face.n1),
                                      std::max(face.n0,face.n1) });
            if (it != edge_phys.end()) {
                auto jt = tag_to_bc.find(it->second);
                face.bc = (jt != tag_to_bc.end()) ? jt->second : BoundaryType::WALL;
            }
            else {
                face.bc = BoundaryType::WALL;
            }
        }
    }
    return mesh;
}
//CHECK: BUILD_FACES
void cell_geometry(Mesh& mesh) {
    for (auto& cell : mesh.cells) {
        const Node& n0 = mesh.nodes[cell.node_ids[0]];
        const Node& n1 = mesh.nodes[cell.node_ids[1]];
        const Node& n2 = mesh.nodes[cell.node_ids[2]];
        double area2 = (n1.x - n0.x) * (n2.y - n0.y) -
            (n2.x - n0.x) * (n1.y - n0.y);
        cell.vol = 0.5 * std::abs(area2);
        cell.cx = (n0.x + n1.x + n2.x) / 3.0;
        cell.cy = (n0.y + n1.y + n2.y) / 3.0;
    }
}

void build_faces(Mesh& mesh) {
    using Edge = std::pair<int, int>;
    std::map<Edge, std::vector<int>> edge_to_cells;

    // Ñîáèðàåì ð¸áðà âñåõ òðåóãîëüíèêîâ
    for (size_t ci = 0; ci < mesh.cells.size(); ++ci) {
        const auto& cell = mesh.cells[ci];
        int n0 = cell.node_ids[0];
        int n1 = cell.node_ids[1];
        int n2 = cell.node_ids[2];
        Edge e1 = { std::min(n0, n1), std::max(n0, n1) };
        Edge e2 = { std::min(n1, n2), std::max(n1, n2) };
        Edge e3 = { std::min(n2, n0), std::max(n2, n0) };
        edge_to_cells[e1].push_back(ci);
        edge_to_cells[e2].push_back(ci);
        edge_to_cells[e3].push_back(ci);
    }

    mesh.faces.clear();
    for (const auto& kv : edge_to_cells) {
        const Edge& e = kv.first;
        const std::vector<int>& owners = kv.second;
        Face f;
        f.n0 = e.first;
        f.n1 = e.second;
        // Ãåîìåòðèÿ ãðàíè
        double x0 = mesh.nodes[f.n0].x, y0 = mesh.nodes[f.n0].y;
        double x1 = mesh.nodes[f.n1].x, y1 = mesh.nodes[f.n1].y;
        double dx = x1 - x0, dy = y1 - y0;
        f.length = std::hypot(dx, dy);
        f.mx = (x0 + x1) * 0.5;
        f.my = (y0 + y1) * 0.5;
        // Íà÷àëüíàÿ íîðìàëü (ëåâàÿ)
        f.nx = dy / f.length;
        f.ny = -dx / f.length;

        if (owners.size() == 2) {
            f.left = owners[0];
            f.right = owners[1];
            // Êîððåêòèðóåì íîðìàëü: îíà äîëæíà ñìîòðåòü îò left ê right
            double dcx = mesh.cells[f.right].cx - mesh.cells[f.left].cx;
            double dcy = mesh.cells[f.right].cy - mesh.cells[f.left].cy;
            if (f.nx * dcx + f.ny * dcy < 0) {
                f.nx = -f.nx;
                f.ny = -f.ny;
            }
        }
        else {
            f.left = owners[0];
            f.right = -1;
            // Äëÿ ãðàíèöû íîðìàëü íàïðàâëÿåì íàðóæó
            double cx = mesh.cells[f.left].cx;
            double cy = mesh.cells[f.left].cy;
            double dot = f.nx * (f.mx - cx) + f.ny * (f.my - cy);
            if (dot < 0) {
                f.nx = -f.nx;
                f.ny = -f.ny;
            }
        }

        int fid = mesh.faces.size();
        mesh.faces.push_back(f);
        mesh.cells[f.left].face_ids.push_back(fid);
        if (f.right >= 0)
            mesh.cells[f.right].face_ids.push_back(fid);
    }
}
