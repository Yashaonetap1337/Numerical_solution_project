// mesh_io.hpp
#pragma once
#include "types.h"
#include <map>
#include <string>

// Загрузка сетки из .msh v2.2. tag_to_bc – соответствие тегов физгрупп типам BC.
Mesh load_msh(const std::string& filename,
    const std::map<int, BoundaryType>& tag_to_bc);

// Построение граней из списка ячеек (треугольников)
void build_faces(Mesh& mesh);

// Вычисление площади и центроидов ячеек
void cell_geometry(Mesh& mesh);