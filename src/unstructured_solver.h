#pragma once
#include "types.h"

// Основная функция для неструктурированной сетки
void run_unstructured(Mesh& mesh, const Config& cfg, const std::string& output_filename);

// Вспомогательные функции
void init_mesh(Mesh& mesh, const Config& cfg);
void compute_residuals_unstructured(Mesh& mesh, const Config& cfg);
double compute_dt_unstructured(const Mesh& mesh, const Config& cfg);
void update_solution_unstructured(Mesh& mesh, double dt, const Config& cfg);
void time_step_unstructured(Mesh & mesh, double dt, const Config & cfg);
