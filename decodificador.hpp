#ifndef DECODIFICADOR_HPP
#define DECODIFICADOR_HPP

#include "scp_cs_data.hpp"
#include <vector>
#include <utility>

double decodificar(std::vector<float> genes, const SCPCSInstance& instancia);
std::set<int> busca_local_remocao(std::set<int> solucao_inicial, SCPCSInstance& instancia, double& custo_inicial);
double calcular_custo_solucao(const std::set<int>& subconjuntos_selecionados, const SCPCSInstance& instancia);
std::set<int> decodificar_para_solucao(std::vector<float> genes, const SCPCSInstance& instancia);
#endif // DECODIFICADOR_HPP