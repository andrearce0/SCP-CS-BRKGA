#include "decodificador.hpp"
#include "scp_cs_data.hpp"
#include <set>
#include <vector>
#include <algorithm>
#include <limits>
#include <numeric>

using namespace std;

//recebe os genes de um cromossomo, constroi uma solucao e retorna o custo dela
double decodificar(std::vector<float> genes, const SCPCSInstance& instancia) {
    int m = instancia.num_elementos;
    int n = instancia.num_subconjuntos;

    //vetor booleano para rastrear elementos cobertos
    std::vector<bool> elementos_cobertos_mask(m, false);
    int elementos_cobertos_count = 0; 

    //rastreia a solução sendo construída
    std::set<int> subconjuntos_selecionados;
    std::vector<bool> ja_processado(n, false); //rastreia subconjuntos já avaliados/usados

    //acumulador de custo (retorno da funcao)
    double custo_total_acumulado = 0.0;

    //gene_prioridades é um vetor de pares (gene, indice subconjunto) 
    std::vector<std::pair<float, int>> gene_prioridades;
    gene_prioridades.reserve(n);
    for (int j = 0; j < n; ++j) {
        gene_prioridades.push_back({genes[j], j});
    }
    //ordenar gene_prioridade em ordem decrescente pela prioridade do subconjunto
    std::sort(gene_prioridades.begin(), gene_prioridades.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });
    
    //tamanho da LCR (ex: 20% da população, no mínimo 1)
    const int TAMANHO_LCR = std::max(1, (int)(n * 0.20)); 

    //loop principal: Continua enquanto a cobertura não for total
    while (elementos_cobertos_count < m) {
        double melhor_metrica = std::numeric_limits<double>::max(); 
        int melhor_indice = -1;
        double custo_efetivo_do_melhor = 0.0; //armazena o custo do vencedor

        //busca na lcr
        for (int i = 0; i < std::min((int)gene_prioridades.size(), TAMANHO_LCR); ++i) {
            int j = gene_prioridades[i].second; //indice do subconjunto
            if (ja_processado[j]) continue;
            //verifica se o subconjunto adiciona algum novo elemento a solucao atual
            int novos = 0;
            for (int e : instancia.matriz_incidencia[j]) {
                if (!elementos_cobertos_mask[e]) {
                    novos++;
                }
            }
            if (novos == 0) continue; //se o subconjunto nao agrega elemento algum, vá para o proximo da lista

            //calcula os custos de penalidade causados se o subconjunto atual for adicionado
            double penalidade_conf = 0.0;
            for (int l : subconjuntos_selecionados) {
                penalidade_conf += instancia.matriz_conflitos[j][l];
            }
            double custo_efetivo = (double)instancia.custos[j] + penalidade_conf;
            double metrica_gulosa = custo_efetivo / novos; //calcula o quociente do subconjunto

            //se o quociente do subconjunto atual é menor que o menor quociente
            //encontrado até agora, o subconjunto atual é o novo melhor candidato
            if (metrica_gulosa < melhor_metrica) {
                melhor_metrica = metrica_gulosa;
                melhor_indice = j;
                custo_efetivo_do_melhor = custo_efetivo; //salva o custo do subconjunto de menor quociente
            }
        }
        //se nenhum subconjunto da lcr é viavel, realizar a busca no restante dos subconjuntos
        if (melhor_indice == -1) {
            for (int i = TAMANHO_LCR; i < n; ++i) {//a busca é feita da mesma maneira que na lcr 
                int j = gene_prioridades[i].second;
                if (ja_processado[j]) continue;
                int novos = 0;

                for (int e : instancia.matriz_incidencia[j]) {
                    if (!elementos_cobertos_mask[e]) novos++;
                }
                if (novos == 0) continue;

                double penalidade_conf = 0.0;
                for (int l : subconjuntos_selecionados) {
                    penalidade_conf += instancia.matriz_conflitos[j][l];
                }
                double custo_efetivo = (double)instancia.custos[j] + penalidade_conf;
                double metrica_gulosa = custo_efetivo / novos;
                if (metrica_gulosa < melhor_metrica) {
                    melhor_metrica = metrica_gulosa;
                    melhor_indice = j;
                    custo_efetivo_do_melhor = custo_efetivo;
                }
            }
        }
        if (melhor_indice == -1) {
            break; //se nao foi encontrado nenhum candidato viavel, interromper a funcao
        }
        //adiciona o custo total do melhor candidato ao custo corrente da solucao 
        custo_total_acumulado += custo_efetivo_do_melhor;

        //adiciona os novos elementos do candidato selecionado a cobertura atual da solucao  
        for (int e : instancia.matriz_incidencia[melhor_indice]) {
            if (!elementos_cobertos_mask[e]) {
                elementos_cobertos_mask[e] = true;
                elementos_cobertos_count++;
            }
        }
        //adiciona à lista de selecionados (para o calculo de conflito da proxima iteracao)
        subconjuntos_selecionados.insert(melhor_indice);
        ja_processado[melhor_indice] = true; 
    }
    return custo_total_acumulado;
}

double calcular_custo_solucao(const std::set<int>& subconjuntos_selecionados, SCPCSInstance& instancia){
    double custo_total = 0.0;

    //vetor de indices
    std::vector<int> selecionados(subconjuntos_selecionados.begin(), subconjuntos_selecionados.end());

    //soma dos custos dos subconjuntos
    for (int j : selecionados) {
        custo_total += instancia.custos[j];
    }

    //soma das penalidades de conflito
    //como a matriz de conflitos é simetrica, não precisa ser percorrida por completo
    for (size_t i = 0; i < selecionados.size(); ++i) {
        int sub_i = selecionados[i];
        const auto& linha_conflito = instancia.matriz_conflitos[sub_i];

        for (size_t j = i + 1; j < selecionados.size(); ++j) {
            int sub_j = selecionados[j];
            custo_total += linha_conflito[sub_j];
        }
    }
    return custo_total;
}

//a busca local é aplicada somente uma vez: no melhor individuo da ultima populacao
//ela serve para remover subconjuntos redundantes, que nao possuem, exclusivamente, nenhum elemento 
std::set<int> busca_local_remocao(std::set<int> solucao_inicial, SCPCSInstance& instancia, double& custo_inicial) {
    std::set<int> solucao_atual = solucao_inicial;
    bool mudanca_feita = true;
    double custo_atual = custo_inicial;

    //converte o set para vector e ordena pelo custo (tentar remover o mais caro primeiro)
    std::vector<int> candidatos(solucao_inicial.begin(), solucao_inicial.end());
    //ordena de forma decrescente
    std::sort(candidatos.begin(), candidatos.end(), 
        [&](int a, int b) {
            return instancia.custos[a] > instancia.custos[b];
        });

    while(mudanca_feita){
        mudanca_feita = false;

        //tentamos remover cada subconjunto (do mais caro para o mais barato)
        for(int indice_sub_j : candidatos){
            //verifica se o subconjunto j ainda está na solução
            if (solucao_atual.find(indice_sub_j) == solucao_atual.end()) {
                continue; //ja foi removido anteriormente
            }
            bool pode_remover = true;

            //para cada elemento i coberto por j...
            for(int elemento_i : instancia.matriz_incidencia[indice_sub_j]){
                //verificador de cobertura
                bool coberto_por_substituto = false;

                //itera sobre subconjuntos k que cobrem i (usando lista_incidencia)
                for(int subconjunto_k : instancia.lista_incidencia[elemento_i]){
                    //se K não for o próprio J E K estiver na solução, achar um substituto
                    if (subconjunto_k != indice_sub_j && solucao_atual.count(subconjunto_k)){
                        coberto_por_substituto = true;
                        break; 
                    }
                }
                //se nao for coberto por nenhum substituto, j é essencial
                if(!coberto_por_substituto){
                    pode_remover = false;
                    break; //'j' é essencial para manter a cobertura total, não pode ser removido
                }

            }

            if(pode_remover){
                double custo_remocao_delta = instancia.custos[indice_sub_j];
                //subtrair as penalidades de conflito perdidas
                for (int sub_k : solucao_atual) {
                    if (sub_k != indice_sub_j) {
                        custo_remocao_delta += instancia.matriz_conflitos[indice_sub_j][sub_k];
                    }
                }
                //o movimento é sempre de melhoria, pois estamos só removendo.
                
                //aplicar remocao
                solucao_atual.erase(indice_sub_j);
                custo_atual -= custo_remocao_delta; //atualiza o custo global
                mudanca_feita = true; 

                break; //recomeca o while para re-avaliar a nova solução
            } 
        }
    }

    custo_inicial = custo_atual;

    return solucao_atual;

}

//versao alternativa da funcao decodificar, que retorna o conjunto de subconjuntos selecionados
//utilizada para exibir quais os subconjuntos da melhor solucao encontrada ao final
std::set<int> decodificar_para_solucao(std::vector<float> genes, const SCPCSInstance& instancia) {
    int m = instancia.num_elementos;
    int n = instancia.num_subconjuntos;

    std::vector<bool> elementos_cobertos_mask(m, false);
    int elementos_cobertos_count = 0; 
    std::set<int> subconjuntos_selecionados;
    std::vector<bool> ja_processado(n, false); 

    std::vector<std::pair<float, int>> gene_prioridades;
    gene_prioridades.reserve(n);
    for (int j = 0; j < n; ++j) {
        gene_prioridades.push_back({genes[j], j});
    }
    std::sort(gene_prioridades.begin(), gene_prioridades.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    const int TAMANHO_LCR = std::max(1, (int)(n * 0.20)); 

    while (elementos_cobertos_count < m) {
        double melhor_metrica = std::numeric_limits<double>::max(); 
        int melhor_indice = -1;
        //busca na LCR
        for (int i = 0; i < std::min((int)gene_prioridades.size(), TAMANHO_LCR); ++i) {
            int j = gene_prioridades[i].second; 
            if (ja_processado[j]) continue;

            int novos = 0;
            for (int e : instancia.matriz_incidencia[j]) {
                if (!elementos_cobertos_mask[e]) novos++;
            }
            if (novos == 0) continue; 

            double penalidade_conf = 0.0;
            for (int l : subconjuntos_selecionados) {
                penalidade_conf += instancia.matriz_conflitos[j][l];
            }
            double custo_efetivo = (double)instancia.custos[j] + penalidade_conf;
            double metrica_gulosa = custo_efetivo / novos;

            if (metrica_gulosa < melhor_metrica) {
                melhor_metrica = metrica_gulosa;
                melhor_indice = j;
            }
        } 
        //busca fora da lcr
        if (melhor_indice == -1) {
            for (int i = TAMANHO_LCR; i < n; ++i) { 
                int j = gene_prioridades[i].second;
                if (ja_processado[j]) continue;

                int novos = 0;
                for (int e : instancia.matriz_incidencia[j]) {
                    if (!elementos_cobertos_mask[e]) novos++;
                }
                if (novos == 0) continue;

                double penalidade_conf = 0.0;
                for (int l : subconjuntos_selecionados) {
                    penalidade_conf += instancia.matriz_conflitos[j][l];
                }
                double custo_efetivo = (double)instancia.custos[j] + penalidade_conf;
                double metrica_gulosa = custo_efetivo / novos;

                if (metrica_gulosa < melhor_metrica) {
                    melhor_metrica = metrica_gulosa;
                    melhor_indice = j;
                }
            }
        }
        if (melhor_indice == -1) break; 
         
        for (int e : instancia.matriz_incidencia[melhor_indice]) {
            if (!elementos_cobertos_mask[e]) {
                elementos_cobertos_mask[e] = true;
                elementos_cobertos_count++;
            }
        }
        subconjuntos_selecionados.insert(melhor_indice);
        ja_processado[melhor_indice] = true; 
    } 
    //retorna o conjunto dos subconjuntos selecionados
    return subconjuntos_selecionados; 
}