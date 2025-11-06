#include "scp_cs_data.hpp"
#include "decodificador.hpp"
#include <iostream>
#include <string>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <future>
#include <chrono>

#define MIN_VALUE 0.01 //valor minimo para um gene
#define MAX_VALUE 0.99 //valor maximo para um gene
#define RHO 0.7 //probabilidade de um gene de individuo de elute ser escolhido

using namespace std;

struct Cromossomo{
    std::vector<float> genes;
    double fitness = 0;
};

static std::mt19937_64& get_rng() {
    static std::mt19937_64 rng(std::random_device{}());
    return rng;
}

//Gerar número aleatório entre MIN_VALUE e MAX_VALUE
float gerar_numero_aleatorio(){
    std::uniform_real_distribution<double> dist(MIN_VALUE, MAX_VALUE);
    return dist(get_rng());
}

int gerar_indice_aleatorio(int max_exclusive) { //Gera índice entre 0 e max_exclusive-1
    std::uniform_int_distribution<int> dist_int(0, max_exclusive - 1);
    return dist_int(get_rng());
}

//gera uma populacao inicial de cromossomos de maneira aleatoria
vector<Cromossomo> gerar_populacao_inicial(SCPCSInstance instancia, int tamanho_populacao){
    vector<Cromossomo> populacao;

    for(int contador_populacao = 0; contador_populacao < tamanho_populacao; contador_populacao++){
        Cromossomo novo_cromossomo;
        novo_cromossomo.genes.resize(instancia.num_subconjuntos);

        for(int gene_idx = 0; gene_idx < instancia.num_subconjuntos; gene_idx++){
            novo_cromossomo.genes[gene_idx] = gerar_numero_aleatorio();
        }

        populacao.push_back(novo_cromossomo);
    }
    return populacao;
}

void aplicar_fitness_paralela(vector<Cromossomo>& populacao, int indice_inicio_novos, const SCPCSInstance& instancia) {
    //Cada future representa um thread que está calculando um fitness
    std::vector<std::future<double>> futuros_fitness;

    int num_novos = populacao.size() - indice_inicio_novos; //Calcular fitness apenas dos novos
    futuros_fitness.reserve(num_novos); //Reserva espaço
    //lançar as tarefas em paralelo
    for (size_t i = indice_inicio_novos; i < populacao.size(); ++i) {
        //std::async agenda a função 'decodificar' para rodar em um thread separado
        futuros_fitness.push_back(
            std::async(
                std::launch::async, //força a execução em um novo thread
                decodificar,        //a funcao a ser chamada no thread
                populacao[i].genes,   //genes é copiado
                std::cref(instancia)) //instancia é passada por referência constante
        );
    }
    //o programa principal espera cada thread terminar
    int fut_idx = 0; //indice para o vetor de futuros
    for (size_t i = indice_inicio_novos; i < populacao.size(); ++i) {
        //.get() espera o thread [fut_idx] terminar e retorna o 'double' (fitness)
        populacao[i].fitness = futuros_fitness[fut_idx].get();
        fut_idx++;
    }
}

double brkga(SCPCSInstance& instancia, int tamanho_elite, int tamanho_populacao, int num_geracoes, float percentual_mutantes){
    vector<Cromossomo> populacao = gerar_populacao_inicial(instancia, tamanho_populacao);

    aplicar_fitness_paralela(populacao, 0, instancia);
    std::sort(populacao.begin(), populacao.end(), [](const Cromossomo& a, const Cromossomo& b) {
        return a.fitness < b.fitness;
    });

    int cont_geracao = 0;

    auto start = std::chrono::steady_clock::now(); //inicio contagem do tempo

    while(cont_geracao < num_geracoes){
        //selecionar elite
        vector<Cromossomo> elite(populacao.begin(), populacao.begin() + tamanho_elite);

        //copiar elite para nova populacao
        vector<Cromossomo> nova_populacao = elite;

        //salva o índice onde os indivíduos novos (não-elite) comecarao
        int indice_inicio_novos = nova_populacao.size();

        //determinar porcentagem de mutantes na populacao e gerá-los efetivamente
        vector<Cromossomo> mutantes = gerar_populacao_inicial(instancia, (int)(tamanho_populacao * percentual_mutantes));

        nova_populacao.insert(nova_populacao.end(), mutantes.begin(), mutantes.end());

        int tamanho_atual = nova_populacao.size();

        //implementacao do crossover
        for(int i = 0; i < tamanho_populacao - tamanho_atual; i++){
            //pais escolhidos aleatoriamente (pai1 é escolhido dentro do conjunto de elite)
            const Cromossomo pai1 = elite[gerar_indice_aleatorio(tamanho_elite)]; 
            const Cromossomo pai2 = populacao[tamanho_elite + gerar_indice_aleatorio(tamanho_populacao - tamanho_elite)];

            Cromossomo filho;
            filho.genes.resize(instancia.num_subconjuntos);

            for(int j = 0; j < instancia.num_subconjuntos; j++){
                float probabilidade = gerar_numero_aleatorio();

                if(probabilidade <= RHO){ //se o numero aleatorio gerado é menor ou igual a rho, filho herda gene do pai da elite 
                    filho.genes[j] = pai1.genes[j];
                } else {
                    filho.genes[j] = pai2.genes[j];
                }
            }
            nova_populacao.push_back(filho);
        }
        populacao = nova_populacao;

        //chama a funcao paralela, que só avalia os novos (individuos copiados nao precisam ser re-avaliados)
        aplicar_fitness_paralela(populacao, indice_inicio_novos, instancia);

        //ordenar a populacao em ordem decrescente de fitness
        std::sort(populacao.begin(), populacao.end(), [](const Cromossomo& a, const Cromossomo& b) {
            return a.fitness < b.fitness;
        });
        cont_geracao++;
    }
    double melhor_fitness_bruto = populacao[0].fitness;
    cout << "Melhor solucao (Antes da Busca Local): " << melhor_fitness_bruto << std::endl;
    
    //obter o conjunto da melhor solucao encontrada
    std::set<int> melhor_solucao_bruta = decodificar_para_solucao(populacao[0].genes, instancia);

    //aplica a Busca Local (que atualiza o custo por referencia)
    double custo_refinado = melhor_fitness_bruto; //passa o custo bruto como ponto de partida
    std::set<int> solucao_refinada = busca_local_remocao(melhor_solucao_bruta, instancia, custo_refinado);

    cout << "Melhor solucao (Pos Busca Local): " << custo_refinado << std::endl;

    //contagem tempo
    auto end = std::chrono::steady_clock::now();
    auto elapsed = end - start;
    cout << "Tempo total de execucao: "
            << elapsed/std::chrono::seconds(1)
            << " s" << std::endl;

    cout << "Melhor conjunto encontrado e composto pelos subconjuntos: " << endl;
    for(auto elemento : solucao_refinada)
        cout << elemento + 1 << " ";
    cout << endl << endl;

    return custo_refinado; //custo refinado pela busca local
}

int main(){
    SCPCSInstance inst;

    vector<Cromossomo> populacao;

    string nome_arquivo = "instancias//scpcyc08-3.txt"; //nome do arquivo que contem a instancia
    int k_threshold = 1; //valor k (tolerancia de elementos em comum)
    int tamanho_populacao = 128; //numero de individuos da populacao
    int tamanho_elite = 24; //numero de individuos da elite
    int num_testes = 5; //numero de vezes que a instancia sera executada com a configuracao determinada
    int num_geracoes = 500; //numero de geracoes que cada execução terá
    float percentual_mutantes = 0.2; //percentual de mutantes na populacao

    ler_instancia_scpcs(nome_arquivo, inst, k_threshold);
    calcular_custos_conflito(inst, k_threshold);

    cout << "Teste" << endl;
    cout << "instancia: " << nome_arquivo << endl;
    cout << "populacao: " << tamanho_populacao << endl;
    cout << "tamanho elite: " << tamanho_elite << endl;
    cout << "numero de geracoes: " << num_geracoes << endl;
    cout << "k: " << k_threshold << endl;

    double melhor_solucao = std::numeric_limits<int>::max(), media = 0.0, resultado = 0.0; 
    for(int i = 0; i < num_testes; i++){
        resultado = brkga(inst, tamanho_elite, tamanho_populacao, num_geracoes, percentual_mutantes);
        media += resultado;
        if(resultado < melhor_solucao)
            melhor_solucao = resultado;
    }
    cout << "Melhor solucao encontrada entre todas as execucoes: " << melhor_solucao << endl;
    cout << "Media de todas as solucoes encontradas: " << media / num_testes << endl;

    system("pause");
}