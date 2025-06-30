#include <iostream>
#include <cmath>

int main (){
    float EI = 1.;
    float L = 1.;
    float k = 3.*EI/L*L*L;
    float m = 1.;
    float posicao = 0.;
    float velocidade = 0.;
    float forca = 1.;
    float tempo = 20.;
    float flinha = 0.;
    float passo = 0.1;
    float w = sqrt(k/m);
    float aceleracao = 0.;//(forca - k*posicao) / m;
    float klinha = k + ((4*m)/(passo*passo));

    std::cout << "Klinha: " << klinha << '\n';
    std::cout << "Sem amortecimento: " << '\n';   

    for (float i = 0.; i <= tempo; i = i + passo){
        flinha = forca + m*(aceleracao+ 4*velocidade/passo + 4*posicao/(passo*passo));
        float posicaoaux = flinha/klinha;
        float velocidadeaux = -velocidade + 2*(posicaoaux - posicao)/passo;
        float aceleracaoaux = -aceleracao + 4*(posicaoaux - posicao - velocidade*passo)/(passo*passo);
        aceleracao = aceleracaoaux;
        velocidade = velocidadeaux;
        posicao = posicaoaux;
        std::cout << posicao <<  " " << velocidade << " " << aceleracao << " " << flinha << '\n';    
    }

    std::cout << "Com amorteciemnto: " << '\n';   

    posicao = 0.;
    velocidade = 0.;
    float xi = 0.05;
    float c = 2*xi*m*w;
    aceleracao = (forca - k*posicao - c*velocidade) / m;
    klinha = k + ((4*m)/(passo*passo)) + (2*c)/passo;

    for (float i = 0.; i <= tempo; i = i + passo){
        flinha = forca + m*(aceleracao+ 4*velocidade/passo + 4*posicao/(passo*passo)) + c*(velocidade + (2*posicao)/passo);
        float posicaoaux = flinha/klinha;
        float velocidadeaux = -velocidade + 2*(posicaoaux - posicao)/passo;
        float aceleracaoaux = -aceleracao + 4*(posicaoaux - posicao - velocidade*passo)/(passo*passo);
        aceleracao = aceleracaoaux;
        velocidade = velocidadeaux;
        posicao = posicaoaux;
        std::cout << posicao << '\n';    
    }

    return 0;
}