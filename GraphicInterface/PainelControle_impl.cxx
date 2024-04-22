#include <iostream>
#include <FL/Fl_File_Chooser.H>
#include "PainelControle.h"

using namespace std;

void PainelControle::iniciar() {
    cout << "Iniciando interface homem-maquina." << endl;
    this->win->show();
    this->apresentarMenuPrincipal();
}

void PainelControle::apresentarMenuImprimirArquivo() {
    cout << "Apresentando o menu de impressao de arquivos." << endl;
    this->grpMenuPrincipal->hide();
    this->grpImprimirArquivo->show();

    // inicializa os campos e op��es com valores padr�o
    this->outArquivo->value("");
    // this->inCopias->value(1);
    this->btColorido->value(1);
    this->btPeB->value(0);
    this->btImprimir->deactivate();
}

void PainelControle::apresentarMenuPrincipal() {
    cout << "Apresentando o menu principal." << endl;
    this->grpMenuPrincipal->show();
    this->grpImprimirArquivo->hide();
}

void PainelControle::imprimir() {
    string cores;

    if (this->btColorido->value() == 1) {
        cores = "cores";
    }
    else {
        cores = "preto e branco";
    }

    cout << "Imprimindo arquivo " << this->outArquivo->value()
        //  << ", em " << this->inCopias->value() << " copias"
         << ", em " << cores << endl;

    this->apresentarMenuPrincipal();
}

void PainelControle::escolherArquivo() {
    const char * arquivo = fl_file_chooser("Arquivo para imprimir", "*.*", "", 0);

    this->outArquivo->value(arquivo);

    if (arquivo == NULL) {
        this->btImprimir->deactivate();
    }
    else {
        this->btImprimir->activate();
    }
}
