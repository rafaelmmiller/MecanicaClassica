/******************************************************
*            PGF 5005 - Mecânica Clássica             *
*              Primeiro Estudo Dirigido               *
*               2º Semestre de 2015                   *
*                                                     *
*               http://goo.gl/7PvBxi                  *
* Programa 1 - Equações de Euler para pêndulo simples *
*                                                     *
* Aluno: Rafael Mendonça Miller   NUSP.:7581818       *
*         e-mail: rafael.miller@usp.br                *
*******************************************************/

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

double hamiltoniana(double p, double q) {
  return (p*p)/2 - cos(q);
}

double dHamiltoniana_dp(double p) {
  return p;
}

double dHamiltoniana_dq(double q) {
  return -sin(q);
}

double q_nMaisUm(double q_n, double deltaT, double p_n) {
  return q_n + deltaT*dHamiltoniana_dp(p_n);
}

double p_nMaisUm(double p_n, double deltaT, double q_n) {
  return p_n - deltaT*dHamiltoniana_dq(q_n);
}

int main() {

  ofstream dadosLibracao_q_01, dadosLibracao_p_01, dadosLibracao_e_01;
  ofstream dadosLibracao_q_001, dadosLibracao_p_001, dadosLibracao_e_001;
  ofstream dadosLibracao_q_0001, dadosLibracao_p_0001, dadosLibracao_e_0001;
  ofstream dadosLibracao_pvq_01, dadosLibracao_pvq_001, dadosLibracao_pvq_0001;

  dadosLibracao_q_01.open("out/dadosLibracao_q_01.dat"); dadosLibracao_q_001.open("out/dadosLibracao_q_001.dat"), dadosLibracao_q_0001.open("out/dadosLibracao_q_0001.dat");
  dadosLibracao_p_01.open("out/dadosLibracao_p_01.dat"); dadosLibracao_p_001.open("out/dadosLibracao_p_001.dat"), dadosLibracao_p_0001.open("out/dadosLibracao_p_0001.dat");
  dadosLibracao_e_01.open("out/dadosLibracao_e_01.dat"); dadosLibracao_e_001.open("out/dadosLibracao_e_001.dat"), dadosLibracao_e_0001.open("out/dadosLibracao_e_0001.dat");
  dadosLibracao_pvq_01.open("out/dadosLibracao_pvq_01.dat"), dadosLibracao_pvq_001.open("out/dadosLibracao_pvq_001.dat"), dadosLibracao_pvq_0001.open("out/dadosLibracao_pvq_0001.dat");

  double T = 0.0;
  double deltaT = 0.1;

  vector<double> Q_n (99999);
  vector<double> P_n (99999);
  vector<double> E_n (99999);

  // Condições Iniciais para Libração:
  P_n[0] = 0.25;
  Q_n[0] = 1.0;
  E_n[0] = hamiltoniana(P_n[0], Q_n[0]);

  for(int i = 0; T < 25; i++) {

    Q_n[i+1] = q_nMaisUm(Q_n[i], deltaT, P_n[i]);
    P_n[i+1] = p_nMaisUm(P_n[i], deltaT, Q_n[i]);
    E_n[i+1] = hamiltoniana(P_n[i], Q_n[i]);

    // Imprime os dados nos arquivos
    dadosLibracao_q_01 << T << "\t" << Q_n[i] << endl;
    dadosLibracao_p_01 << T << "\t" << P_n[i] << endl;
    dadosLibracao_e_01 << T << "\t" << E_n[i] << endl;
    dadosLibracao_pvq_01 << Q_n[i] << "\t" << P_n[i] << endl;

    T+=deltaT;
  }

  // Diminuindo deltaT e reajustando as Condições Iniciais:
  T = 0.0;
  deltaT = 0.01;

  P_n[0] = 0.25;
  Q_n[0] = 1.0;
  E_n[0] = hamiltoniana(P_n[0], Q_n[0]);

  for(int i = 0; T < 25; i++) {

    Q_n[i+1] = q_nMaisUm(Q_n[i], deltaT, P_n[i]);
    P_n[i+1] = p_nMaisUm(P_n[i], deltaT, Q_n[i]);
    E_n[i+1] = hamiltoniana(P_n[i], Q_n[i]);

    dadosLibracao_q_001 << T << "\t" << Q_n[i] << endl;
    dadosLibracao_p_001 << T << "\t" << P_n[i] << endl;
    dadosLibracao_e_001 << T << "\t" << E_n[i] << endl;
    dadosLibracao_pvq_001 << Q_n[i] << "\t" << P_n[i] << endl;

    T+=deltaT;
  }

  // Diminuindo deltaT e reajustando as Condições Iniciais:
  T = 0.0;
  deltaT = 0.001;

  P_n[0] = 0.25;
  Q_n[0] = 1.0;
  E_n[0] = hamiltoniana(P_n[0], Q_n[0]);

  for(int i = 0; T < 25; i++) {

    Q_n[i+1] = q_nMaisUm(Q_n[i], deltaT, P_n[i]);
    P_n[i+1] = p_nMaisUm(P_n[i], deltaT, Q_n[i]);
    E_n[i+1] = hamiltoniana(P_n[i], Q_n[i]);

    dadosLibracao_q_0001 << T << "\t" << Q_n[i] << endl;
    dadosLibracao_p_0001 << T << "\t" << P_n[i] << endl;
    dadosLibracao_e_0001 << T << "\t" << E_n[i] << endl;
    dadosLibracao_pvq_0001 << Q_n[i] << "\t" << P_n[i] << endl;

    T+=deltaT;
  }

  // Agora redefinimos os arquivos e ajustamos parâmetros para a Rotação

  ofstream dadosRotacao_q_01, dadosRotacao_p_01, dadosRotacao_e_01;
  ofstream dadosRotacao_q_001, dadosRotacao_p_001, dadosRotacao_e_001;
  ofstream dadosRotacao_q_0001, dadosRotacao_p_0001, dadosRotacao_e_0001;
  ofstream dadosRotacao_pvq_01, dadosRotacao_pvq_001, dadosRotacao_pvq_0001;

  dadosRotacao_q_01.open("out/dadosRotacao_q_01.dat"); dadosRotacao_q_001.open("out/dadosRotacao_q_001.dat"), dadosRotacao_q_0001.open("out/dadosRotacao_q_0001.dat");
  dadosRotacao_p_01.open("out/dadosRotacao_p_01.dat"); dadosRotacao_p_001.open("out/dadosRotacao_p_001.dat"), dadosRotacao_p_0001.open("out/dadosRotacao_p_0001.dat");
  dadosRotacao_e_01.open("out/dadosRotacao_e_01.dat"); dadosRotacao_e_001.open("out/dadosRotacao_e_001.dat"), dadosRotacao_e_0001.open("out/dadosRotacao_e_0001.dat");
  dadosRotacao_pvq_01.open("out/dadosRotacao_pvq_01.dat"), dadosRotacao_pvq_001.open("out/dadosRotacao_pvq_001.dat"), dadosRotacao_pvq_0001.open("out/dadosRotacao_pvq_0001.dat");


  // Condições Iniciais para a Rotação:
  T = 0.0;
  deltaT = 0.1;

  P_n[0] = 3.0;
  Q_n[0] = 0.0;
  E_n[0] = hamiltoniana(P_n[0], Q_n[0]);

  for(int i = 0; T < 25; i++) {

    Q_n[i+1] = q_nMaisUm(Q_n[i], deltaT, P_n[i]);
    P_n[i+1] = p_nMaisUm(P_n[i], deltaT, Q_n[i]);
    E_n[i+1] = hamiltoniana(P_n[i], Q_n[i]);

    // Imprime os dados nos arquivos
    dadosRotacao_q_01 << T << "\t" << Q_n[i] << endl;
    dadosRotacao_p_01 << T << "\t" << P_n[i] << endl;
    dadosRotacao_e_01 << T << "\t" << E_n[i] << endl;
    dadosRotacao_pvq_01 << Q_n[i] << "\t" << P_n[i] << endl;

    T+=deltaT;
  }

  // Diminuindo deltaT e reajustando as Condições Iniciais:
  T = 0.0;
  deltaT = 0.01;

  P_n[0] = 3.0;
  Q_n[0] = 0.0;
  E_n[0] = hamiltoniana(P_n[0], Q_n[0]);

  for(int i = 0; T < 25; i++) {

    Q_n[i+1] = q_nMaisUm(Q_n[i], deltaT, P_n[i]);
    P_n[i+1] = p_nMaisUm(P_n[i], deltaT, Q_n[i]);
    E_n[i+1] = hamiltoniana(P_n[i], Q_n[i]);

    dadosRotacao_q_001 << T << "\t" << Q_n[i] << endl;
    dadosRotacao_p_001 << T << "\t" << P_n[i] << endl;
    dadosRotacao_e_001 << T << "\t" << E_n[i] << endl;
    dadosRotacao_pvq_001 << Q_n[i] << "\t" << P_n[i] << endl;

    T+=deltaT;
  }

  // Diminuindo deltaT e reajustando as Condições Iniciais:
  T = 0.0;
  deltaT = 0.001;

  P_n[0] = 3.0;
  Q_n[0] = 0.0;
  E_n[0] = hamiltoniana(P_n[0], Q_n[0]);

  for(int i = 0; T < 25; i++) {

    Q_n[i+1] = q_nMaisUm(Q_n[i], deltaT, P_n[i]);
    P_n[i+1] = p_nMaisUm(P_n[i], deltaT, Q_n[i]);
    E_n[i+1] = hamiltoniana(P_n[i], Q_n[i]);

    dadosRotacao_q_0001 << T << "\t" << Q_n[i] << endl;
    dadosRotacao_p_0001 << T << "\t" << P_n[i] << endl;
    dadosRotacao_e_0001 << T << "\t" << E_n[i] << endl;
    dadosRotacao_pvq_0001 << Q_n[i] << "\t" << P_n[i] << endl;

    T+=deltaT;
  }




  return 0;
}
