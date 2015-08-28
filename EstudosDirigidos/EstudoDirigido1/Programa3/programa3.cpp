/******************************************************************************
*                            PGF 5005 - Mecânica Clássica                     *
*                                Primeiro Estudo Dirigido                     *
*                                  2º Semestre de 2015                        *
*                                                                             *
*                                   http://goo.gl/7PvBxi                      *
* Programa 3 - Método de Euler simplético para a Hamiltoniana de Hénon-Helies *
*                                                                             *
*                   Aluno: Rafael Mendonça Miller   NUSP.:7581818             *
*                           e-mail: rafael.miller@usp.br                      *
*******************************************************************************/

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

double hamiltoniana(double p1, double p2, double q1, double q2) {
  return (1./2)*(p1*p1 + p2*p2 + q1*q1 + q2*q2) + q1*q1*q2 - (1./3)*q2*q2*q2;
}

double dHamiltoniana_dp1(double p1) {
  return p1;
}

double dHamiltoniana_dp2(double p2) {
  return p2;
}

double dHamiltoniana_dq1(double q1, double q2) {
  return q1 + 2*q1*q2;
}

double dHamiltoniana_dq2(double q1, double q2) {
  return q2 + q1*q1 - q2*q2;
}

double q1_nMaisUm(double q1_n, double deltaT, double p1_n) {
  return q1_n + deltaT*dHamiltoniana_dp1(p1_n);
}

double q2_nMaisUm(double q2_n, double deltaT, double p2_n) {
  return q2_n + deltaT*dHamiltoniana_dp2(p2_n);
}

double p1_nMaisUm(double p1_n, double deltaT, double q1_n, double q2_n) {
  return p1_n - deltaT*dHamiltoniana_dq1(q1_n, q2_n);
}

double p2_nMaisUm(double p2_n, double deltaT, double q1_n, double q2_n) {
  return p2_n - deltaT*dHamiltoniana_dq2(q1_n, q2_n);
}

double p2s(double p1_n, double p2_n, double q1_n, double q2_n, double DeltaQ1) {
    return p2_n - DeltaQ1*(dHamiltoniana_dq2(q1_n, q2_n)/dHamiltoniana_dp1(p1_n));
}

double q2s(double p1_n, double p2_n, double q2_n, double DeltaQ1) {
    return q2_n + DeltaQ1*(dHamiltoniana_dp2(p2_n)/dHamiltoniana_dp1(p1_n));
}

int main() {

  double deltaT = 0.0001;

  double Q1_n;
  double Q1_nMaisUm;

  double Q2_n;
  double Q2_nMaisUm;

  double P1_n;
  double P1_nMaisUm;

  double P2_n;
  double P2_nMaisUm;

  double E_n;

  double Q2s;
  double P2s;


  /***************************** E = 0.08333: *********************************/

  vector <double> CondicaoInicial_Q2 (30, 0.0);
  vector <double> CondicaoInicial_P2 (30, 0.0);

  CondicaoInicial_P2[0] = 0.40824; // Separatriz
  CondicaoInicial_P2[1] = 0.1;
  CondicaoInicial_P2[2] = 0.27515; // Separatriz
  CondicaoInicial_P2[3] = 0.22;
  CondicaoInicial_P2[4] = -0.22;

  CondicaoInicial_Q2[5] = 0.35;
  CondicaoInicial_Q2[6] = -0.2;
  CondicaoInicial_Q2[7] = 0.3;
  CondicaoInicial_Q2[8] = -0.29;


  for(int j = 0; j < 9; j ++) {

    ofstream dadosMapaDePoincare;
    string file_name = "out/dadosMapaDePoincare";
    file_name += to_string(j);
    file_name += ".dat";
    dadosMapaDePoincare.open(file_name);

    Q1_n = 0.0;

    P2_n = CondicaoInicial_P2[j];
    Q2_n = CondicaoInicial_Q2[j];

    cout << P2_n << "\t" << Q2_n << "\t" << file_name << endl;
    E_n = 0.08333;
    P1_n  = sqrt(2*(E_n - P2_n*P2_n*(1./2) - Q1_n*Q1_n*(1./2) - Q2_n*Q2_n*(1./2) - Q1_n*Q1_n*Q2_n + (1./3)*Q2_n*Q2_n*Q2_n));

    for(double T = 0.0; T < 20000; T+=deltaT) {

      Q1_nMaisUm = q1_nMaisUm(Q1_n, deltaT, P1_n);
      Q2_nMaisUm = q2_nMaisUm(Q2_n, deltaT, P2_n);

      P1_nMaisUm = p1_nMaisUm(P1_n, deltaT, Q1_nMaisUm, Q2_nMaisUm);
      P2_nMaisUm = p2_nMaisUm(P2_n, deltaT, Q1_nMaisUm, Q2_nMaisUm);

      if(Q1_nMaisUm*Q1_n < 0 && P1_n > 0) {

        Q2s = q2s(P1_n, P2_n, Q2_n, -Q1_n);
        P2s = p2s(P1_n, P2_n, Q1_n, Q2_n, -Q1_n);

        dadosMapaDePoincare << Q2s << "\t" << P2s << endl;
      }

      Q1_n = Q1_nMaisUm;
      Q2_n = Q2_nMaisUm;
      P1_n = P1_nMaisUm;
      P2_n = P2_nMaisUm;
    }

  }

  /***************************** E = 0.125: ***********************************/

  ofstream q1PeriodicaE0125,
           p1PeriodicaE0125,
           q2PeriodicaE0125,
           p2PeriodicaE0125;

  q1PeriodicaE0125.open("out/q1PeriodicaE0125.dat");
  p1PeriodicaE0125.open("out/p1PeriodicaE0125.dat");
  q2PeriodicaE0125.open("out/q2PeriodicaE0125.dat");
  p2PeriodicaE0125.open("out/p2PeriodicaE0125.dat");

  ofstream q1CaoticaE0125,
           p1CaoticaE0125,
           q2CaoticaE0125,
           p2CaoticaE0125;

  q1CaoticaE0125.open("out/q1CaoticaE0125.dat");
  p1CaoticaE0125.open("out/p1CaoticaE0125.dat");
  q2CaoticaE0125.open("out/q2CaoticaE0125.dat");
  p2CaoticaE0125.open("out/p2CaoticaE0125.dat");

  CondicaoInicial_P2[9] = 0.40824;
  CondicaoInicial_P2[10] = 0.1; // Quasi-Periódica ou Separatriz?
  CondicaoInicial_P2[11] = 0.22;
  CondicaoInicial_P2[12] = -0.22;

  CondicaoInicial_Q2[13] = 0.4; // Caótico
  CondicaoInicial_Q2[14] = -0.2;
  CondicaoInicial_Q2[15] = 0.3; // Centro da ilha
  CondicaoInicial_Q2[16] = 0.5;
  CondicaoInicial_Q2[17] = 0.56;

  for(int j = 9; j < 18; j ++) {

    ofstream dadosMapaDePoincare;
    string file_name = "out/dadosMapaDePoincare";
    file_name += to_string(j);
    file_name += ".dat";
    dadosMapaDePoincare.open(file_name);

    Q1_n = 0.0;

    P2_n = CondicaoInicial_P2[j];
    Q2_n = CondicaoInicial_Q2[j];

    cout << P2_n << "\t" << Q2_n << "\t" << file_name << endl;
    E_n = 0.125;
    P1_n  = sqrt(2*(E_n - P2_n*P2_n*(1./2) - Q1_n*Q1_n*(1./2) - Q2_n*Q2_n*(1./2) - Q1_n*Q1_n*Q2_n + (1./3)*Q2_n*Q2_n*Q2_n));

    for(double T = 0.0; T < 20000; T+=deltaT) {

      Q1_nMaisUm = q1_nMaisUm(Q1_n, deltaT, P1_n);
      Q2_nMaisUm = q2_nMaisUm(Q2_n, deltaT, P2_n);

      P1_nMaisUm = p1_nMaisUm(P1_n, deltaT, Q1_nMaisUm, Q2_nMaisUm);
      P2_nMaisUm = p2_nMaisUm(P2_n, deltaT, Q1_nMaisUm, Q2_nMaisUm);

      if(Q1_nMaisUm*Q1_n < 0 && P1_n > 0) {

        Q2s = q2s(P1_n, P2_n, Q2_n, -Q1_n);
        P2s = p2s(P1_n, P2_n, Q1_n, Q2_n, -Q1_n);

        dadosMapaDePoincare << Q2s << "\t" << P2s << endl;

        // Para construir os gráficos de q1 x t, q2 x t, p1 x t, p2 x t:
        if(j == 11) {
          q1PeriodicaE0125 << T << "\t" << Q1_nMaisUm << endl;
          p1PeriodicaE0125 << T << "\t" << P1_nMaisUm << endl;
          q2PeriodicaE0125 << T << "\t" << Q2_nMaisUm  << endl;
          p2PeriodicaE0125 << T << "\t" << P2_nMaisUm << endl;
        }

        if(j == 13) {
          q1CaoticaE0125 << T << "\t" << Q1_nMaisUm << endl;
          p1CaoticaE0125 << T << "\t" << P1_nMaisUm << endl;
          q2CaoticaE0125 << T << "\t" << Q2_nMaisUm << endl;
          p2CaoticaE0125 << T << "\t" << P2_nMaisUm << endl;
        }
      }

      Q1_n = Q1_nMaisUm;
      Q2_n = Q2_nMaisUm;
      P1_n = P1_nMaisUm;
      P2_n = P2_nMaisUm;
    }

  }

  /**************************** E = 0.16667: ***********************************/

  CondicaoInicial_P2[18] = 0.1;
  CondicaoInicial_P2[19] = 0.22;
  CondicaoInicial_P2[20] = -0.22;

  CondicaoInicial_Q2[21] = -0.2;
  CondicaoInicial_Q2[22] = 0.3;
  CondicaoInicial_Q2[23] = 0.5;

  for(int j = 18; j < 24; j ++) {

    ofstream dadosMapaDePoincare;
    string file_name = "out/dadosMapaDePoincare";
    file_name += to_string(j);
    file_name += ".dat";
    dadosMapaDePoincare.open(file_name);

    Q1_n = 0.0;

    P2_n = CondicaoInicial_P2[j];
    Q2_n = CondicaoInicial_Q2[j];

    cout << P2_n << "\t" << Q2_n << "\t" << file_name << endl;
    E_n = 0.16667;
    P1_n  = sqrt(2*(E_n - P2_n*P2_n*(1./2) - Q1_n*Q1_n*(1./2) - Q2_n*Q2_n*(1./2) - Q1_n*Q1_n*Q2_n + (1./3)*Q2_n*Q2_n*Q2_n));

    for(double T = 0.0; T < 20000; T+=deltaT) {

      Q1_nMaisUm = q1_nMaisUm(Q1_n, deltaT, P1_n);
      Q2_nMaisUm = q2_nMaisUm(Q2_n, deltaT, P2_n);

      P1_nMaisUm = p1_nMaisUm(P1_n, deltaT, Q1_nMaisUm, Q2_nMaisUm);
      P2_nMaisUm = p2_nMaisUm(P2_n, deltaT, Q1_nMaisUm, Q2_nMaisUm);

      if(Q1_nMaisUm*Q1_n < 0 && P1_n > 0) {

        Q2s = q2s(P1_n, P2_n, Q2_n, -Q1_n);
        P2s = p2s(P1_n, P2_n, Q1_n, Q2_n, -Q1_n);

        dadosMapaDePoincare << Q2s << "\t" << P2s << endl;
      }

      Q1_n = Q1_nMaisUm;
      Q2_n = Q2_nMaisUm;
      P1_n = P1_nMaisUm;
      P2_n = P2_nMaisUm;
    }

  }

  return 0;
}
