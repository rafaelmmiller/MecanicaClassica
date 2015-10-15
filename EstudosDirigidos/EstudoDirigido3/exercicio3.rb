#########################################################################
# PGF 5005 - Mecânica Clássica
# Professor Iberê L. Caldas
# Alunos: Elion Hack       NUSP.:
#         Rafael M. Miller NUSP.: 7581818
#
#   Terceiro Estudo Dirigido
# http://web.if.usp.br/controle/sites/web.if.usp.br.controle/
#                                   files/EstudoDirigido3.pdf
#########################################################################

# Variáveis Globais:
$deltaT = 0.001

def Hamiltoniana (q1, q2, p1, p2, f)
  return (1/2)*( p1**2 + p2**2 ) + f * q2 + (1/2) * (1 - f - Math.sqrt(q1**2 + (1 - q2)**2))**2
end

umArquivo_mapaDePoincare = File.open("mapaDePoincare.dat","w")

q2 = 0.0
p1 = 0.0
q1 = 0.05

f = 0.2
E = 0.03
p2 = Math.sqrt(2*E - 2 * f * q2 - (1 - f - Math.sqrt(q1**2 + (1 - q2)**2))**2 - p1**2)

$T = 0.01

while $T < 10.0 do
  q1_n = q1 + $deltaT * p1
  q2_n = q2 + $deltaT * p2

  p1_n = p1 + $deltaT * ((q1 * (1 - f - Math.sqrt(q1**2 + (1 - q2)**2)))/Math.sqrt(q1**2 + (1 - q2)**2))
  p2_n = p2 - $deltaT * (f + ((1 - f - Math.sqrt(q1**2 + (1 - q2)**2)) * (1 - q2))/Math.sqrt(q1**2 + (1 - q2)**2))

  if (q1_n * q1 < 0) and (p1 > 0)

    q2s = q2_n - q2_n*p2_n/p1_n
    p2s = p2_n + p2_n*(f + ((-1 + f + Math.sqrt(q1**2 + (-1 + q2)**2)) * (-1 + q2))/Math.sqrt(q1**2 + (-1 + q2)**2))/p1

    umArquivo_mapaDePoincare.syswrite("#{q2s}\t#{p2s}\n")
  end

  q1 = q1_n
  q2 = q2_n
  p1 = p1_n
  p2 = p2_n
  $T += $deltaT

end

# Perturbado
=begin
while $T < 10.0 do
  psi_n = psi + $deltaT * ($omega0 - (j/(3 * $m * $l * $l))) * Math.sin(psi) ** 4
  j_n = j + $deltaT * ((2 * j * j)/(3 * $m * $l * $l)) * (Math.sin(psi_n) ** 3) * Math.cos(psi_n)

  umArquivo_perturbado.syswrite("#{$T}\t#{psi_n}\n")
  umArquivo_pertubado_EspacoDeFase.syswrite("#{Phi(j_n,psi_n)}\t#{P_phi(j_n,psi_n)}\n")

  psi = psi_n
  j = j_n
  $T += $deltaT
end


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
=end
