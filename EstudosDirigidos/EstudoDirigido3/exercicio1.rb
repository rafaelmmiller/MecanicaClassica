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

# Pêndulo Simples: Movimento de Libração

# Variáveis Globais:
$m = 0.1
$l = 1
$g = 10
$omega0 = Math.sqrt($g/$l)
$deltaT = 0.001


def H0 (j)
  return j * $omega0
end

def H1 (j, psi)
  return -((j * j)/(6 * $m * $l * $l)) * Math.sin(psi) ** 4
end

def Hamiltoniana (j, psi)
  return H0(j) + H1(j,psi) - $m * $l * $l * $omega0 * $omega0
end

def Phi (j, psi)
  return Math.sqrt((2*j)/($m * $l**2 * $omega0)) * Math.sin(psi)
end

def P_phi (j, psi)
  return Math.sqrt(2 * j * $m * $l**2 * $omega0) * Math.cos(psi)
end

def Omega (i)
  return $omega0 - i/(8 * $m * $l**2)
end

def J_perturbacao (theta0, i0, t)
  return i0 + (i0**2)/(12 * $m * ($l**2) * $omega0) * ((1/4) * Math.cos(4 * Omega(i0) * t + 4 * theta0) - Math.cos(2 * Omega(i0) * t + 2 * theta0))
end

def Psi_perturbacao (theta0, i0, t)
  return Omega(i0) * t + theta0 + (i0/(12 * $m * ($l**2) * $omega0)) * (Math.sin(2 * Omega(i0) * t + 2 * theta0) - ((1/8) * Math.sin(4 * Omega(i0) * t + 4 * $omega0)))
end

# Para o exercício 1.13
umArquivo_numerico = File.open("numericoLibracao3.dat","w")
umArquivo_perturbado = File.open("perturbadoLibracao3.dat","w")

# Para o exercício 1.14
umArquivo_numerico_EspacoDeFase = File.open("numericoEspacoDeFaseLibracao3.dat","w")
umArquivo_pertubado_EspacoDeFase = File.open("perturbadoEspacoDeFaseLibracao3.dat","w")

$T = 0.01
i0 = 0.8
theta0 = Math::PI/2

psi = Psi_perturbacao(theta0, i0, $T)
j = J_perturbacao(theta0,i0,$T)

phi = Phi(j,psi)
pPhi = P_phi(j,psi)

puts Hamiltoniana(j,psi)

while $T < 10.0 do
  phi_n = phi + $deltaT * pPhi/($m * $l**2 )
  pPhi_n = pPhi - $deltaT * $m * $l**2 * $omega0**2 * Math.sin(phi_n)

  umArquivo_numerico.syswrite("#{$T}\t#{phi_n}\n")
  umArquivo_perturbado.syswrite("#{$T}\t#{Phi(J_perturbacao(theta0,i0,$T),Psi_perturbacao(theta0,i0,$T))}\n")

  umArquivo_pertubado_EspacoDeFase.syswrite("#{Phi(J_perturbacao(theta0,i0,$T),Psi_perturbacao(theta0,i0,$T))}\t#{P_phi(J_perturbacao(theta0,i0,$T),Psi_perturbacao(theta0,i0,$T))}\n")
  umArquivo_numerico_EspacoDeFase.syswrite("#{phi_n}\t#{pPhi_n}\n")

  phi = phi_n
  pPhi = pPhi_n
  $T += $deltaT
end

#Para o exercicio 1.15
umArquivo_perturbado_Periodo = File.open("perturbadoPeriodo3.dat","w")

k = -1.0
i0 = 0.001

while k < 1 do
  omega = Omega(i0)
  t = $omega0 / omega
  k = $omega0 * i0 - (i0**2)/(16.0*$m*$l**2) - $m* $l**2* $omega0**2
  i0 = i0 + 0.001
  umArquivo_perturbado_Periodo.syswrite("#{k}\t#{t}\n")
end
