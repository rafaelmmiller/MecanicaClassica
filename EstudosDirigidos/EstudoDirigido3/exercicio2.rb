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

# Pêndulo Simples: Movimento de Rotação

# Variáveis Globais:
$m = 0.1
$l = 1
$g = 10
$omega0 = Math.sqrt($g/$l)
$deltaT = 0.001


def Hamiltoniana (phi, pPhi)
  return pPhi**2/(2 * $m * $l**2) - $m * $l**2 * $omega0**2 * Math.cos(phi)
end

def Phi (theta0, i0, t)
  return theta0 + (i0 * t)/($m * $l**2) + (($m**2 * $l**4 * $omega0**2)/i0**2)*Math.sin(theta0 + (i0 * t)/($m * $l**2))
end

def P_phi (theta0, i0, t)
  return i0 + (($m**2 * $l**4 * $omega0**2)/i0)*Math.cos(theta0 + (i0 * t)/($m * $l**2))
end

def Omega (i)
  return $omega0 - i/(8 * $m * $l**2)
end


# Para o exercício 1.13
umArquivo_numerico = File.open("numericoRotacao4.dat","w")
umArquivo_perturbado = File.open("perturbadoRotacao4.dat","w")

# Para o exercício 1.14
umArquivo_numerico_EspacoDeFase = File.open("numericoEspacoDeFaseRotacao4.dat","w")
umArquivo_pertubado_EspacoDeFase = File.open("perturbadoEspacoDeFaseRotacao4.dat","w")

$T = 0.01
i0 = 0.5
theta0 = Math::PI/10

phi = Phi(theta0, i0, $T)
pPhi = P_phi(theta0 , i0, $T)

puts Hamiltoniana(phi,pPhi)

while $T < 10.0 do
  phi_n = phi + $deltaT * pPhi/($m * $l**2 )
  pPhi_n = pPhi - $deltaT * $m * $l**2 * $omega0**2 * Math.sin(phi_n)

  umArquivo_numerico.syswrite("#{$T}\t#{phi_n}\n")
  umArquivo_perturbado.syswrite("#{$T}\t#{Phi(theta0,i0,$T)}\n")

  umArquivo_pertubado_EspacoDeFase.syswrite("#{Phi(theta0,i0,$T)}\t#{P_phi(theta0,i0,$T)}\n")
  umArquivo_numerico_EspacoDeFase.syswrite("#{phi_n}\t#{pPhi_n}\n")

  phi = phi_n
  pPhi = pPhi_n
  $T += $deltaT
end

=begin
#Para o exercicio 1.15
umArquivo_perturbado_Periodo = File.open("perturbadoPeriodo.dat","w")

k = -1.0
i0 = 0.001

while k < 1 do
  omega = Omega(i0)
  t = $omega0 / omega
  k = $omega0 * i0 - (i0**2)/(16.0*$m*$l**2) - $m* $l**2* $omega0**2
  i0 = i0 + 0.001
  umArquivo_perturbado_Periodo.syswrite("#{k}\t#{t}\n")
end
=end
