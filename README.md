# Distribuição de probablidade de Weibull

Ao executar cargo run, o programa apresentará:

A quantidade de corpos de prova ensaiados.

A tensão média observada experimentalmente.

O valor estimado de k, parâmetro de forma da Weibull.

O valor estimado de λ, parâmetro de escala em MPa.

A probabilidade acumulada de ruptura em cada tensão consultada.

As tensões associadas a 5% e 10% de probabilidade acumulada de ruptura.

A função fit() estima os parâmetros a partir dos dados fornecidos; portanto, os resultados serão diferentes dos  obtidos com os parâmetros ilustrativos anteriores.

##  Como interpretar os parâmetros

k elevado: geralmente indica menor dispersão relativa da resistência, no modelo Weibull.


k baixo: geralmente indica maior dispersão relativa.


λ: é a tensão na qual a probabilidade acumulada de ruptura é de aproximadamente 63,2%, pois 
F(λ)=1−exp(−1)

 Um cuidado importante: a tensão correspondente a 5% de ruptura estimada pelo modelo não é automaticamente um limite de projeto seguro. Para isso, podem ser necessários intervalos de confiança, tamanho amostral adequado e os critérios da norma técnica aplicável.

Limitação do ajuste atual: ele pressupõe que todos os corpos de prova romperam e que os resultados são independentes. Se houver ensaios interrompidos sem ruptura, será necessário adaptar a verossimilhança para tratar dados censurados.

Os parâmetros k e λ são ajustados automaticamente.