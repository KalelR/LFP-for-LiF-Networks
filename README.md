# README for LFP directory.
This dir contains research about properties of LFP. Main focus is probabily on simulated networks.

# Programs
1. calcPSD: teste da biblioteca de PSD com seno
2. calcFFT teste do FFT
3. geraNoise: geraNoise com o OU Process do DiffEqs e calcula o PSD
4. geraNoise2: geraNoise com a implementaçao mais manual
5. LIF: codigo para a rede, tentando replicar a do paper
6. generatePoissonSpikeTrains: to generate ...; tests the result to make sure it is Poissonian. Eventually I should make better tests, as done in the Poisson Model of Spike Generation pdf, from which I got the algorithm 
7. generateConnections: gera a matriz de conexoes e os arquivos com os vetores I e E


# LIF
Cuidado com o dt na geraçao do ruido!
## Version
1. 
 - 0.1: implementei como diz no paper; Isyn ta dando 0 atualmente, precisa checar.
- 0.2: onde parei antes do RK4. Tava testando as condutancias
2.
    - 2.0: Runge Kutta RK4
    - 2.1: com metodo 2 de calculo do acoplamento (considerando todos os spikes passados tbm). Agora parece replicar o paper! Testei qtos spikes precisa considerar: de 30 para 50 nao muda tnato, mas menos que isso muda bastante.

# cpp
## Version
2. Rede desacoplada


## tests
1. Aumentei a condutancia devido as correntes externas e plotei o potencial com as correntes em 1, 10, 20 para os neuronios. Da pra ver bem a mudança na derivada do potencial.
2. Teste 1 sem aumento da condutancia
3. Resumo dos testes com conexoes, gerado pelo LIF_teste_conexoes.jl

# GenerateNoise
Gerei ruido para usar como no paper atraves do noiseProcess de Ornstein-Uhlenbeck do DiffEqs. 
Acertei as unidades (constantes em ms), freq em Hz.
A média dá prox de zero
O PSD não dá como o descrito no paper (sei lá!)

Por via das duvidas, tambem implementei de outa forma o processo e os resultados sao bem proximos, entao duvido que seja problema com o DiffEqs. Claro, eu posso ter feito merda na interpretaçao das constantes ou algo assim. 
Tbm posso ter feito merda no calculo do PSD.
O calculo do PSD foi feito com a biblioteca do julia. Eu testei com seno e deu certo.


# Database
http://crcns.org/data-sets/vc/pvc-3/about/?searchterm=local%20field%20potential
http://crcns.org/data-sets/ac/ac-2/?searchterm=local%20field%20potential
http://crcns.org/data-sets/methods/ieeg-1/about-ieeg-1