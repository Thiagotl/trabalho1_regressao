
# pacotes uteis
library(hnp) # pacote para envelope simulado
library(lmtest) # teste reset
library(car) # para teste de multicolinearidade (fatores de inflacao de variancia)
library(tseries) # teste de Jarque-Bera
## Variaveis:
# peso (y): peso da criancas ao nascer (gramas)
# pnatal (x1): variavel binaria que informa se a mae fez (1) ou nao (2) pre-natal
# pesom (x2): peso da mae (Kg)
# renda (x3): renda familiar percapta da mae (reais)
# A ideia eh modelar y em funcao das demais


## Importando os dados
dados<- read.table("dados-trabalho1.txt",h=T) 


## analise descritiva
summary(dados) # medidas descritivas
cor(dados) # correlacao
plot(dados) # diagrama de dispersao

## ajustando o modelo
fit <- lm(peso~pnatal+pesom+renda,data=dados) # ajustando o modelo
summary(fit)
# Verifique R2, R2-ajustado, teste F e testes t
# Mas note que soh podemos "acreditar" nesses testes apos analise de 
# diagnostico e os varios testes estudados. Veja a seguir.

fit_backup<-fit # para guardar esse primeiro fit

##### Analise de influencia
n<-dim(dados)[1] # tamanho da amostra

# com a seguinte funcao se obtem varias medidas de influencia
influence.measures(fit)

# Alavancagem
hatvalues(fit)
h_bar<-fit$rank / n
limite<-2*h_bar
abline(plot(hatvalues(fit),ylab="Alavancagem"), 
       col="red", h=limite,lty=2)
which(hatvalues(fit)>limite)

# DFFIT
dffits(fit)
limite<-2*sqrt(fit$rank / n)
abline(plot(dffits(fit),ylab="DFFITS"), 
       col="red", h=c(-limite,limite),lty=2)
which(abs(dffits(fit))>limite)

# DFBETA
dfbetas(fit) # cada beta tem seu DF

dfb1<-dfbetas(fit)[,1]
dfb2<-dfbetas(fit)[,2]
dfb3<-dfbetas(fit)[,3]
dfb4<-dfbetas(fit)[,4]

limite<-2/sqrt(n)
abline(plot(dfb1,ylab="DFBETA 1"), 
       col=c("red","blue","red"), h=c(-limite,0,limite),lty=c(2,1,2))

abline(plot(dfb2,ylab="DFBETA 2"), 
       col=c("red","blue","red"), h=c(-limite,0,limite),lty=c(2,1,2))

# distancia de Cook
cooks.distance(fit)
limite<-4/(n-fit$rank )
abline(plot(cooks.distance(fit),ylab="Distancia de Cook"), 
       col="red", h=limite,lty=2)


# residuo
residuo <- rstudent(fit) # residuo studentizado

plot(residuo,type='p',pch="+",main="Residuos",xlab="indices") # plota os residuos do modelo
abline(h=c(-2,0,2),lty=3) # inclui linhas horizontais no grafico

which(abs(residuo)>3)

hist(residuo) # histograma dos residuos

# envelope simulado baseado nos residuos studentizados
hnp(fit,resid.type="student",halfnormal = F) # envelope simulado 

## A observacao 45 eh bastante influente. Vamos investigar:
dados[45,]
summary(dados)
boxplot(peso~pnatal, data=dados)
# Ha algo muito estranho nesse individuo 45. Uma das rendas mais altas com o menor peso da crianca.
# Vamos retira-lo dos dados, pois e reestimar.


# Retirando observacao possivelmente influente e reestimando 

dados45<-dados[-45,]

## ajustando o modelo
fit <- lm(peso~pnatal+pesom+renda,data=dados45) # ajustando o modelo
summary(fit)
# Verifique R2, R2-ajustado, teste F e testes t
# Mas note que soh podemos "acreditar" nesses testes apos analise de 
# diagnostico e os varios testes estudados. Veja a seguir.

summary(fit_backup) 

# Comparando os dois modelos, perceba que muda muito. Ou seja, 
# a observacao 45 eh realmente influente.

##### Analise de influencia
n<-dim(dados45)[1] # tamanho da amostra

# com a seguinte funcao se obtem varias medidas de influencia
influence.measures(fit)

# Alavancagem
hatvalues(fit)
h_bar<-fit$rank / n
limite<-2*h_bar
abline(plot(hatvalues(fit),ylab="Alavancagem"), 
       col="red", h=limite,lty=2)

# DFFIT
dffits(fit)
limite<-2*sqrt(fit$rank / n)
abline(plot(dffits(fit),ylab="DFFITS"), 
       col="red", h=c(-limite,limite),lty=2)


# DFBETA
dfbetas(fit) # cada beta tem seu DF

dfb1<-dfbetas(fit)[,1]
dfb2<-dfbetas(fit)[,2]
dfb3<-dfbetas(fit)[,3]
dfb4<-dfbetas(fit)[,4]

limite<-2/sqrt(n)
abline(plot(dfb1,ylab="DFBETA 1"), 
       col=c("red","blue","red"), h=c(-limite,0,limite),lty=c(2,1,2))

abline(plot(dfb2,ylab="DFBETA 2"), 
       col=c("red","blue","red"), h=c(-limite,0,limite),lty=c(2,1,2))


# distancia de Cook
cooks.distance(fit)
limite<-4/(n-fit$rank )
abline(plot(cooks.distance(fit),ylab="Distancia de Cook"), 
       col="red", h=limite,lty=2)


# residuo
residuo <- rstudent(fit) # residuo studentizado

plot(residuo,type='p',pch="+",main="Residuos",xlab="indices") # plota os residuos do modelo
abline(h=c(-2,0,2),lty=3) # inclui linhas horizontais no grafico

hist(residuo) # histograma dos residuos

# envelope simulado baseado nos residuos studentizados
hnp(fit,resid.type="student",halfnormal = F) # envelope simulado 


# A observacao 43 apareceu "levemente" em algumas medidas. Mas eh possivel 
# verificar que sua exclusao nao influencia no ajuste. 
# Envelope estah muito bom e nao ha fortes evidencias de algo incorreto. Entao seguimos.




##### Testando as suposições do modelo

## [S0] O modelo estah corretamente especificado
## [S1] A media dos erros eh zero
## [s2] Homoscedasticidade dos erros
## [S3] Nao autocorrelacao 
## [S4] Ausencia de Multicolinearidade
## [S5] Normalidade dos erros

## Obs.: Para testes de hipoteses, se p-value < alpha (5%) 
## entao rejeita a hipotese nula (H0)

## Testa [S0]
## Teste RESET de especificacao
## H0: O modelo estah corretamente especificado
resettest(fit)

## Testa [S1]
## Teste t para a média dos errros
## H0: média dos erros eh igual a zero
t.test(resid(fit),mu=0,alternative="two.sided")

## Testa [s2]
## Teste de Bressch-Pagan (Koenker) de Heteroscedasticidade
## H0: erros sao homoscedasticos
bptest(fit, studentize = TRUE)

## Testa [S3]
## Teste de Durbin-Watson de autocorrelacao
## H0: : Nao hah autocorrelacao 
dwtest(fit)
acf(rstudent(fit))

## Testa [S4]
## Usa Fatores de Inflacao de Variancia para detectar multicolinearidade
## Regra de bolso: vif > 10 indica multicolinearidade. vif=1 seria o ideal.
vif(fit)

## Testa [S5]
## Teste Jarque-Bera de Normalidade
## H0: Os erros possuem distribuicao normal
jarque.bera.test(resid(fit))





## Agora com o modelo checado, com boas evidencias de que as suposicoes estao 
## satisfeitas, eh possivel fazer inferencias e interpretacoes. 
summary(fit)
# Quanto cada covariavel influencia na media de y?
# O quanto a mae ter feito pre-natal influencia no peso da crianca (em media)?

## Para fazer predicao, fazemos
predict(fit) # valores preditos na amostra usada na estimacao

# Predicao com novos dados
novos_dados <- data.frame(cbind(
  pnatal=1, pesom=70, renda=1000
))
predict(fit, newdata=novos_dados) # com novos valores (fora da amostra)

## Obs. Se seu banco de dados tiver muitas covariaveis candidatas, 
## voce pode utilizar a funcao step() para selecionar, via AIC, 
## um modelo candidato mais enxuto. 

step(fit) # seleciona modelo baseado no AIC e stepwise
# neste caso o modelo eh ok e o step nao excluiu covariaveis

# Obs.: dados simulados -> pode/deve não representar a realidade
