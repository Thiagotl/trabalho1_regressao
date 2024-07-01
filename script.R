library(hnp) # pacote para envelope simulado
library(lmtest) # teste reset
library(car) # para teste de multicolinearidade (fatores de inflacao de variancia)
library(tseries) # teste de Jarque-Bera
library(gtsummary)
library(dplyr)



#data set
dados<- read.table("dados-trabalho1.txt",h=T) 
View(dados)
summary(dados)
dados |> 
  select(y,x2,x3,x4,x5,x6,x7) |> 
  tbl_summary(type = all_continuous() ~ "continuous2",
              statistic = all_continuous()~c("{mean}","{median}","{min}", "{max}"))
summary(dados)
cor(dados) # ha multicolinearidade entre x3 e x4
plot(dados)


## ajustando o modelo
fit <- lm(y~x1+x2+x3+x4+x5+x6+x7,data=dados) # ajustando o modelo
#summary(fit)

#tbl_regression(fit)
#teste F: Pelo teste F se rejeitou HO, ou seja, pelo menos um beta é diferente de zero.
#teste t: para o modelo em estudo x1 e x6 foram significativos para o modelo, visto que o p-valor foi menor que 0.05

##TESTE DAS SUPOSIÇÕES DO MODELO

## Testa [S0]
## Teste RESET de especificacao
## H0: O modelo esta corretamente especificado

resettest(fit) # como p-valor > alfa, então não se rejeita H0 e o modelo esta corretamente especifcado

## Testa [S1]
## Teste t para a média dos erros
## H0: média dos erros eh igual a zero
t.test(resid(fit),mu=0,alternative="two.sided") # como o p-valor foi de 1 > que alfa, entao nao se rejeita H0
                                                # média dos erros é igual a 0 

## Testa [s2]
## Teste de Bressch-Pagan (Koenker) de Heteroscedasticidade
## H0: erros sao homoscedasticos
bptest(fit, studentize = TRUE) #como o p-valor foi de 0.07121 > 0.05, entao nao se rejeita H0, logo o erros são homoscedasticos


## Testa [S3]
## Teste de Durbin-Watson de autocorrelacao
## H0: : Nao hah autocorrelacao 
dwtest(fit)
acf(rstudent(fit))

## Testa [S4]
## Usa Fatores de Inflacao de Variancia para detectar multicolinearidade
## Regra de bolso: vif > 10 indica multicolinearidade. vif=1 seria o ideal.
vif(fit) #Através dos vifs podemos obsevar que x1, x2, x6 e x7 temos o ideal

## Testa [S5]
## Teste Jarque-Bera de Normalidade
## H0: Os erros possuem distribuicao normal
jarque.bera.test(resid(fit)) # como p-valor > que 0.05, não rejeitamos H0


## Agora com o modelo checado, com boas evidencias de que as suposicoes estao 
## satisfeitas, eh possivel fazer inferencias e interpretacoes. 
summary(fit)


## Para fazer predicao, fazemos
predict(fit) # valores preditos na amostra usada na estimacao

# Predicao com novos dados
novos_dados <- data.frame(cbind( ))
predict(fit, newdata=novos_dados) # com novos valores (fora da amostra)



##### ANALISE DE INFLUENCIA
n<-dim(dados)[1] # tamanho da amostra

# com a seguinte funcao se obtem varias medidas de influencia
influence.measures(fit)  # LINHAS COM POSSIVEIS INFLUENCIAS 18,23,64,63

# Alavancagem
hatvalues(fit)
h_bar<-fit$rank / n
limite<-2*h_bar
abline(plot(hatvalues(fit),ylab="Alavancagem"), 
       col="red", h=limite,lty=2)
which(hatvalues(fit)>limite) # potencialmente influentes e podem merecer uma 
                             #análise mais detalhada para verificar se são outliers ou se têm um impacto significativo no modelo ajustado.

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
dfb5<-dfbetas(fit)[,5]
dfb6<-dfbetas(fit)[,6]
dfb7<-dfbetas(fit)[,7]

limite<-2/sqrt(n)
abline(plot(dfb1,ylab="DFBETA 1"), 
       col=c("red","blue","red"), h=c(-limite,0,limite),lty=c(2,1,2))

abline(plot(dfb2,ylab="DFBETA 2"), 
       col=c("red","blue","red"), h=c(-limite,0,limite),lty=c(2,1,2))

abline(plot(dfb3,ylab="DFBETA 3"), 
       col=c("red","blue","red"), h=c(-limite,0,limite),lty=c(2,1,2))

abline(plot(dfb4,ylab="DFBETA 4"), 
       col=c("red","blue","red"), h=c(-limite,0,limite),lty=c(2,1,2))

abline(plot(dfb5,ylab="DFBETA 5"), 
       col=c("red","blue","red"), h=c(-limite,0,limite),lty=c(2,1,2))

abline(plot(dfb6,ylab="DFBETA 6"), 
       col=c("red","blue","red"), h=c(-limite,0,limite),lty=c(2,1,2))

abline(plot(dfb7,ylab="DFBETA 7"), 
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
                                            #Acima no gráfico do envelope dos resíduos, 
                                            #temos eles "studentizados" e ordenados. 
                                            #Esse envelope são bandas de confiança e esperamos que pelo menos 90% dos resíduos estejam entre essas bandas. Essas bandas são feitas através de simulações de Monte Carlo. 
                                            #Esse envelope indica, principalmente, que a suposição distribucional está correta. Ou seja, a distribuição normal é correta para modelar esses dados.



# Pelo VIF , temos que remover x3 e x4, porém temos x4 possui uma correlaçao maior com y do que x3 coom y
# e tambem temos que pela distancia de cook a observacao 23 e a mais influente 





dados23<-dados[-23,]
#head(dados23, 25)
# remocao de x3 divido a alta a correlacao alta com x4
fit2<-lm(y~x1+x2+x4+x5+x6+x7,data=dados)

summary(fit2)

#tenho que x1, x4 e x6 foram significativos
n<-dim(dados23)[1] # tamanho da amostra



## Obs. Se seu banco de dados tiver muitas covariaveis candidatas, 
## voce pode utilizar a funcao step() para selecionar, via AIC, 
## um modelo candidato mais enxuto. 

step(fit) # seleciona modelo baseado no AIC e stepwise
# neste caso o modelo eh ok e o step nao excluiu covariaveis

# Obs.: dados simulados -> pode/deve não representar a realidade





