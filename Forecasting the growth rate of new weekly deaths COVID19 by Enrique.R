library(tidyverse)
library(lubridate)
library(quantmod)
library(xts)
library(tseries)
library(zoo)
library(reshape)
library(readr)
library(dplyr)
library(forecast)
library(seasonal)
library(tsibble)
library(urca)
library(uroot)
library(gridExtra)
library(Metrics)
library(scales)



########################### II.	METODOLOGÍA.#################################

#                         A. Base de datos: 
#Descargando la data sobre covid mundial
data <- read_csv("https://covid.ourworldindata.org/data/owid-covid-data.csv")

#                         B. Preparación de los datos:

#Filtrando los datos de Mexico
data <- data[data$location=="Mexico",]
#view(data)

#Filtrando datos de interes:
covid_mx <- select(data, date, new_deaths)
covid_mx <- covid_mx[1:880,] ## Limita la fecha hasta el día 29/may/22

################################################################################ 
# NOTA IMPORTATE:                                                              #
# Es importante aclarar que dicho código permite actualizar los datos          #
# históricos de manera automática a la fecha de consulta, en este sentido      #
# permite reajustar el modelo ARMA(1,4)-ARCH(1), y generar nuevos pronósticos  #
# correspondientes a las próximas cinco semanas a partir de la fecha de        #
# consulta, para lograr esto FAVOR DE COMENTAT LA LINEA 35.                    #
# Los pronosticos se imprimiran en consola al ejecutar la linea 179.           # 
################################################################################


#Limpiando valores nulos
print("Número de valores nulos:")
print(sum(is.na(covid_mx$new_deaths)))

i=0;
for (i in 1:nrow(covid_mx)){
  for (name in colnames(covid_mx)){
    if (i==1 & name!="date"){
      covid_mx[1,name]=0;
    }else{
      if(is.na(covid_mx[i,name])==TRUE  & is.na(covid_mx[i+1,name]==TRUE) ){
        covid_mx[i,name]=covid_mx[i-1,name]
      }
      if(is.na(covid_mx[i,name]==TRUE)  & is.na(covid_mx[i+1,name])==FALSE ){
        covid_mx[i,name]= round(((covid_mx[i-1,name])+(covid_mx[i+1,name]))/2)
      }
      
    }
    
  }
}


# Cambiando tiempo de muestreo de diario a semanal 
covid_mx <- covid_mx %>%
  mutate(week=tsibble::yearweek(covid_mx$date))

covid_week <- covid_mx %>%
  group_by(week) %>%
  summarise(new_deaths=round(mean(new_deaths)))

covid_week <- covid_week %>%
  mutate(date=seq(as.Date(covid_mx$date[1]),as.Date(covid_mx$date[length(covid_mx$date)]),by="week"))

covid_mx<-covid_week

#Calculo de la tasa de crecimiento
covid_mx <- covid_mx %>%
  mutate(covid_mx, TC=c(NaN,diff(log(covid_mx$new_deaths))))

#limpieza de valores nulos e Inf en la Tasa de crecimiento
covid_mx$TC<- ifelse(covid_mx$TC==Inf | covid_mx$TC==-Inf,NaN,covid_mx$TC)
covid_mx<-na.omit(covid_mx) 


# Graficando tasa de crecimiento nuevos fallecimientos:
ggplot(data=covid_mx, aes(x=date  , y=TC))+
  geom_line(colour="black")+
  scale_x_date(date_breaks = "6 month",
               date_labels = " %m-%Y")+
  scale_y_continuous(breaks = seq(-2,2,0.5))+
  theme_minimal()+
  labs(title="Tasa de crecimiento nuevos fallecimentos semanales por COVID-19",
       subtitle = "en Mexico del años 2020 al 2022",
       caption = "Elaboración EHL con con datos de Our World in Data COVID-19",
       x="Fecha",
       y="Tasa de crecimiento")
#                                   C.Modelado

#                                   C.1) Estacionariedad

# ADF 
adf.test(covid_mx$TC)

# Phillips-Perron 
pp.test(covid_mx$TC)

# KPSS 
kpss.test(covid_mx$TC, null='Level')


#                     C.2) Identificación de la volatilidad

covid_mx$TC_sq <- (covid_mx$TC)^2

# Graficando tasa de crecimiento al cuadrado de nuevos fallecimientos:
ggplot(data=covid_mx, aes(x=date  , y=TC_sq))+
  geom_line(colour="black")+
  scale_x_date(date_breaks = "6 month",
               date_labels = " %m-%Y")+
  theme_minimal()+
  labs(title="Tasa de crecimiento al cuadrado de nuevos fallecimentos diarios \npor COVID-19",
       subtitle = "en Mexico del años 2020 al 2022",
       caption = "Elaboración propia con con datos de Our World in Data COVID-19",
       x="Fecha",
       y="Tasa de crecimiento al cuadrado")


## Prueba de autocorrelación de los retornos al cuadrado

Box.test(covid_mx$TC_sq, lag=round(log(length(covid_mx$TC_sq))), type = 'Ljung-Box')

#Correlación serial

## Relacionamos el ACF con la prueba anterior 
library(forecast)
Acf(covid_mx$TC_sq, lag.max = 20, plot = F) %>% 
  ggplot2::autoplot()+
  labs(title = 'Función de autocorrelación simple',
       subtitle = 'Serie: Tasa de crecimiento al cuadrado de nuevos fallecimentos diarios por COVID-19')

# Prueba de ARCH a los retornos financieros 

library(FinTS)

ArchTest(covid_mx$TC, lags = 1, demean = TRUE)


#             C.3) Estimación del modelo ARMA-ARCH
library(fGarch)
library(rugarch)
library(forecast)
# Covertir el DataFrame en ts con frecuencia diaria
covid_ts<-ts(data=covid_mx$TC,
                 start = c(2020,1),
                 frequency = 52)
# ARMA-GARCH 
#Coeficientes ARMA
auto.arima(covid_ts, max.d = 0)

#Modelo hibrido
armagarch<-garchFit(~arma(1,4)+garch(1,0), trace = FALSE, data=covid_ts)
armagarch

#                             D. Evaluación del modelo:
summary(armagarch)

plot(armagarch,which=3)
fcast_armagarch <- predict(armagarch, n.ahead=5, plot=T)
                                    

                                    #E. Pronostico:
fcast_armagarch

#Graficar pronostico 1er semana de junio a la 1er semana de julio

## Pronostico tasa de crecimineto
Estimado_mean <- c(1:6)
Estimado_mean[1] <- covid_ts[113]
Estimado_mean[2:6] <- (fcast_armagarch$meanForecast )
Estimado_mean=Estimado_mean <-ts(data=Estimado_mean,
                                   start = c(2022,9),frequency = 52)
# Pronostico volatilidad
volatilidad <- ts(data = fcast_armagarch$meanError, start = c(2022,10),
                  frequency = 52)   

# Pronostico intervalo de confianza
Estimado_up<- (fcast_armagarch$upperInterval)
Estimado_down<- (fcast_armagarch$lowerInterval)

Estimado_up=Estimado_up <-ts(data=Estimado_up,
                             start = c(2022,10),
                             frequency = 52)
Estimado_down=Estimado_down <-ts(data=Estimado_down,
                                 start = c(2022,10),
                                 frequency = 52)
#Actualizar datos historios a 1er semana de junio a la 1er semana de julio
covid_mx <- select(data, date, new_deaths)
covid_mx <- covid_mx[1:912,] ## Limita la fecha a la 1er semana de julio 2022

i=0;
for (i in 1:nrow(covid_mx)){
  for (name in colnames(covid_mx)){
    if (i==1 & name!="date"){
      covid_mx[1,name]=0;
    }else{
      if(is.na(covid_mx[i,name])==TRUE  & is.na(covid_mx[i+1,name]==TRUE) ){
        covid_mx[i,name]=covid_mx[i-1,name]
      }
      if(is.na(covid_mx[i,name]==TRUE)  & is.na(covid_mx[i+1,name])==FALSE ){
        covid_mx[i,name]= round(((covid_mx[i-1,name])+(covid_mx[i+1,name]))/2)
      }
      
    }
    
  }
}

######## Cambiando tiempo de muestreo de diario a semanal ################

covid_mx <- covid_mx %>%
  mutate(week=tsibble::yearweek(covid_mx$date))

covid_week <- covid_mx %>%
  group_by(week) %>%
  summarise(new_deaths=round(mean(new_deaths)))

covid_week <- covid_week %>%
  mutate(date=seq(as.Date(covid_mx$date[1]),as.Date(covid_mx$date[length(covid_mx$date)]),by="week"))

covid_mx<-covid_week

#Calculo de la tasa de crecimiento
covid_mx <- covid_mx %>%
  mutate(covid_mx, TC=c(NaN,diff(log(covid_mx$new_deaths))))

#limpieza de valores nulos e Inf en la Tasa de crecimiento
covid_mx$TC<- ifelse(covid_mx$TC==Inf | covid_mx$TC==-Inf,NaN,covid_mx$TC)
covid_mx<-na.omit(covid_mx) 

covid_ts<-ts(data=covid_mx$TC,
             start = c(2020,1),
             frequency = 52)

## Graficar pronostico vs historicos

ggplot()+
  geom_line(data=covid_ts, aes(x=time(covid_ts), y=covid_ts, 
                                   colour='TC historica'),size=0.8)+
  geom_line(data=volatilidad,aes(x=time(volatilidad),y=volatilidad,colour='Pronóstico volatilidad \nARMA(1,4)-ARCH(1)\n '),size=0.8)+
  geom_line(data=Estimado_mean,aes(x=time(Estimado_mean),y=Estimado_mean,
                                    colour='Pronóstico puntual TC\nARMA(1,4)-ARCH(1)\n '),size=0.8)+
  geom_ribbon(aes(x=time(Estimado_down), ymin=Estimado_down, 
                  ymax=Estimado_up,fill = "Intervalo de confianza al\n95% pronóstico puntual TC"), alpha = 0.3)+
  scale_y_continuous(breaks = seq(-2,2,(0.5)))+
  scale_colour_manual(values=c('red','green','steelblue'))+
  scale_fill_manual("",values="grey12")+
  labs(title = "Pronóstico tasa de crecimiento fallecimientos semanales por COVID-19 en México",
       caption = "Elaboración EHL con datos de Our World in Data COVID-\n TC=tasa de crecimiento fallecimientos semanales por COVID-19 en México",
       x='Fecha',
       y='Tasa de crecimiento',
       colour='Series')
