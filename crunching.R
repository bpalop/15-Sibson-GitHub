library("kSamples", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
setwd("~/Documents/git/15-Sibson")
datos <- read.csv("./rdos/rdos.csv",header=T,sep=" ")
pesos <- levels(datos$Peso)
## [1] "Cdistancia"  [2]"Cdistancia2" [3]"CmediaA"     [4]"CmediaA2"    [5]"CmediaG"     
## [6] "CmediaG2"    [7] "Csibson"    [8]"Edistancia"  [9]"Edistancia2" [10]"EmediaA"     
## [11]"EmediaA2"    [12]"EmediaG"    [13] "EmediaG2"   [14]"Esibson" 

vectoresE <- {}
nombresE <- {}
vectoresC <- {}
nombresC <- {}
for (i in pesos){
  lista <- datos[which(datos$Peso==i),3:ncol(datos)]
  lista <- as.vector(t(lista))
  lista <- lista[!is.na(lista)]
  if (substr(i,1,1)=='E') {
    vectoresE=cbind(vectoresE,lista)
    nombresE=c(nombresE,i)
  }
  else{
    vectoresC=cbind(vectoresC,lista)
    nombresC=c(nombresC,i)
  }
}
colnames(vectoresE) <- nombresE
colnames(vectoresC) <- nombresC
vectoresE=as.data.frame(vectoresE)
vectoresC=as.data.frame(vectoresC)

resultadosVectores={}
for (i in 1:length(nombresE)){
  lista={}
  for (j in 1:length(nombresC)){
    valor=ad.test(vectoresE[,i], vectoresC[,j],method="asymptotic",dist=FALSE)  
    lista=cbind(lista,valor$ad)
  }
  resultadosVectores=rbind(resultadosVectores,lista)
}
resultadosVectores <- as.data.frame(resultadosVectores)
#columnas etiquetamos con experimento+valores del test y quitamos espacios blancos
colnames(resultadosVectores) <- gsub(" ", "", paste(rep(nombresC,each=3), colnames(valor$ad)),fixed=TRUE)
#filas etiquetamos con experimento+ version1 y version2 del test y quitamos espacios blancos
rownames(resultadosVectores) <- gsub(" ", "", paste(rep(nombresE,each=2),c("ver1","ver2")),fixed=TRUE)

#En este punto, los colnames son:
#[1] "CdistanciaAD"              "CdistanciaT.AD"            "Cdistanciaasympt.P-value" 
#[4] "Cdistancia2AD"             "Cdistancia2T.AD"           "Cdistancia2asympt.P-value"
#[7] "CmediaAAD"                 "CmediaAT.AD"               "CmediaAasympt.P-value"    
#[10] "CmediaA2AD"                "CmediaA2T.AD"              "CmediaA2asympt.P-value"   
#[13] "CmediaGAD"                 "CmediaGT.AD"               "CmediaGasympt.P-value"    
#[16] "CmediaG2AD"                "CmediaG2T.AD"              "CmediaG2asympt.P-value"   
#[19] "CsibsonAD"                 "CsibsonT.AD"               "Csibsonasympt.P-value"    
# De esta manera, los p-valores del asymptotic, que es el valor que nos gusta, están en las columnas:
# resultadosVectores[,seq(3,length(resultadosVectores[1,]),by=3)]
# Estudiamos los resultados:
# 1) Buscamos en la columna Csibson los valores que destacan, e.d., las medidas que para el punto extremo
#    más se parecen a esta. 
#    Nos sale que son la distancia2 y la mediaA2
# 2) Ahora buscamos en la diagonal esas dos, para ver si la medida para el punto extremo se parece a la 
#    del central. El camino es medidaCentral se parece a medidaExtema que se parece a sibsonCentral.
# 3) En esas dos, la dist2 se queda en las que no funcionan bien, pero la mediaA2 se mantiene como candidata.

# Visualmente, se aprecia esta similitud:
par(mfrow=c(1,3))
hist(vectoresC$Csibson)
hist(vectoresE$EmediaA2)
hist(vectoresC$CmediaA2)

# Pasamos entonces a estudiar, en lugar del vector de pesos tomado como una lista plana, la media de pesos
# que se obtiene con unos y otros métodos, así como las desviaciones típicas.

mediasE <- {}
mediasC <- {}
desviacionesE <- {}
desviacionesC <- {}
for (i in pesos){
  lista <- datos[which(datos$Peso==i),3:ncol(datos)]
  desviaciones <- apply(lista,1,sd,na.rm=TRUE)
  lista <- rowMeans(lista,na.rm=TRUE)
  if (substr(i,1,1)=='E') {
    mediasE <- cbind(mediasE,lista)
    desviacionesE<-cbind(desviacionesE,desviaciones)
  }
  else{
    mediasC=cbind(mediasC,lista)
    desviacionesC<-cbind(desviacionesC,desviaciones)
  }
}
colnames(mediasE) <- nombresE
colnames(mediasC) <- nombresC
mediasE=as.data.frame(mediasE)
mediasC=as.data.frame(mediasC)

colnames(desviacionesE) <- nombresE
colnames(desviacionesC) <- nombresC
mediasE=as.data.frame(desviacionesE)
mediasC=as.data.frame(desviacionesC)


resultadosMedias <- {}
resultadosDesviaciones <- {}
for (i in 1:length(nombresE)){
  lista <- {}
  listaSD <- {}
  for (j in 1:length(nombresC)){
    valor <- ad.test(mediasE[,i], mediasC[,j],method="asymptotic",dist=FALSE)  
    valorSD <- ad.test(desviacionesE[,i], desviacionesE[,j],method="asymptotic",dist=FALSE)  
    lista <- cbind(lista,valor$ad)
    listaSD <- cbind(listaSD,valorSD$ad)
  }
  resultadosMedias=rbind(resultadosMedias,lista)
  resultadosDesviaciones=rbind(resultadosDesviaciones,listaSD)
}
resultadosMedias <- as.data.frame(resultadosMedias)
#columnas etiquetamos con experimento+valores del test y quitamos espacios blancos
colnames(resultadosMedias) <- gsub(" ", "", paste(rep(nombresC,each=3), colnames(valor$ad)),fixed=TRUE)
#filas etiquetamos con experimento+ version1 y version2 del test y quitamos espacios blancos
rownames(resultadosMedias) <- gsub(" ", "", paste(rep(nombresE,each=2),c("ver1","ver2")),fixed=TRUE)

resultadosDesviaciones <- as.data.frame(resultadosDesviaciones)
#columnas etiquetamos con experimento+valores del test y quitamos espacios blancos
colnames(resultadosDesviaciones) <- gsub(" ", "", paste(rep(nombresC,each=3), colnames(valor$ad)),fixed=TRUE)
#filas etiquetamos con experimento+ version1 y version2 del test y quitamos espacios blancos
rownames(resultadosDesviaciones) <- gsub(" ", "", paste(rep(nombresE,each=2),c("ver1","ver2")),fixed=TRUE)

# Al hacer el print de estas dos:
rm <- resultadosMedias[,seq(1,length(resultadosMedias[1,]),by=3)]
rd <- resultadosDesviaciones[,seq(1,length(resultadosDesviaciones[1,]),by=3)]
# Vemos que todo cuadra en las medias!!
# En las medias de pesos por punto, nos sale que a SibsonC se le parecen las que más EmediaA y EmediaA2 (aprox. 5)
# Entonces, miramos la diagonal y vemos que EmediaA2 vs CmediaA2 da 11, también bajito.

# Pero no tanto en las desviaciones:
# Tenemos una cosa curiosa: el Anderson-Darling test devuelve un p-valor de 1 para todas
# las parejas de la misma medida, es decir, reconoce que son la misma distribución
# y toda la diagonal da un p-valor de 1.
# El problema es que la pareja Csibson vs EmediaA2, que es la que mejor iba, ya no es un buen valor :(


