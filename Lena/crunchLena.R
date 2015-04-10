library("tiff", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
setwd("~/Documents/tests")
lena<-readTIFF("lena1.tiff")

pinta2X <- function(imagenRaster){
  #Pinta en pantalla la imagen Raster en tamaño doble de la resolucion original y con margenes
  #256x256 es la imagen
  width<-256
  height<-256
  margin<-10
  require(grDevices)
  op <- par(bg = "thistle")
  plot(c(0, 2*(width+margin)), c(0,2*(height+margin)), type = "n", xlab = "", ylab = "")
  rasterImage(imagenRaster,margin,margin,2*width+margin,2*height+margin,interpolate=FALSE)
}

pinta <- function(imagenRaster){
  #Pinta en pantalla la imagen Raster tal cual es la matriz
  #256x256 es la imagen
  width<-256
  height<-256
  #margin<-10
  require(grDevices)
  op <- par(bg = "thistle")
  plot(c(0, width), c(0,height), type = "n", xlab = "", ylab = "")
  rasterImage(imagenRaster,0,0,width,height,interpolate=FALSE)
}


pintapunto <- function(x,y,color){
  points(x,y,col=color,pch=19)
}

#Test we can print using mini-image
lenita<-as.raster(head(lena[,1:5]))
pinta(lenita) #Vemos que lenita tiene una matriz de 5x6 y se ve pixelada

#And now the real Lena
pinta(lena)

#Vamos a hacer un sampleo de lena para dejarla más sparse
sparsemat <- function(ratio,imagen){
  width <- 256
  height <- 256
  size <- ratio*width*height/100             #length of random number vectors
  set.seed(1) 
  x <- runif(size)*width          # generate samples from uniform distribution (0.0, 1.0)
  y <-runif(size)*height 
  df <-data.frame(x,y)
  copiaImagen <- matrix(0,width,height) #matriz de su mismo tamanyo
  vals <- c() #Vamos a guardar la columna de valores asignados al pixel
  #ratio es el valor porcentual de cuantos me quedo de cada muestra
  for (point in 1:length(df[,1])){
    valorEnMatriz <- imagen[as.integer(df[point,1]),as.integer(df[point,2])]
    copiaImagen[as.integer(df[point,1]),as.integer(df[point,2])] <- valorEnMatriz
    #vals<- c(vals,valorEnMatriz)
  }
  
  pinta(copiaImagen)
  #return(cbind(df,vals))
  return(copiaImagen)
}

isNeighbour <- function(vecs,row,point){
  numVec <- vecs[row,1]
  isNeigh <- FALSE
  k=2
  while (!isNeigh && k<=numVec+2){
    if (vecs[row,k]==point) isNeigh <- TRUE
    k <- k +1
  }
  return(isNeigh)
}

listNeighbours <- function(nube3D){
  tt <- delaunayn(nube3D[,1:2])
  vecinos <- matrix(-1,nrow = nrow(nube3D), ncol = 30)
  for (i in 1:nrow(vecinos)){
    vecinos[i,1]=0
  }
  for (i in 1:nrow(tt)){ 
    tresPuntos <- tt[i,]
    p1=tresPuntos[1]
    p2=tresPuntos[2]
    p3=tresPuntos[3]
    for (j in seq(1:3)){
      numVec <- vecinos[p1,1]
      if (!isNeighbour(vecinos,p1,p2)){
        vecinos[p1,numVec+2]=p2
        numVec <- numVec +1
      }
      if (!isNeighbour(vecinos,p1,p3)){
        vecinos[p1,numVec+2]=p3
        numVec <- numVec +1
      }
      vecinos[p1,1]=numVec
      aux<-p1
      p1<-p2
      p2<-p3
      p3<-aux
    }
  }
  return(vecinos)
}

getPoints3D <- function(imagen){
  puntos<-c()
  for (i in 1:length(imagen[1,])){
    for (j in 1:length(imagen[,1])){
      if(imagen[i,j]!=0){
        puntos<-rbind(puntos,c(i,j,imagen[i,j]))
      }
    }
  }
  return(puntos)
}


lenaBaja90<-sparsemat(90,lena)
pinta(lenaBaja90)
puntos3D90<-getPoints3D(lenaBaja90)

lenaBaja10<-sparsemat(10,lena)
pinta(lenaBaja10)
puntos3D10<-getPoints3D(lenaBaja10)
#adyacencias<-listNeighbours(puntos3D10) //Esto no hace falta porque el Voronoi lo hacemos en CGAL

write.csv(puntos3D10,file="lena10.csv")
write.csv(puntos3D90,file="lena90.csv")



pintapunto(puntos3D10[1,1],puntos3D10[1,2],puntos3D10[1,3])

