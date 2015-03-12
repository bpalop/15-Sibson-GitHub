// bpalop, 2015.03.12

// A partir de los conjuntos de datos almacenados en ./data/conj* calculamos Voronois
// Tomamos un punto "central" en el Voronoi, el (0,0) y analizamos sus vecinos naturales
// Hacemos igual para un punto "extremo"

// Para cada medida, generamos un fichero que tiene una linea por experimento
// con los pesos de cada vecino natural segun la formula propuesta

// ./rdos/rdosC* indica que el fichero es del punto central
// ./rdos/rdosE* indica que el fichero es del punto extremo

// *sibson pondera peso de cada vecino segun área que le roba = sibN (valor normalizado 0-1)
// *distancia pondera con 1/distancia = idistN (1o calcula y normaliza 0-1)
// *distancia2 pondera con 1/distancia^2 = idist2N (1o calcula y normaliza 0-1)
// *mediaA pondera con sibN+idistN/2 (media de los valores ya normalizados y se re-normaliza)
// *mediaG pondera con sqrt(sibN*idistN) (partiendo de valores ya normalizados y se re-normaliza)

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <vector>
#include <fstream>
#include <iterator>
#include <string.h>

//Numero de conjuntos que usamos
#define NF 200
//Numero de ficheros de salida
#define NFS 20
//Numero de puntos que tomamos de cada fichero (para rodear bien al centro)
//#define NP 100
#define NP 50

//Paso con el que vamos metiendo el punto extremo hasta que entra en el 
//cierre convexo
#define INCREMENTO 0.0001

using namespace CGAL;

typedef Simple_cartesian<float>             K;
typedef K::Point_2                          Point2;
typedef K::FT                               flotante;
typedef std::vector<Point2>                 Vector2;
typedef std::istream_iterator<Point2>       istream_iterator;
typedef CGAL::Delaunay_triangulation_2<K>   Delaunay;
typedef std::vector< std::pair< Point2, K::FT > > Point_coordinate_vector;

Delaunay CargarPuntos(int fich){

//Get the pointset in file ./data/conjXXX, where XXX=fich
//Return a Delaunay triangulation with that pointset

  char nombreFich[100]="./data/conj";
  int longNombre=strlen(nombreFich);
  //Generate filename
  nombreFich[longNombre]='0'+fich/100;
  nombreFich[longNombre+1]='0'+(fich%100)/10;
  nombreFich[longNombre+2]='0'+fich%10;
  nombreFich[longNombre+3]='\0';
  //std::cout<<nombreFich<<std::endl;

  // Prepare a vector for NP points.
  Vector2 puntos;
  puntos.reserve(NP);
  
  std::ifstream fichero;
  fichero.open(nombreFich);
  if(!fichero.is_open()){
    std::cout<<"Error de fichero de entrada\n";
    exit(-1);
  }

  // Get the points and fill the vector
  std::copy(istream_iterator(fichero), istream_iterator(), std::back_inserter(puntos));


  fichero.close();

  // All points are now in the vector puntos
  Delaunay dt;
  //for(int i=0;i<puntos.size();i++){
  for(int i=0;i<NP;i++){
    dt.insert(puntos[i]);
  }

  return(dt);
}

Point2 BuscarExtremo(Delaunay dt){
  //Coordinates computations for the extreme point so
  //that it falls inside the convex hull
  Point2 extremo(100.0, 0);
  int veces=1;
  Point_coordinate_vector coords;
  CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>,K::FT, bool> 
    totalArea = CGAL::natural_neighbor_coordinates_2(dt, extremo, std::back_inserter(coords));
  do{
    totalArea = CGAL::natural_neighbor_coordinates_2(dt, extremo, std::back_inserter(coords));
    if(!totalArea.third){
      extremo=Point2(100.0-veces*INCREMENTO,0);
      veces++;
    }
  }while(!totalArea.third);

  return(extremo);
}

int Menu(){
  int opcion;
  do{
    printf("\n Indica las salidas a generar para los puntos central y extremo\n");
    printf("\n Sale una linea por punto\n");
    printf("1.- Peso sibson \t\t\t\t (./rdos/rdos[C-E]sibson.csv) \n");
    printf("2.- Inverso de distancia \t\t\t\t (./rdos/rdos[C-E]distancia.csv) \n");
    printf("3.- Inverso de distancia cuadrado \t\t\t\t (./rdos/rdos[C-E]distancia2.csv) \n");
    printf("4.- Media aritmetica sibson, inverso de distancia \t\t\t\t (./rdos/rdos[C-E]mediaA.csv) \n");
    printf("5.- Media aritmetica sibson, inverso de distancia cuadrado \t\t\t\t (./rdos/rdos[C-E]mediaA2.csv) \n");
    printf("6.- Media geometrica sibson, inverso de distancia \t\t\t\t (./rdos/rdos[C-E]mediaG.csv) \n");
    printf("7.- Media geometrica sibson, inverso de distancia cuadrado \t\t\t\t (./rdos/rdos[C-E]mediaG2.csv) \n");
    printf("8.- Generar todos los ficheros\n");

/* Version anterior. Guardo por ahora
    printf("2.- Peso 1/distancia \t\t\t\t (./rdos/ppDistancia1[C-E].csv) \n");
    printf("3.- Peso 1/distancia^2 \t\t\t\t (./rdos/ppDistancia2[C-E].csv) \n");
    printf("4.- M. Aritmetica pesos (Psib+Pd^-1) \t\t (./rdos/ppAritmetica1[C-E].csv) \n");
    printf("5.- M. Aritmetica pesos (Psib+Pd^-2) \t\t (./rdos/ppAritmetica2[C-E].csv) \n");
    printf("6.- M. Geometrica pesos (Psib*Pd^-1) \t\t (./rdos/ppGeometrica1[C-E].csv) \n");
    printf("7.- M. Geometrica pesos (Psib*Pd^-2) \t\t (./rdos/ppGeometrica2[C-E].csv) \n");
    printf("8.- Generar todos los ficheros\n");
*/

    printf("\n0.- Salir\n\n");
    printf("Opcion:  ");
    scanf("%d",&opcion);
  }while((opcion<0)||(opcion>9));

  return(opcion);
}

void AbrirFicheros(int opcion,std::ofstream *ficheros){


  if(opcion==1 || opcion==8){
    ficheros[0].open("./rdos/rdosCsibson.csv");
    ficheros[1].open("./rdos/rdosEsibson.csv");
    if((!ficheros[0].is_open())||(!ficheros[1].is_open())){
      printf("Error de creacion de ficheros!!!\n");
      exit(-1);
    }
    for (int i=1;i<=9;i++){  ficheros[0]<<"N"<<i<<" "; }
    ficheros[0]<<std::endl;
    for (int i=1;i<=14;i++){ ficheros[1]<<"N"<<i<<" "; }
    ficheros[1]<<std::endl;
  }
  if(opcion==2 || opcion==8){
    ficheros[2].open("./rdos/rdosCdistancia.csv");
    ficheros[3].open("./rdos/rdosEdistancia.csv");
    if((!ficheros[2].is_open())||(!ficheros[3].is_open())){
      printf("Error de creacion de ficheros!!!\n");
      exit(-1);
    }
    for (int i=1;i<=9;i++){  ficheros[2]<<"N"<<i<<" "; }
    ficheros[2]<<std::endl;
    for (int i=1;i<=14;i++){ ficheros[3]<<"N"<<i<<" "; }
    ficheros[3]<<std::endl;
  }
  if(opcion==3 || opcion==8){
    ficheros[4].open("./rdos/rdosCdistancia2.csv");
    ficheros[5].open("./rdos/rdosEdistancia2.csv");
    if((!ficheros[4].is_open())||(!ficheros[5].is_open())){
      printf("Error de creacion de ficheros!!!\n");
      exit(-1);
    }
    for (int i=1;i<=9;i++){  ficheros[4]<<"N"<<i<<" "; }
    ficheros[4]<<std::endl;
    for (int i=1;i<=14;i++){ ficheros[5]<<"N"<<i<<" "; }
    ficheros[5]<<std::endl;
  }

  if(opcion==4 || opcion==8){
    ficheros[6].open("./rdos/rdosCmediaA.csv");
    ficheros[7].open("./rdos/rdosEmediaA.csv");
    if((!ficheros[6].is_open())||(!ficheros[7].is_open())){
      printf("Error de creacion de ficheros!!!\n");
      exit(-1);
    }
    for (int i=1;i<=9;i++){  ficheros[6]<<"N"<<i<<" "; }
    ficheros[6]<<std::endl;
    for (int i=1;i<=14;i++){ ficheros[7]<<"N"<<i<<" "; }
    ficheros[7]<<std::endl;
  }
  if(opcion==5 || opcion==8){
    ficheros[8].open("./rdos/rdosCmediaA2.csv");
    ficheros[9].open("./rdos/rdosEmediaA2.csv");
    if((!ficheros[8].is_open())||(!ficheros[9].is_open())){
      printf("Error de creacion de ficheros!!!\n");
      exit(-1);
    }
    for (int i=1;i<=9;i++){  ficheros[8]<<"N"<<i<<" "; }
    ficheros[8]<<std::endl;
    for (int i=1;i<=14;i++){ ficheros[9]<<"N"<<i<<" "; }
    ficheros[9]<<std::endl;
  }
  if(opcion==6 || opcion==8){
    ficheros[10].open("./rdos/rdosCmediaG.csv");
    ficheros[11].open("./rdos/rdosEmediaG.csv");
    if((!ficheros[10].is_open())||(!ficheros[11].is_open())){
      printf("Error de creacion de ficheros!!!\n");
      exit(-1);
    }
    for (int i=1;i<=9;i++){  ficheros[10]<<"N"<<i<<" "; }
    ficheros[10]<<std::endl;
    for (int i=1;i<=14;i++){ ficheros[11]<<"N"<<i<<" "; }
    ficheros[11]<<std::endl;
  }
  if(opcion==7 || opcion==8){
    ficheros[12].open("./rdos/rdosCmediaG2.csv");
    ficheros[13].open("./rdos/rdosEmediaG2.csv");
    if((!ficheros[12].is_open())||(!ficheros[13].is_open())){
      printf("Error de creacion de ficheros!!!\n");
      exit(-1);
    }
    for (int i=1;i<=9;i++){  ficheros[12]<<"N"<<i<<" "; }
    ficheros[12]<<std::endl;
    for (int i=1;i<=14;i++){ ficheros[13]<<"N"<<i<<" "; }
    ficheros[13]<<std::endl;
  }
}

void CerrarFicheros(int opcion,std::ofstream *ficheros){
  //close all outputfiles that are open
  for(int i=0;i<NFS;i++){
    if(ficheros[i].is_open()){
      ficheros[i].close();
      printf("Cerrando el fichero %d\n",i);
    }
  }
}

int main() {

  int opcion=Menu();

  //Abrimos todos los ficheros de salida para ir metiendo los datos
  std::ofstream ficheros[NFS];
  AbrirFicheros(opcion,ficheros);

  //Vamos abriendo conjunto a conjunto y cargando la triangulacion
  for(int fich=0;fich<NF;fich++){
    Delaunay dt=CargarPuntos(fich);

    //Punto central y extremo
    Point2 observa[2];
    observa[0]=Point2(0,0);
    observa[1]=BuscarExtremo(dt);

    //For each observed point, compute the corresponding values
    for(int i=0;i<2;i++){
      //Locate the point in the Delaunay Triangulation
      Point_coordinate_vector coords;
      CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>,K::FT, bool> 
        totalArea = CGAL::natural_neighbor_coordinates_2(dt, observa[i], std::back_inserter(coords));
      if(!totalArea.third){
        std::cout << "The coordinate computation was not successful." << std::endl;
        std::cout << "The point (" << observa[i] << ") lies outside the convex hull." << std::endl;
        return 1;
      }
      //Store distance-1 and distance-2 for latter use
      //Store weights sibson, distance-1 and distance-2 for latter use
      std::vector<flotante> pesoSibson;
      std::vector<flotante> pesoDistancia1;
      std::vector<flotante> pesoDistancia2;

      std::vector<flotante> distancia1;
      std::vector<flotante> distancia2;
      flotante acumuladorDistancia1=0.0;
      flotante acumuladorDistancia2=0.0;


// *sibson pondera peso de cada vecino segun área que le roba = sibN (valor normalizado 0-1)
// *distancia pondera con 1/distancia = idistN (1o calcula y normaliza 0-1)
  // we will store these distances in a vector and accumulate its total as we proceed
// *distancia2 pondera con 1/distancia^2 = idist2N (1o calcula y normaliza 0-1)
  // we will store these distances in a vector and accumulate its total as we proceed
// *mediaA pondera con sibN+idistN/2 (media de los valores ya normalizados y se re-normaliza)
// *mediaG pondera con sqrt(sibN*idistN) (partiendo de valores ya normalizados y se re-normaliza)



      for(int c=0;c<coords.size();c++){
        distancia1.push_back(1/sqrt(squared_distance(coords[c].first,observa[i])));
        acumuladorDistancia1+=distancia1[c];
        distancia2.push_back(1/squared_distance(coords[c].first,observa[i]));
        acumuladorDistancia2+=distancia2[c];

	//bpr debug
	if(squared_distance(coords[c].first,observa[i])<1){
	  std::cout<<squared_distance(coords[c].first,observa[i])<<std::endl;
	  std::cout<<"Punto estudiado con dist al cuadrado menor que 1: "<<observa[i]<<std::endl;
	  std::cout<<"Vecinos: "<<std::endl;
	  for(int k=0;k<coords.size();k++){
	    std::cout<<coords[k].first<<std::endl;
	  }
	}
	//bpr

      }
      //We iterate now on the neighbours and divide by total to obtain relative weight
      for(int c=0;c<coords.size();c++){
        pesoSibson.push_back(coords[c].second/totalArea.second);
        pesoDistancia1.push_back(distancia1[c]/acumuladorDistancia1);
        pesoDistancia2.push_back(distancia2[c]/acumuladorDistancia2);
      }

      //Now that we have normalized values for pesoSibson, pesoDistancia1 and pesoDistancia2,
      //we can compute average and geometric average for each combination
      std::vector<flotante> mediaA;
      std::vector<flotante> mediaA2;
      flotante acumuladorMediaA=0.0;
      flotante acumuladorMediaA2=0.0;
      std::vector<flotante> mediaG;
      std::vector<flotante> mediaG2;
      flotante acumuladorMediaG=0.0;
      flotante acumuladorMediaG2=0.0;
      for(int c=0;c<coords.size();c++){
        mediaA.push_back((pesoSibson[c]+pesoDistancia1[c]/2.0));
        acumuladorMediaA+=mediaA[c];
        mediaA2.push_back((pesoSibson[c]+pesoDistancia2[c]/2.0));
        acumuladorMediaA2+=mediaA2[c];
        mediaG.push_back(sqrt(pesoSibson[c]*pesoDistancia1[c]));
        acumuladorMediaG+=mediaG[c];
        mediaG2.push_back(sqrt(pesoSibson[c]*pesoDistancia2[c]));
        acumuladorMediaG2+=mediaG2[c];
      }
      //With stored averages, we normalize again
      for (int c=0;c<coords.size();c++){
        mediaA[c]=mediaA[c]/acumuladorMediaA;
        mediaA2[c]=mediaA2[c]/acumuladorMediaA2;
        mediaG[c]=mediaG[c]/acumuladorMediaG;
        mediaG2[c]=mediaG2[c]/acumuladorMediaG2;
      } 

      //Move all these data to the files

      if((opcion==1)||(opcion==8)){
        int numfich=0;
        for(int c=0;c<coords.size();c++){
          ficheros[2*numfich+i]<<pesoSibson[c]<<" ";
        }
        ficheros[2*numfich+i]<<"\n";
      }
     
      if((opcion==2)||(opcion==8)){
        int numfich=1;
        for(int c=0;c<coords.size();c++){
          ficheros[2*numfich+i]<<pesoDistancia1[c]<<" ";
        }
        ficheros[2*numfich+i]<<"\n";
      }
      
      if((opcion==3)||(opcion==8)){
        int numfich=2;
        for(int c=0;c<coords.size();c++){
          ficheros[2*numfich+i]<<pesoDistancia2[c]<<" ";
        }
        ficheros[2*numfich+i]<<"\n";
      }
       
      if((opcion==4)||(opcion==8)){
        int numfich=3;
        for(int c=0;c<coords.size();c++){
          ficheros[2*numfich+i]<<mediaA[c]<<" ";
        }
        ficheros[2*numfich+i]<<"\n";
      }
        
      if((opcion==5)||(opcion==8)){
        int numfich=4;
        for(int c=0;c<coords.size();c++){
          ficheros[2*numfich+i]<<mediaA2[c]<<" ";
        }
        ficheros[2*numfich+i]<<"\n";
      }
 
      if((opcion==6)||(opcion==8)){
        int numfich=5;
        for(int c=0;c<coords.size();c++){
          ficheros[2*numfich+i]<<mediaG[c]<<" ";
        }
        ficheros[2*numfich+i]<<"\n";
      }
  
      if((opcion==7)||(opcion==8)){
        int numfich=6;
        for(int c=0;c<coords.size();c++){
          ficheros[2*numfich+i]<<mediaG2[c]<<" ";
        }
        ficheros[2*numfich+i]<<"\n";
      }

    }

  }

  CerrarFicheros(opcion,ficheros);

}

