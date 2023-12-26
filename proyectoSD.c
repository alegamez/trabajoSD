#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <string.h> //para el strtok
#include <stdlib.h> //para el atoi
#include <time.h>

#define NMEDIDAS 24 //Cantidad de medidas o dimensiones de los datos en el problema
#define MASTERPID 0
void hazProblema(int splitSize, int rest, int nMedidas, int nDias, int pid, int prn, char *archivo);
float calcularDistancia(float dia1[], float dia2[], int nMedidas);
void copiarDias(float destino[], float origen[], int nMedidas);
void sustituirPrimero(float diaMejor[][NMEDIDAS], float diaActual[], float mejorDistancia[], float distancia, int nDias);
void sustituirSegundo(float diaMejor[][NMEDIDAS], float diaActual[], float mejorDistancia[], float distancia, int nDias);
void calcularMedia(float media[], float dia1[], float dia2[], int nMedidas);
void imprimirVector(float v[], int size);
float calcularMAPE(float real[], float prediccion[], int h, int pid);

int main(int argc, char *argv[])
{

    int prn; // Numero de procesos
    int pid; // Identificador del proceso
    int splitSize, rest; //Tamaño del las divisiones y su resto
    int nDias; //Numero total de das en el conjunto de datos. Se obtiene leyendo el primer valor del archivo y se usa para determinar la cantidad de dias en el conjunto de datos
    int nMedidas; // Cantidad de medidas o dimensiones de los datos. Se obtiene leyendo el segundo valor del archivo y se usa para dar tamaño a los vectores y realizar operaciones con las medidas
    int salto; //Cantidad de caracteres que hay que saltarse al leer el archivo. Se usa para ajustarse al formato del archivo en operaciones de lectura
    int j;  //Indice en bucles y operaciones de lectura y procesamiento de datos
    char *buff = (char *)malloc(1000 * sizeof(char));  // asignamos memoria dinamica al buffer
    MPI_File fh; // Manipulador de archivo mpi para manejar operaciones de entrada y salida en archivos distribuidos entre multiples procesos
    char archivo[] = "datos_1x.txt";
    char *datos;
    MPI_Offset offset;  // Almacena el desplazamiento dentro del archivo. Se usa para operaciones de lectura o escritura en el archivo para indicar la posicion desde la que hay que leer o escribir los datos
    float mejoresDistancias[2]; // Almacenar las mejores distancias calculadas. Cada elemento representa una de las dos mejores distancias encontradas en algun contexto concreto del programa
    // Mide el tiempo de ejecucion del programa
    double t1;  //Tiempo de inicio de la medicion
    double t2; // Tiempo de detecion de la medicion
    double tiempo = 0; // Tiempo total de ejecucion (t2-t1)
    float media[splitSize][NMEDIDAS]; // Media generada para cada prediccion en el codigo
    float error[splitSize]; // Almacenar los errores asociados a cada prediccion. Hay un error apra cada conjunto de datos
    float MAPE_global; // Almacena el error error global del conjunto de datos completo. Evalua la precision de un modelo de prediccion

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    MPI_Comm_size(MPI_COMM_WORLD, &prn);

    // Lee la primera línea del fichero que tiene: número de días y número de datos/día
    MPI_File_open(MPI_COMM_WORLD, archivo, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read(fh, buff, num, MPI_CHAR, NULL);
    //fh --> manipudaro del archivo, archivo desde el que leo
    //buff --> Donde almaceno los datos leidos
    //num --> numero de elementos que deben leerse del archivo
    //MPI_CHAR --> Tipo de dato que se lee
    //NULL --> Vista completa dle archivo, sin restricciones

    datos = strtok(buff, " "); //Divide en una cadena en partes mas pequeñas basadas en el delimitador
    // Guardo en datos la primera entrada del archivo antes de el primer espacio en blanco
    nDias = atoi(datos); // Convierte en un entero el valor anterior

    datos = strtok(NULL, " "); // Indica que se debe continuar desde la ultima posicion, por lo tanto sera la seguna entrada
    nMedidas = atoi(datos); // Convierto en un entero el valor anterior

    // Calcula numero de dias para cada proceso
    splitSize = 1001 / prn;
    rest = 1001 % prn;

    printf("[%d] Dias: %d, Horas: %d\n", pid, nDias, nMedidas);
    MPI_File_close(&fh); // Cierro el archivo porque ya he leido todo su contenido

    // for (int i = 0; i < 100; i++)
    //   {
    MPI_Barrier(MPI_COMM_WORLD); // Sincronizo la ejecucion de todos los procesos del comunicador
    if (pid == MASTERPID){
        t1 = clock(); //Guardo el instante en el que se comienza a ejecutar
    }
    //Ejecuto el programa
    hazProblema(splitSize, rest, nMedidas, nDias, pid, prn, archivo);
    //Finalizo la sincronizacion de los procesos
    MPI_Barrier(MPI_COMM_WORLD);
    if (pid == MASTERPID){
        t2 = clock(); // Guardo el instante en el que finaliza el programa
        tiempo = (double)(t2 - t1) / CLOCKS_PER_SEC; // Tiempo de ejecucion
        // (t2-t1) diferencia entre tiempo de fin y tiempo de inicio
        // CLOCKS_PER_SEC constante que representa el numero de ciclos de reloj por segundo del ssitema
        //(t2-t1)/CLOCKS_PER_SEC --> Convierte el tiempo de ciclos de reloj a segundos
        printf("Tiempo: %f\n", tiempo);
    }

    // Cada proceso MPI guarda sus resultados localmente
    // (asumimos que cada proceso contribuye a los tres ficheros de salida)
    FILE *prediccionesFile; 
    FILE *mapeFile;
    FILE *tiempoFile;

    // El proceso maestro abre los tres ficheros de salida en modo de escritura
    if (pid == MASTERPID){
        prediccionesFile = fopen("Predicciones.txt", "w");
        mapeFile = fopen("MAPE.txt", "w");
        tiempoFile = fopen("Tiempo.txt", "w");
    }else{
        //Cualqueir otro proceso abre el archivo pero solo si ya existe, añadiendo el nuevo contenido al final
        prediccionesFile = fopen("Predicciones.txt", "a");
        mapeFile = fopen("MAPE.txt", "a");
        tiempoFile = fopen("Tiempo.txt", "a");
    }

    // Llamada a la función hazProblema que ahora devuelve media y error
    hazProblema(splitSize, rest, nMedidas, nDias, pid, prn, archivo);

    // Cada proceso MPI escribe sus resultados en los ficheros
    for (int i = 0; i < splitSize; i++){
        // Escribir predicciones en Predicciones.txt
        for (int j = 0; j < nMedidas; j++){
            fprintf(prediccionesFile, "%.1f ", media[i][j]);
        }
        fprintf(prediccionesFile, "\n"); // Despues de escrbir todos los valores de una fila, añado un salto de linea

        // Escribir MAPE en MAPE.txt
        fprintf(mapeFile, "%.1f\n", error[i]); // Escribo el contenido de los errores en mape.txt
    }

    // El proceso 0 escribe el tiempo y otros detalles en Tiempo.txt
    if (pid == MASTERPID){
        fprintf(tiempoFile, "Tiempo de ejecución: %f segundos\n", tiempo);
        fprintf(tiempoFile, "Archivo procesado: %s\n", archivo);
        fprintf(tiempoFile, "MAPE del conjunto de datos completo: %.1f\n", MAPE_global); // Calcula esto según tu lógica
        fprintf(tiempoFile, "Número de procesos MPI utilizados: %d\n", prn);
    }

    // Cerrar los ficheros
    fclose(prediccionesFile);
    fclose(mapeFile);
    fclose(tiempoFile);

    MPI_Finalize();
}

void hazProblema(int splitSize, int rest, int nMedidas, int nDias, int pid, int prn, char *archivo){
    // Almacena los 1001 valores a los que se les busca el knn
    float datos[splitSize][nMedidas];
    float resto[rest][nMedidas];
    // Dia que se esta procesando de todo el archivo
    float diaActual[nMedidas];
    // Distancia euclidea calculada
    float distancia;
    // Mejores distancias para cada medida de las 1001
    float mejoresDistancias[splitSize][2];
    float mejoresDistanciasResto[rest][2];
    // Datos asociados a las mejores distancias de las 1001
    float diasMejores[splitSize][2][nMedidas];
    float diasMejoresResto[rest][2][nMedidas];

    // Media generada para cada prediccion
    float media[splitSize][nMedidas];
    float error[splitSize];

    // Buffer usado
    char *buff;
    char *aux;
    int i, j, h, salto, tam;
    MPI_File fh;

    MPI_File_open(MPI_COMM_WORLD, archivo, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    // Lee los 1001 dias finales
    h = 0;
    // #pragma omp parallel for private(salto, tam, tid, aux, j) firstprivate(h)
    for (i = pid * splitSize; i < pid * splitSize + splitSize; i++){
        if (i < 999){
            salto = -(192 * 2 + 193 * (1001 - i - 2));
            tam = 191;
        }else if (i == 999){
            tam = 192;
            salto = -(192 * (1001 - i));
        }else{
            tam = 191;
            salto = -(192 * (1001 - i) - 1);
        }
        MPI_File_seek(fh, salto, MPI_SEEK_END);
        MPI_File_read(fh, buff, tam, MPI_CHAR, NULL);
        j = 0;
        aux = strtok(buff, ",");
        while (aux != NULL){
            datos[h][j] = atof(aux);
            aux = strtok(NULL, ",");
            j++;
        }
        h++;
        // printf("\n");
    }

    // 1. Hace que la vista del fichero se salte la cabecera
    MPI_File_set_view(fh, num, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
    // 2. lectura de 192 caracteres hasta la línea nDias - 1001
    // i representa el numero del dia que se lee (desde 0 hasta el dia anterior al correspondiente)

    for (i = 0; i < nDias - (splitSize * prn - pid * splitSize + rest); i++){
        // Se lee cada dia y se hace la distancia euclidea con los dias leidos
        MPI_File_read(fh, buff, 191, MPI_CHAR, NULL);
        // j controla el vector del dia que se esta leyendo
        j = 0;
        aux = strtok(buff, ",");
        while (aux != NULL){
            diaActual[j] = atof(aux);
            aux = strtok(NULL, ",");
            j++;
        }

        // Se calculan las distancias euclideas para cada dia de la matriz datos (no incluye las ultimas 1000 filas)
        // m es cada una de las medidas del dia
        // #pragma omp parallel for
        for (int m = 0; m < splitSize; m++){
            distancia = calcularDistancia(datos[m], diaActual, nMedidas);

            if (i == 0){
                // Inicializamos la primera y segunda mejores medidas al primer valor
                copiarDias(diasMejores[m][0], diaActual, nMedidas);
                mejoresDistancias[m][0] = distancia;
                copiarDias(diasMejores[m][1], diaActual, nMedidas);
                mejoresDistancias[m][1] = distancia;
            }else{

                // printf("\nDistancia calculada: %f < %f < %f ?\n", distancia, mejoresDistancias[m][0], mejoresDistancias[m][1]);
                if (distancia < mejoresDistancias[m][0]){
                    sustituirPrimero(diasMejores[m], diaActual, mejoresDistancias[m], distancia, nMedidas);
                }else if (distancia < mejoresDistancias[m][1]){
                    sustituirSegundo(diasMejores[m], diaActual, mejoresDistancias[m], distancia, nMedidas);
                }
            }
        }
    }

    // TODO: Hacer las cuentas para el resto

    for (int i = 0; i < splitSize; i++){
        calcularMedia(media[i], diasMejores[i][0], diasMejores[i][1], nMedidas);
        //   printf("Mejor Distancia de %d: %.1f\n", pid * splitSize + i, mejoresDistancias[i][0]);
    }

    for (int i = 0; i < splitSize; i++){
        error[i] = calcularMAPE(media[i], datos[i], nMedidas, pid);
        // printf("[%d] Error de %d: %.1f\n", pid, pid * splitSize + i, error[i]);
    }

    MPI_File_close(&fh);
}

float calcularDistancia(float dia1[], float dia2[], int nMedidas){
     float result = 0;
    int i;

    // #pragma omp parallel for private(i) reduction(+:result)
    for (i = 0; i < nMedidas; i++){
        result += pow(dia1[i] - dia2[i], 2);
    }
    return sqrt(result);
}

void sustituirPrimero(float diaMejor[][NMEDIDAS], float diaActual[], float mejorDistancia[], float distancia, int nDias){
    sustituirSegundo(diaMejor, diaMejor[0], mejorDistancia, mejorDistancia[0], nDias);
    mejorDistancia[0] = distancia;
    copiarDias(diaMejor[0], diaActual, nDias);
}

void sustituirSegundo(float diaMejor[][NMEDIDAS], float diaActual[], float mejorDistancia[], float distancia, int nDias){
    mejorDistancia[1] = distancia;
    copiarDias(diaMejor[1], diaActual, nDias);
}

void calcularMedia(float media[], float dia1[], float dia2[], int nMedidas){
    for (int i = 0; i < nMedidas; i++){
        media[i] = (dia1[i] + dia2[i]) / 2;
    }
}

void copiarDias(float destino[], float origen[], int nMedidas){
    for (int i = 0; i < nMedidas; i++){
        destino[i] = origen[i];
    }
}

void imprimirVector(float v[], int size){
    for (int i = 0; i < size; i++){
        printf("%f ", v[i]);
    }
    printf("\n");
}

float calcularMAPE(float real[], float prediccion[], int h, int pid){
    float res = 0;

    for (int i = 0; i < h; i++){
        res += (fabs((real[i] - prediccion[i]) / real[i]));
    }
    res = 100 * res / h;

    return res;
}