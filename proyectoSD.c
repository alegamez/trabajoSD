#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <string.h> //para el strtok
#include <stdlib.h> //para el atoi
#include <time.h>

#define k 2
#define NMEDIDAS 24
#define num 8

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

    int prn, pid, splitSize, rest, value, nDias, nMedidas, salto, j;
    char *buff = (char *)malloc(1000 * sizeof(char));  // asignamos memoria
    MPI_File fh;
    char archivo[] = "datos_1x.txt";
    char *datos;
    MPI_Offset offset;
    float mejoresDistancias[2];
    double t1, t2, tiempo = 0;
    float media[splitSize][NMEDIDAS];
    float error[splitSize];
    float MAPE_global;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    MPI_Comm_size(MPI_COMM_WORLD, &prn);

    // Lee la primera línea del fichero que tiene: número de días y número de datos/día
    MPI_File_open(MPI_COMM_WORLD, archivo, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read(fh, buff, num, MPI_CHAR, NULL);

    datos = strtok(buff, " ");
    nDias = atoi(datos);

    datos = strtok(NULL, " ");
    nMedidas = atoi(datos);

    // Calcula numero de dias para cada proceso
    splitSize = 1001 / prn;
    rest = 1001 % prn;

    printf("[%d] Dias: %d, Horas: %d\n", pid, nDias, nMedidas);
    MPI_File_close(&fh);

    // for (int i = 0; i < 100; i++)
    //   {
    MPI_Barrier(MPI_COMM_WORLD);
    if (pid == 0)
    {
        t1 = clock();
    }
    hazProblema(splitSize, rest, nMedidas, nDias, pid, prn, archivo);
    MPI_Barrier(MPI_COMM_WORLD);
    if (pid == 0)
    {
        t2 = clock();
        tiempo = (double)(t2 - t1) / CLOCKS_PER_SEC;
        printf("Tiempo: %f\n", tiempo);
    }

    // Cada proceso MPI guarda sus resultados localmente
    // (asumimos que cada proceso contribuye a los tres ficheros de salida)
    FILE *prediccionesFile;
    FILE *mapeFile;
    FILE *tiempoFile;

    if (pid == 0)
    {
        prediccionesFile = fopen("Predicciones.txt", "w");
        mapeFile = fopen("MAPE.txt", "w");
        tiempoFile = fopen("Tiempo.txt", "w");
    }
    else
    {
        prediccionesFile = fopen("Predicciones.txt", "a");
        mapeFile = fopen("MAPE.txt", "a");
        tiempoFile = fopen("Tiempo.txt", "a");
    }

    // Llamada a la función hazProblema que ahora devuelve media y error
    hazProblema(splitSize, rest, nMedidas, nDias, pid, prn, archivo);

    // Cada proceso MPI escribe sus resultados en los ficheros
    for (int i = 0; i < splitSize; i++)
    {
        // Escribir predicciones en Predicciones.txt
        for (int j = 0; j < nMedidas; j++)
        {
            fprintf(prediccionesFile, "%.1f ", media[i][j]);
        }
        fprintf(prediccionesFile, "\n");

        // Escribir MAPE en MAPE.txt
        fprintf(mapeFile, "%.1f\n", error[i]);
    }

    // El proceso 0 escribe el tiempo y otros detalles en Tiempo.txt
    if (pid == 0)
    {
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

void hazProblema(int splitSize, int rest, int nMedidas, int nDias, int pid, int prn, char *archivo)
{
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
    for (i = pid * splitSize; i < pid * splitSize + splitSize; i++)
    {
        if (i < 999)
        {
            salto = -(192 * 2 + 193 * (1001 - i - 2));
            tam = 191;
        }
        else if (i == 999)
        {
            tam = 192;
            salto = -(192 * (1001 - i));
        }
        else
        {
            tam = 191;
            salto = -(192 * (1001 - i) - 1);
        }
        MPI_File_seek(fh, salto, MPI_SEEK_END);
        MPI_File_read(fh, buff, tam, MPI_CHAR, NULL);
        j = 0;
        aux = strtok(buff, ",");
        while (aux != NULL)
        {
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

    for (i = 0; i < nDias - (splitSize * prn - pid * splitSize + rest); i++)
    {
        // Se lee cada dia y se hace la distancia euclidea con los dias leidos
        MPI_File_read(fh, buff, 191, MPI_CHAR, NULL);
        // j controla el vector del dia que se esta leyendo
        j = 0;
        aux = strtok(buff, ",");
        while (aux != NULL)
        {
            diaActual[j] = atof(aux);
            aux = strtok(NULL, ",");
            j++;
        }

        // Se calculan las distancias euclideas para cada dia de la matriz datos (no incluye las ultimas 1000 filas)
        // m es cada una de las medidas del dia
        // #pragma omp parallel for
        for (int m = 0; m < splitSize; m++)
        {
            distancia = calcularDistancia(datos[m], diaActual, nMedidas);

            if (i == 0)
            {
                // Inicializamos la primera y segunda mejores medidas al primer valor
                copiarDias(diasMejores[m][0], diaActual, nMedidas);
                mejoresDistancias[m][0] = distancia;
                copiarDias(diasMejores[m][1], diaActual, nMedidas);
                mejoresDistancias[m][1] = distancia;
            }
            else
            {

                // printf("\nDistancia calculada: %f < %f < %f ?\n", distancia, mejoresDistancias[m][0], mejoresDistancias[m][1]);
                if (distancia < mejoresDistancias[m][0])
                {
                    sustituirPrimero(diasMejores[m], diaActual, mejoresDistancias[m], distancia, nMedidas);
                }
                else if (distancia < mejoresDistancias[m][1])
                {
                    sustituirSegundo(diasMejores[m], diaActual, mejoresDistancias[m], distancia, nMedidas);
                }
            }
        }
    }

    // TODO: Hacer las cuentas para el resto

    for (int i = 0; i < splitSize; i++)
    {
        calcularMedia(media[i], diasMejores[i][0], diasMejores[i][1], nMedidas);
        //   printf("Mejor Distancia de %d: %.1f\n", pid * splitSize + i, mejoresDistancias[i][0]);
    }

    for (int i = 0; i < splitSize; i++)
    {
        error[i] = calcularMAPE(media[i], datos[i], nMedidas, pid);
        // printf("[%d] Error de %d: %.1f\n", pid, pid * splitSize + i, error[i]);
    }

    MPI_File_close(&fh);
}

float calcularDistancia(float dia1[], float dia2[], int nMedidas)
{
     float result = 0;
    int i;

    // #pragma omp parallel for private(i) reduction(+:result)
    for (i = 0; i < nMedidas; i++)
    {
        result += pow(dia1[i] - dia2[i], 2);
    }
    return sqrt(result);
}

void sustituirPrimero(float diaMejor[][NMEDIDAS], float diaActual[], float mejorDistancia[], float distancia, int nDias)
{
    sustituirSegundo(diaMejor, diaMejor[0], mejorDistancia, mejorDistancia[0], nDias);
    mejorDistancia[0] = distancia;
    copiarDias(diaMejor[0], diaActual, nDias);
}

void sustituirSegundo(float diaMejor[][NMEDIDAS], float diaActual[], float mejorDistancia[], float distancia, int nDias)
{
    mejorDistancia[1] = distancia;
    copiarDias(diaMejor[1], diaActual, nDias);
}

void calcularMedia(float media[], float dia1[], float dia2[], int nMedidas)
{
    for (int i = 0; i < nMedidas; i++)
    {
        media[i] = (dia1[i] + dia2[i]) / 2;
    }
}

void copiarDias(float destino[], float origen[], int nMedidas)
{
    for (int i = 0; i < nMedidas; i++)
    {
        destino[i] = origen[i];
    }
}

void imprimirVector(float v[], int size)
{
    for (int i = 0; i < size; i++)
    {
        printf("%f ", v[i]);
    }
    printf("\n");
}

float calcularMAPE(float real[], float prediccion[], int h, int pid)
{
    float res = 0;

    for (int i = 0; i < h; i++)
    {
        res += (fabs((real[i] - prediccion[i]) / real[i]));
    }
    res = 100 * res / h;

    return res;
}