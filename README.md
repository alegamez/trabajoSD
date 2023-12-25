# Proyecto de Predicción de Series Temporales en Sistemas Distribuidos

¡[@alegamez](https://github.com/alegamez)! 👋 ¡Proyecto de sistemas distribuidos para la predicción de series temporales! 🚀

## Introducción

En este proyecto, implementaremos un algoritmo de predicción de series temporales. ¿Qué es una serie temporal? ¡Es como una historia cronológica de datos! 📈 En nuestro caso, estamos buscando predecir el futuro basándonos en patrones del pasado.

### Objetivo

Desarrollar un algoritmo que prediga los valores futuros de una serie temporal utilizando técnicas de aprendizaje automático y evaluar su rendimiento en términos de escalabilidad.

## Descripción del Método

¡Vamos a hacer esto simple y efectivo! 🚀 Nuestra solución se basará en encontrar los k vecinos más cercanos. Para entenderlo mejor, aquí hay un ejemplo rápido:

Supongamos que tenemos mediciones horarias de una serie temporal. Queremos predecir los valores para el próximo día. Entonces, buscaremos en nuestro historial de datos los días más similares al actual y miraremos qué pasó después en esos días. Tomaremos la media de esos valores para realizar la predicción. Simple, ¿verdad?

## Ejemplo Gráfico

¡Mira esta figura que ilustra cómo funciona! 📊

![Proceso de Predicción](link_a_la_imagen.png)

## Estudio de la Escalabilidad

No solo queremos un algoritmo efectivo, ¡queremos que sea rápido sin importar el tamaño de los datos! ⚡ Vamos a probar nuestra solución con conjuntos de datos de diferentes tamaños para ver cómo se desempeña. Estudiaremos el tiempo de procesamiento y ajustaremos el número de procesos e hilos para optimizar el rendimiento.

## Medición de la Calidad de la Predicción

Queremos asegurarnos de que nuestras predicciones sean precisas. Utilizaremos la métrica MAPE (Mean Absolute Percentage Error) para evaluar la calidad de nuestras predicciones. ¡Aim for the stars, pero también apunta a un MAPE bajo! 🌟

### Fórmula MAPE

$$
MAPE = \frac{1}{n} \sum_{i=1}^n \left| \frac{A_i - F_i}{A_i} \right| \times 100
$$


## Contacto

Si tienes alguna pregunta, inquietud o simplemente quieres compartir tus ideas brillantes, ¡no dudes en contactarme en [alejandrogamezct@gmail.com](mailto:alejandrogamezct@gmail.com)!

¡Vamos a hacer de este proyecto algo asombroso! 🚀💻🌟
