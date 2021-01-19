# PSTH_online

## Versión 2 (javi)
## 19/01/2020

Este programa genera PSTH y otros gráficos para los *experimentos de presentación de estímulos auditivos*. Puede hacer estos gráficos independientemente de la cantidad de estimulos. Para que corra solo hace falta hacer correr *hace_todo.m* . 

Esta estructurado de tal manera que son varios scripts corridos en serie que necesitan uno del otro. Esta forma de organizar código no es recomendable. A futuro lo mejor es hacer funciones y un solo script. 

### Para correr se necestia estos archivos en el mismo directorio:

- Archivos generandos INTAN en formato 'One File per Signal Type'

- Archivo .txt con el orden de los estímulos presentados en formato vector columna. La primer columna es un numero id del estimulo y la segunda es el nombre del archivo de audio separados por una tabulación. En el nombre del archivo debe figurar la palabra 'estimulos'. 

- Los mismos archivos de audio de los estímulos presentados
