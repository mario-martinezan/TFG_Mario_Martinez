# Adaptar el modelo geometrico a elipsoides

El modelo geometroco que se está usando del articulo de XX lo que hace es crear un teselado hexagonal sobre la superfiece de la esfera teniendo en cuenta la separación angulare entre los centros de cada hexagono. 
Esta es la distribución de puntos de manera equidistante mas optima segun el tames problem y eules (explicar mas).
Para adaptar este tesalado a un elipsoide se me ocurre:

- Cambiar el teselado. Si queremos que el tesalado se adapte a la forma no simetria del elipside, seria conveniente cmabiar los exagonos por otra forma. Esto es un problema matematico que no voy a intentar pq no soy tan listo, además de que dejaria hueco a colocar más elipsoides, cosa q no queremos.
- Mantener el teselado. Manteniendo el teselado hexagonal, aumentaremos los "huecos" entre elipsoides, pero podremos resolverlo de la misma manera que con esferas

Manteniendo el teselado hay varios cambios que se pueden hacer para adaptar el código con el modelo esferico y al calculo de los momentos de inercia y sus correspondinetes radios.
  - Cambiar el caclulo de los momentos de inercia. Aqui tengo dudas. Lo que he hecho aqui es que en cada frame de tiempo, se calculan los momentos de inercia principles. Se asigna el mayor al eje X y se calcula el radio mayor. Luego se calcula el otro radio usando el momento de inercia mas pequeño. Así obtenemos el radio mayor y el radio menor del elipsoide.
  - La separcion entre centros de hexagonos será de dos veces el radio mayor del elipsoide
