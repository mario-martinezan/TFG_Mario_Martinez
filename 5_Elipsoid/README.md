# Como adaptar el modelo geometrico a elipsoides (explicación para enterarme)

El modelo geométrico que se está usando del articulo de M. Soloviev et al. (2021) crea un teselado esferico colocando las esferas pequeñas en el centro de hexagonos sobre la superfiece de la esfera grande teniendo en cuenta la separación angular entre los centros de cada hexagono. 
Este modelo es una solución asintótica del problema de Tammes. Esta disposición minimiza la distancia entre los puntos distribuidos en la esfera y maximiza la densidad de empaquetamiento.
Para adaptar este tesalado a elipsoides en vez de esferas se me ocurre:

- Cambiar el teselado. Si queremos que el tesalado se adapte a la forma no simetria en todas las direcciones del elipside, seria conveniente cambiar los hexagonos por otra forma, ya que manteniendo la forma hexagonal, se crearan "huecos" y los elipsoides no estarán en contancto entre ellos. Esto es un problema de empaquetamiento más complicado, ademas que aumentaria el número de elipsoides que se podrian colocar sobre la superficie, algo que no interesa.
- Mantener el teselado. Manteniendo el teselado hexagonal, aumentaremos los "huecos" entre elipsoides, pero podremos resolverlo de la misma manera que con esferas.

Manteniendo el teselado hay varios cambios que se pueden hacer para adaptar el código con el modelo esferico y al calculo de los momentos de inercia y sus correspondinetes radios.
  - Cambiar el caclulo de los momentos de inercia. Lo que he hecho aqui es que en cada frame de tiempo, se calculan los momentos de inercia principles. Se asigna el mayor de estos tres momentos al eje X y se calcula el radio mayor a partir de $$I_x = \dfrac{2}{5} M b^2$$ donde $b$ será el radio mayor. Luego se calcula el otro radio usando el momento de inercia mas pequeño y la formula $$I_y= \dfrac{1}{5} M (a^2+b^2)$$ donde $a$ es el radio menor. Así obtenemos el radio mayor y el radio menor del elipsoide (Aqui tengo dudas sobre si son correctas estas formulas). Inertia radius at equilibrium (nm) = [  3.9303248056934716  ,  2.712936151092221  ,  1.9330036911570772  ]
  - En relación a la siguiente figura,
    
     <img src="https://github.com/user-attachments/assets/c7164393-5cdc-4b08-b4f5-e4315c7a82b7" width="200px">

    R2 ya no será igual al segmento CB. CB pasa a ser el radio mayor (b) y R2 será el radio emnos (a). Por lo tanto el calculo de $\alpha$ cambia a $\sin\alpha = \dfrac{b}{R1+a}$ (duda: sigue siendo el sinus o es la tangente??? las graficas estan hechas con el sinus). Como se ve en la figura, al mantener en teselado no tiene sentido usar un elipsoide asimetrico en los tres angulos, por lo que usaré el radio principal de inercia más pequeño.

En relación a estos cambios, he hecho dos calculos. Uno que solo hace el segundo cambio (es decir, adaptando a la forma del elipsoide pero calculando los radios con las formulas de los momentos de inercia de una esfera dura) y un con los dos cambios (forma elipsoidal y radios calculados con los momentos de inercia del elipsoide escritos arriba).




