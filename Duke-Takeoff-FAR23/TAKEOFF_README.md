# Duke Takeoff FAR 23

Calculo preliminar y trazable de la distancia de despegue del Beechcraft Duke
60 mediante el procedimiento de la clase *Performance - Despegue y
Aterrizaje*, complementado por Roskam, capitulo 10.

## Ejecucion

```bash
cd /home/ignacio/Work/TPI_Mecanica_del_Vuelo/Duke-Takeoff-FAR23
python takeoff_duke.py
```

Comprobaciones del caso de verificacion, separado de la interfaz editable:

```bash
python -m unittest -v
```

El primer bloque de `takeoff_duke.py` contiene las condiciones del aeropuerto;
el segundo, la condicion del avion; y el tercero, las opciones del modelo. No
hay dependencias externas: se usa solamente la biblioteca estandar de Python.

## Metodo

La nomenclatura es la de la clase:

```text
s_TO = s_G + s_A = s_NGR + s_R + s_TR + s_CL
```

- `s_NGR`: ecuacion longitudinal de la clase, Ec. (4), e integracion por punto
  medio de la Ec. (8). La velocidad de aire comienza en la velocidad de viento
  de frente, porque el avion parte con velocidad de suelo nula. El metodo
  alternativo evalua la aceleracion en `0.74 V_R`, Ecs. (9)-(10), y solo admite
  viento nulo tal como en las diapositivas.
- `s_R`: velocidad media entre `V_R` y `V_LOF`, tiempo de rotacion de 1 s para
  avion liviano y termino de viento escrito como en la Ec. (11).
- `s_TR`: arco circular, `Delta CL` estadistico, radio `R_TR` y geometria de las
  Ecs. (12)-(19).
- `s_CL`: altura restante dividida por `tan(gamma_c)`, Ec. (20). Si los 50 ft
  se alcanzan dentro de la transicion, `s_CL=0` y se usa la Ec. (22) sin la
  aproximacion de angulo pequeno.
- FAR 23: para un avion normal multiengine se calcula hasta el obstaculo de
  50 ft y se exige `V_50 >= 1.2 V_S1` y `V_50 >= 1.1 V_MC` (Roskam, tabla
  10.1). No se aplica el factor 1.15 de FAR 25; por eso el factor sobre la
  distancia fisica es 1.0.

Las estimaciones preliminares `V_R=1.10 V_S` y `V_LOF=1.15 V_S` proceden de
Roskam 10.3.3.1. `V_50=1.20 V_S` procede de su tabla FAR 23. El codigo no es
una herramienta de certificacion ni reemplaza el AFM.

## Datos y limites conocidos

Los valores `M`, `S`, `CL_max`, `CD0` y `K` corresponden a `Duke.py`. La curva
local de potencia usa los mismos ajustes cuadraticos por MAP de
`motor_tio541_2900rpm.py` y es valida entre 0 y 22 kft. La conversion a empuje
usa el mismo modelo Solies que `eficiencia_helice.py`, incluido el limite de
empuje estatico.

Dos entradas siguen pendientes de datos Duke especificos:

- `CL_g` y `CD_g` con configuracion de despegue y efecto suelo. El caso de
  verificacion usa `CL_g=0.40`, tomado como valor representativo del ejemplo de
  clase, y evalua `CD_g` con la polar preliminar del repositorio.
- `V_MC` certificada para completar los dos chequeos de velocidad FAR 23.

`Duke.py` contiene tanto `CL_max=1.27` para el avion como `Clmax=1.70` en el
bloque de planta alar equivalente. Se usa el primero y NO se interpreta el
segundo como `CL_max,TO` del avion completo.

## Convenciones

- SI internamente: m, s, kg, N, W y Pa.
- Viento de frente positivo y viento de cola negativo.
- Pendiente positiva cuesta arriba, por lo tanto desfavorable.
- Si no se introduce presion, se calcula la presion ISA a la elevacion y la
  densidad con la temperatura real indicada. Si se introduce una presion, su
  altitud de presion equivalente alimenta la curva de potencia del motor.
