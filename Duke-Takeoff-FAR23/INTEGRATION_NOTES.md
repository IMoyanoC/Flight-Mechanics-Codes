# Notas de integracion futura

El repositorio `/home/ignacio/Work/TPI_Mecanica_del_Vuelo/Flight-Mechanics-Codes`
fue inspeccionado exclusivamente en lectura. No se realizo la integracion.

## Correspondencias

| Necesidad | Modulo existente | Interfaz / variables |
|---|---|---|
| Datos Duke | `Duke.py` | `M`, `S`, `CL_max`, `CD0`, `K` |
| Atmosfera ISA | `Climb.py` | `atmosfera_isa(h) -> rho, a` |
| Densidad ISA simple | `eficiencia_helice.py` | `densidad_isa(h)` |
| Potencia de motor | `motor_tio541_2900rpm.py` | `potencia_total_duke(altitud, MAP, unidad_altitud="m")` |
| Eficiencia de helice | `eficiencia_helice.py` | `eficiencia_helice(V, h)` |
| Polar | `Climb.py` / `Duke.py` | `coeficiente_CD(CL)` / `CD0 + K CL^2` |
| Polar de baja velocidad | `Aerodynamics.py` | `polar_LS_func(CL)`; requiere depuracion antes de reutilizar |

`takeoff_model.shaft_power_available` puede sustituirse directamente por una
envoltura de `potencia_total_duke`, convirtiendo hp a W. Hay que pasar
`unidad_altitud="m"` de forma explicita: el valor por defecto de
`potencia_total_duke` es `"ft"`.

`takeoff_model.propeller_efficiency_and_thrust` puede sustituirse por
`eficiencia_helice(V, h)` y `T=eta P_eje/V`, conservando aparte el limite
estatico. Antes de hacerlo debe resolverse si la densidad de la helice se desea
ISA o la densidad real del aeropuerto.

## Hallazgos que requieren decision

- `Duke.py` declara `dp=2.286 m`, mientras `propeller.py` y
  `eficiencia_helice.py` usan `D=1.8796 m`. El proyecto nuevo conserva 1.8796 m
  para ser coherente con el modelo Solies existente; debe confirmarse el
  diametro fisico antes de integrar.
- `Duke.py` declara `CL_max=1.27` y tambien `Clmax=1.70` para la planta alar
  equivalente. No son intercambiables; se eligio `CL_max=1.27` como parametro
  del avion completo.
- No se encontraron `CL_g`, `CD_g` corregidos por efecto suelo ni `V_MC`
  certificada del Duke. Esos datos siguen como `TODO_USER_INPUT`.
- `Climb.py` ya combina atmosfera, potencia, eficiencia y polar, pero es un
  script con calculos y graficos al importar; no conviene importarlo como
  biblioteca sin una refactorizacion futura.
- Los modulos reutilizables dependen de NumPy y el de motor puede usar SciPy.
  El entorno disponible para esta entrega no los tiene instalados; por ello se
  implementaron interfaces locales minimas, sin dependencias externas, con las
  mismas ecuaciones/datos.

## Integracion sugerida

Una integracion posterior deberia inyectar tres funciones en el modelo de
despegue: `atmosfera(h, T, p)`, `potencia_eje_disponible(h, MAP)` y
`empuje_disponible(V, h)`. Esto evita importar los scripts que generan graficos
y permite comparar la implementacion local contra los modulos existentes sin
cambiar las ecuaciones de las cuatro fases.
