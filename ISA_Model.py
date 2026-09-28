"""
MODELO INTERNACIONAL DE ARMOSFERA ESTANDAR (ISA)
Se definen los parámetros atmosféricos más importantes, calculados
según lo establecido en el modelo de atmósfera estandar únicamente 
válido entre 0 y 11000 m (Tropósfera).   

- Los valores con subíndice "_0" corresponden al valor que toma la 
variable a nivel del mar.  

-T y P son variables primarias.
    Si no se provee un valor fijo de los mismos, el modelo supone
    el calculado según estándar para el nivel de altitud, si se hace, 
    utiliza el mismo para el cálculo de las variables derivadas

-rho y mu son variables derivadas.

    Altitud:
        unidades="m"  -> h en metros (POR DEFECTO)
        unidades="ft" -> h en pies

    Temperatura:
        Por defecto, valores introducidos en Kelvin.
        Puede utilizarse T_C para introducir grados Celsius.

    Presión:
        Pa
        
MODO DE USO
Devuelve únicamente los parámetros solicitados.
    
    from ISA_Model import ISA
    
            Parámetros posibles:
            "T"   -> temperatura [K]
            "P"   -> presión [Pa]
            "rho" -> densidad [kg/m³]
            "mu"  -> viscosidad dinámica [Pa s]

        T y P pueden sobrescribirse independientemente.

            Ejemplos:
                ISA().tk(3000, "rho")
                ISA().tk(3000, "T", "P")
                ISA().tk(0, "rho", T_C=25, P=101325)           
    
    rho = ISA(unidades="ft").tk(3000.0, "rho")
"""
from math import isfinite

#CONSTANTES


GAMMA_AIR = 1.4
"""
Rho_0 = 1.225      # Kg/m**3

def rho(h):
    rho = Rho_0 * (1-2.2558e-5 * h)**4.2559
    return rho
"""



class ISA:

    # Constantes ISA
    T_0 = 288.15          # K
    P_0 = 101325.0        # Pa
    R = 287.05287        # J/(kg K)
    g_0 = 9.80665         # m/s²
    LAPSE_T = -0.0065          # K/m

    # Sutherland
    T_REF = 273.15       # K
    MU_REF = 1.716e-5    # Pa s
    S = 110.4            # K

    H_MAX = 11000.0     # m

    def __init__(self, unidades="m"):
        if unidades not in ("m", "ft"):
            raise ValueError("unidades debe ser 'm' o 'ft'.")

        self.unidades = unidades

    def _altura_m(self, h):
        """Convierte la altitud recibida a metros."""
        h = float(h)

        if not isfinite(h):
            raise ValueError("La altitud debe ser un número finito.")

        if self.unidades == "ft":
            h *= 0.3048

        if not 0.0 <= h <= self.H_MAX:
            raise ValueError(
                "El modelo está definido únicamente para "
                "la troposfera: 0 <= h <= 11000 m."
            )

        return h

    @classmethod
    def _temperatura_isa(cls, h):
        """Temperatura ISA [K]."""
        return cls.T_0 + cls.LAPSE_T * h

    @classmethod
    def _presion_isa(cls, h):
        """Presión ISA [Pa]."""
        T = cls._temperatura_isa(h)

        return cls.P_0 * (T / cls.T_0) ** (
            -cls.g_0 / (cls.R * cls.LAPSE_T)
        )

    @classmethod
    def _viscosidad(cls, T):
        """Viscosidad dinámica mediante ley de Sutherland [Pa s]."""
        return (
            cls.MU_REF * (T / cls.T_REF) ** 1.5 * (cls.T_REF + cls.S) / (T + cls.S)
        )

    def tk(
        self,
        h,
        *parametros,
        T=None,
        T_C=None,
        P=None,
    ):


        if not parametros:
            raise ValueError(
                "Debe solicitar al menos un parámetro: "
                "'T', 'P', 'rho' o 'mu'."
            )

        validos = {"T", "P", "rho", "mu"}

        for parametro in parametros:
            if parametro not in validos:
                raise ValueError(
                    f"Parámetro desconocido: {parametro}"
                )

        if T is not None and T_C is not None:
            raise ValueError(
                "Especifique T o T_C, no ambos."
            )

        h_m = self._altura_m(h)

        # --------------------------------------------------
        # Variables primarias
        # --------------------------------------------------

        if T_C is not None:
            T_real = float(T_C) + 273.15

        elif T is not None:
            T_real = float(T)

        else:
            T_real = self._temperatura_isa(h_m)

        if P is not None:
            P_real = float(P)

        else:
            P_real = self._presion_isa(h_m)

        if T_real <= 0:
            raise ValueError(
                "La temperatura absoluta debe ser mayor que 0 K."
            )

        if P_real <= 0:
            raise ValueError(
                "La presión debe ser mayor que 0 Pa."
            )

        # --------------------------------------------------
        # Variables derivadas
        # --------------------------------------------------

        rho = P_real / (self.R * T_real)
        mu = self._viscosidad(T_real)

        resultados = {
            "T": T_real,
            "P": P_real,
            "rho": rho,
            "mu": mu,
        }

        # Un único parámetro -> float
        if len(parametros) == 1:
            return resultados[parametros[0]]

        # Varios parámetros -> diccionario
        return {
            parametro: resultados[parametro]
            for parametro in parametros
        }

    def altitude_from_pressure(self, P):
        """Calcula la altitud a partir de la presión [Pa]."""
        P = float(P)

        if not isfinite(P):
            raise ValueError("La presión debe ser un número finito.")

        if P <= 0:
            raise ValueError("La presión debe ser mayor que 0 Pa.")

        h = self.T_0 / self.LAPSE_T * (
            (P / self.P_0) ** (-self.R * self.LAPSE_T / self.g_0) - 1
        )

        if not 0.0 <= h <= self.H_MAX:
            raise ValueError(
                "La presión corresponde a una altitud fuera del rango "
                "definido por el modelo: 0 <= h <= 11000 m."
            )

        if self.unidades == "ft":
            h /= 0.3048

        return h