# Framework de simulación cuántica para la factorización de números primos en rust.

## Infraestructura

Archivos principales del proyecto (ordenados de mayor a menor granularidad):

- `main.rs`: funcionalidad completa, abstrae la mayor parte de detalles de implementación.
- `shor.rs`: implementa el algoritmo de shor, abstrae el funcionamiento del simulador cuántico.
- `sim.rs`: implementa el simulador cuántico.
- `utils.rs`: implementa varias funciones usadas a lo largo del proyecto.

## Uso

Para obtener los mejores resultados se recomienda compilar el proyecto en modo `release` (también se puede ejecutar con _cargo_ o compilar con valores por defecto pero la ejecución se verá notablemente penalizada):

```bash
cargo build --release
./target/release/tfm <número a factorizar> <opcional | qbits adicionales de precisión (3 por defecto)>
```

## Interpretación de resultados

Cada iteración del algoritmo imprimirá por linea de commandos los siguientes valores:

- `N`: el número a factorizar.
- `a`: el valor pseudoaleatoriamente generado para la iteracion concreta.
- `t`: el valor de precisión especificado.
- Los siguientes valores de cada registro:
  1. Índice del registro (binario).
  2. Estado del registro (número complejo).
  3. Probabilidad del registro (entre 0 y 1).
  4. Posibles resultados derivados de la subrutina cuántica.
  5. **Opcional**: `[~]`.
  6. **Opcional**: `p=<factor 1>, q=<factor 2> [!]`.
